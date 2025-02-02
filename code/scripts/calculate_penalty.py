import csv
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed

# from cachetools import cached, TTLCache
import redis
import requests
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from Bio.Seq import Seq


# @cached(TTLCache(maxsize=100000, ttl=86400))
# def fetch_cds_sequence(ensembl_transcript_id):
#     server = "https://rest.ensembl.org"
#     ext = f"/sequence/id/{ensembl_transcript_id}?type=cds"
#
#     headers = {"Content-Type": "text/plain"}
#     response = requests.get(server + ext, headers=headers)
#
#     if not response.ok:
#         response.raise_for_status()
#         return None
#
#     return response.text.strip()

def fetch_cds_sequence(ensembl_transcript_id):
    """
    Fetches the CDS (Coding DNA Sequence) for a given Ensembl transcript ID.
    Utilizes Redis caching to prevent redundant API calls.

    Parameters:
        ensembl_transcript_id (str): The Ensembl transcript ID (e.g., "ENST00000507110").

    Returns:
        str: The CDS sequence as a string.
    """
    # Connect to Redis
    redis_conn = redis.Redis(host="localhost", port=6379, db=0, decode_responses=True)

    # Try to get the sequence from Redis cache
    cache_key = f"cds_sequence:{ensembl_transcript_id}"
    cached_sequence = redis_conn.get(cache_key)

    if cached_sequence:
        # Sequence found in cache
        redis_conn.expire(cache_key, 60 * 86400)
        redis_conn.close()
        return cached_sequence
    else:
        # Sequence not in cache; fetch from Ensembl REST API
        server = "https://rest.ensembl.org"
        ext = f"/sequence/id/{ensembl_transcript_id}?type=cds"

        headers = {"Content-Type": "text/plain"}
        response = requests.get(server + ext, headers=headers)

        if not response.ok:
            # response.raise_for_status()
            return None

        cds_sequence = response.text.strip()

        # Store the sequence in Redis cache with an expiration time (optional)
        redis_conn.set(cache_key, cds_sequence, ex=60 * 86400)  # Cache expires in 60*24 hours (60*86400 seconds)

        redis_conn.close()

        return cds_sequence


def translate_cds_to_protein(cds_sequence):
    coding_dna = Seq(cds_sequence)
    protein_sequence = coding_dna.translate(to_stop=True)
    return str(protein_sequence)


def apply_mutation(protein_seq, aa_mutation):
    import re

    # Deletion-Insertion (delins) Mutations with Optional Insertion
    delins_pattern = r"p\.([A-Z\*])(\d+)_([A-Z\*])(\d+)delins([A-Z\*]*)$"
    delins_match = re.match(delins_pattern, aa_mutation)
    if delins_match:
        start_aa = delins_match.group(1)
        start_pos = int(delins_match.group(2))
        end_aa = delins_match.group(3)
        end_pos = int(delins_match.group(4))
        insertion = delins_match.group(5)  # May be empty

        # Apply the deletion and insertion (insertion may be empty)
        mutated_seq = protein_seq[:start_pos - 1] + insertion + protein_seq[end_pos:]
        return mutated_seq

    # Deletion Mutations
    del_pattern = r"p\.([A-Z\*])(\d+)_([A-Z\*])(\d+)del$"
    del_match = re.match(del_pattern, aa_mutation)
    if del_match:
        start_aa = del_match.group(1)
        start_pos = int(del_match.group(2))
        end_aa = del_match.group(3)
        end_pos = int(del_match.group(4))
        # Apply the deletion
        mutated_seq = protein_seq[:start_pos - 1] + protein_seq[end_pos:]
        return mutated_seq

    # Duplication (dup)
    dup_pattern = r"p\.([A-Z])(\d+)_([A-Z])(\d+)dup"
    dup_match = re.match(dup_pattern, aa_mutation)
    if dup_match:
        start_pos = int(dup_match.group(2))
        end_pos = int(dup_match.group(4))
        duplicated_segment = protein_seq[start_pos - 1:end_pos]
        mutated_seq = protein_seq[:end_pos] + duplicated_segment + protein_seq[end_pos:]
        return mutated_seq

    # Frameshift mutations (fs)
    frameshift_pattern = r"p\.([A-Z])(\d+)([A-Z])fs\*(\d+)"
    frameshift_match = re.match(frameshift_pattern, aa_mutation)
    if frameshift_match:
        pos = int(frameshift_match.group(2))
        new_aa = frameshift_match.group(3)
        stop_distance = int(frameshift_match.group(4))
        # Simulate frameshift: Replace from position with new aa and add stop codon after stop_distance
        frameshift_seq = new_aa + "X" * (stop_distance - 2) + "*"
        mutated_seq = protein_seq[:pos - 1] + frameshift_seq
        return mutated_seq

    # Nonsense mutations (*)
    nonsense_pattern = r"p\.([A-Z])(\d+)\*"
    nonsense_match = re.match(nonsense_pattern, aa_mutation)
    if nonsense_match:
        pos = int(nonsense_match.group(2))
        # Introduce a stop codon at the position
        mutated_seq = protein_seq[:pos - 1]
        return mutated_seq

    # Missense mutations
    missense_pattern = r"p\.([A-Z])(\d+)([A-Z])$"
    missense_match = re.match(missense_pattern, aa_mutation)
    if missense_match:
        pos = int(missense_match.group(2))
        mut_aa = missense_match.group(3)
        # Apply the amino acid substitution
        mutated_seq = protein_seq[:pos - 1] + mut_aa + protein_seq[pos:]
        return mutated_seq

    # Synonymous mutations (=)
    synonymous_pattern = r"p\.([A-Z])(\d+)="
    synonymous_match = re.match(synonymous_pattern, aa_mutation)
    if synonymous_match:
        # No change in protein sequence
        return protein_seq

    # Unknown effect (p.?)
    unknown_pattern = r"p\.\?"
    unknown_match = re.match(unknown_pattern, aa_mutation)
    if unknown_match:
        # Cannot apply mutation; return original sequence or raise an exception
        print(f"Warning: Mutation effect is unknown for \"{aa_mutation}\". Returning original sequence.")
        return protein_seq  # or raise an exception if preferred

    # If mutation type is not recognized
    raise ValueError(f"Mutation type not recognized or not supported: {aa_mutation}")


def align_sequences(seq1, seq2):
    # Load BLOSUM62 matrix
    blosum62 = substitution_matrices.load("BLOSUM62")

    # Create an aligner object
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = blosum62
    aligner.open_gap_score = -10  # Gap opening penalty
    aligner.extend_gap_score = -0.5  # Gap extension penalty

    # Perform alignment
    alignments = aligner.align(seq1, seq2)
    best_alignment = alignments[0]
    # aligned_seq1 = str(best_alignment.aligned[0])
    # aligned_seq2 = str(best_alignment.aligned[1])

    # Extract aligned sequences
    aligned_seq1, aligned_seq2 = extract_aligned_sequences(best_alignment)

    alignment_score = best_alignment.score
    return aligned_seq1, aligned_seq2, alignment_score


def extract_aligned_sequences(alignment):
    # Convert the alignment object to aligned sequences
    seqA = []
    seqB = []
    positions = alignment.aligned
    s1 = alignment.sequences[0]
    s2 = alignment.sequences[1]

    idx1 = 0
    idx2 = 0
    for (start1, end1), (start2, end2) in zip(*positions):
        # Add matches/mismatches
        while idx1 < start1 and idx2 < start2:
            seqA.append(s1[idx1])
            seqB.append(s2[idx2])
            idx1 += 1
            idx2 += 1
        # Add gaps in seq1
        while idx1 < start1:
            seqA.append(s1[idx1])
            seqB.append("-")
            idx1 += 1
        # Add gaps in seq2
        while idx2 < start2:
            seqA.append("-")
            seqB.append(s2[idx2])
            idx2 += 1
        # Add the aligned block
        for i in range(end1 - start1):
            seqA.append(s1[idx1])
            seqB.append(s2[idx2])
            idx1 += 1
            idx2 += 1

    # Add remaining sequence
    while idx1 < len(s1) and idx2 < len(s2):
        seqA.append(s1[idx1])
        seqB.append(s2[idx2])
        idx1 += 1
        idx2 += 1
    while idx1 < len(s1):
        seqA.append(s1[idx1])
        seqB.append("-")
        idx1 += 1
    while idx2 < len(s2):
        seqA.append("-")
        seqB.append(s2[idx2])
        idx2 += 1

    return "".join(seqA), "".join(seqB)


def blosum_penalty(seq1, seq2, gap_open=-10, gap_extend=-1):
    # Load BLOSUM62 matrix
    blosum62 = substitution_matrices.load("BLOSUM62")

    seq1 = seq1.upper()
    seq2 = seq2.upper()

    total_score = 0
    in_gap1 = False
    in_gap2 = False

    i = 0  # Index for tracking positions
    print("blosum_penalty len seq1: ", len(seq1))
    while i < len(seq1):
        aa1 = seq1[i]
        aa2 = seq2[i]

        if aa1 == "-" and aa2 == "-":
            i += 1
            continue
        elif aa1 == "-":
            if in_gap1:
                total_score += gap_extend
            else:
                total_score += gap_open
                in_gap1 = True
            in_gap2 = False
        elif aa2 == "-":
            if in_gap2:
                total_score += gap_extend
            else:
                total_score += gap_open
                in_gap2 = True
            in_gap1 = False
        else:
            score = blosum62.get((aa1, aa2)) or blosum62.get((aa2, aa1))
            if score is None:
                raise ValueError(f"Invalid amino acid pair: \"{aa1}\", \"{aa2}\"")
            total_score += score
            in_gap1 = False
            in_gap2 = False
        i += 1

    return total_score


def process_mutation(record, apply_mut: bool = True):
    gene_name = record["gene"]
    transcript_id = record["transcript"].split(".")[0]  # Remove version number
    aa_mutation = record["aa_mutation"]

    result = {
        "gene": gene_name,
        "transcript": transcript_id,
        "aa_mutation": aa_mutation,
        "blosum_penalty": None,
        "error": None
    }

    try:
        # Fetch and translate the CDS sequence
        cds_seq = fetch_cds_sequence(transcript_id)
        if not cds_seq:
            error_msg = f"No CDS sequence found for transcript {transcript_id}"
            result["error"] = error_msg
            return result
        print("process_mutation translate start")
        protein_seq = translate_cds_to_protein(cds_seq)
        print("process_mutation translate end")

        # Apply the mutation
        print("process_mutation apply_mutation start")
        mutated_seq = apply_mutation(protein_seq, aa_mutation)
        # mutated_seq = protein_seq
        print("process_mutation apply_mutation end")

        # Align the sequences
        print("process_mutation align_sequences start")
        aligned_ref_seq, aligned_mut_seq, alignment_score = align_sequences(protein_seq, mutated_seq)
        real_aligned_ref_seq, real_aligned_mut_seq, real_alignment_score = align_sequences(protein_seq, protein_seq)
        print("process_mutation align_sequences end")

        # Calculate the penalty
        print("process_mutation blosum_penalty start")
        penalty = blosum_penalty(aligned_ref_seq, aligned_mut_seq)
        real_penalty = blosum_penalty(real_aligned_ref_seq, real_aligned_mut_seq)
        print("process_mutation blosum_penalty end")
        result["blosum_penalty"] = penalty
        result["real_blosum_penalty"] = real_penalty

    except Exception as e:
        result["error"] = str(e)

    return result


def cache_transcript_sequence(transcript_id):
    try:
        fetch_cds_sequence(transcript_id)
        logging.info(f"{transcript_id} done")
    except Exception as e:
        logging.error(f"{transcript_id} Error {e}")


from more_itertools import chunked

if __name__ == "__main__":
    # Read mutations from CSV
    csv_file = "data/unique_combs.csv"  # Replace with your CSV file path

    # old_output_file = "blosum_penalties3.csv"
    # output_file = "data/blosum_penalties_full.csv"

    # Read the mutations into a list of records
    # all_transcripts = set()
    # with open(csv_file, mode="r", newline="") as infile:
    #     reader = csv.DictReader(infile, fieldnames=["gene", "transcript", "aa_mutation"])
    #     for row in reader:
    #         all_transcripts.add(row["transcript"].split(".")[0])

    # true_transcripts = set()
    # try:
    #     with open(old_output_file, mode="r", newline="") as infile:
    #         reader = csv.DictReader(infile, fieldnames=["gene", "transcript", "aa_mutation", "blosum_penalty", "error"])
    #         for row in reader:

    #             if not row["error"] and row["blosum_penalty"]:
    #                 true_transcripts.add(row["transcript"])
    # except FileNotFoundError:
    #     pass
    # remaining_transcripts = all_transcripts.difference(true_transcripts)

    mutations = list()
    redis_conn = redis.Redis(host="localhost", port=6379, db=0, decode_responses=True)
    with open(csv_file, mode="r", newline="") as infile:
        reader = csv.DictReader(infile, fieldnames=["gene", "transcript", "aa_mutation"])
        for row in reader:
            transcript = row["transcript"].split(".")[0]
            if redis_conn.get(f"cds_sequence:{transcript}"):
                mutations.append(row)
    print(len(mutations))

    # with ProcessPoolExecutor(max_workers=16) as executor:
    #     future_to_record = {executor.submit(cache_transcript_sequence, transcript): transcript for transcript in all_transcripts}
    #     for future in as_completed(future_to_record):
    #         print(future.result())
    # executor.shutdown()

    # Open the output CSV file for writing
    chunk_size = 10_000
    for start in range(0, len(mutations), chunk_size):
        if start <= 1_630_000:
            continue
        end = start + chunk_size
        with open(f"data/blosum_penalties_full_{start}_{end}.csv", mode="w", newline="") as outfile:
            fieldnames = ["gene", "transcript", "aa_mutation", "blosum_penalty", "real_blosum_penalty", "error"]
            writer = csv.DictWriter(outfile, fieldnames=fieldnames)

            writer.writeheader()

            # Use ProcessPoolExecutor to process mutations in parallel
            with ProcessPoolExecutor(max_workers=5) as executor:
                # Submit all mutation processing tasks
                future_to_record = {executor.submit(process_mutation, record): record
                                    for record in mutations[start:end]}

                # Collect the results as they complete
                for future in as_completed(future_to_record):
                    result = future.result()
                    gene_name = result["gene"]
                    aa_mutation = result["aa_mutation"]
                    penalty = result.get("blosum_penalty")
                    real_penalty = result.get("real_blosum_penalty")
                    error = result.get("error")

                    if error:
                        print(f"Error processing mutation {aa_mutation} for gene {gene_name}: {error}")
                    else:
                        print(f"Gene: {gene_name}, Mutation: {aa_mutation}, BLOSUM Penalty: {penalty}")

                    # Write the result to the output CSV
                    # if result["transcript"] in remaining_transcripts:
                    writer.writerow(result)
