import csv
import os
from concurrent.futures import as_completed

import psycopg2
from redis import Redis


def read_data_from_file(file_path):
    data = []

    with open(file_path, mode="r", newline="") as infile:
        reader = csv.reader(infile)
        for row in reader:
            if row[5]:
                continue
            else:
                # Assuming the file has columns: gene_name, aa_mutation, transcript, blosum_penalty
                gene_name = row[0]
                transcript = row[1]
                aa_mutation = row[2]
                blosum_penalty = int(float(row[3]))  # Assuming the penalty is in the fourth column
                real_blosum_penalty = int(float(row[4]))  # Assuming the penalty is in the fourth column
                data.append((blosum_penalty, real_blosum_penalty, gene_name, transcript))
    print(" data add ")
    return data


def update_blosum_penalty(file_path):
    from psycopg2.pool import ThreadedConnectionPool
    from psycopg2.extras import RealDictCursor
    from concurrent.futures import ThreadPoolExecutor

    pool = ThreadPoolExecutor(max_workers=10)
    conn_pool = ThreadedConnectionPool(minconn=1,
                                       maxconn=10,
                                       dbname=os.environ["DB_NAME"],
                                       user=os.environ["DB_USER"],
                                       password=os.environ["DB_PASSWORD"],
                                       host=os.environ["DB_HOST"],
                                       port=os.environ["DB_PORT"],
                                       )

    def _update(data, i=0):
        update_query = """
            UPDATE cosmic.cell_line_mutations
            SET aa_mutation_blosum62_penalty = %s, real_penalty = %s
            WHERE base_gene_name = %s AND base_transcript = %s 
            RETURNING id;
        """
        c = conn_pool.getconn()
        cursor = c.cursor()
        for row in data:
            cursor.execute(update_query, row)
        c.commit()
        print(f"task {i} done:", cursor.fetchall())
        conn_pool.putconn(c)
        # print(f"task {i} done")

    # Read data from CSV file
    data_to_update = read_data_from_file(file_path)
    print(len(data_to_update))

    # for count, row in enumerate(data_to_update):
    #     if count < 104728:
    #         continue
    #     else :
    #         cursor.execute(update_query, row)
    #         if count % 100 == 0:
    #             conn.commit()
    #             print("commit:" + str(count))
    # conn.commit()
    chunks = [data_to_update[i:i + 100] for i in range(0, len(data_to_update), 100)]
    tasks = [pool.submit(_update, data=chunk, i=i) for i, chunk in enumerate(chunks)]
    for task in as_completed(tasks):
        task.result()

    # cursor.close()
    # conn.close()


import os
import glob

files = glob.glob("data/blosum_penalties_full_*.csv")

# Call the function with the path to your file
for file in files:
    update_blosum_penalty(file)
# update_transcript_seq_length()
