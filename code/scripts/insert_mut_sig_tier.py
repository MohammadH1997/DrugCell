import json
import os

import psycopg2

# Database connection details
DB_HOST = os.environ.get("DB_HOST", "localhost")
DB_PORT = os.environ.get("DB_PORT", "5432")
DB_NAME = os.environ.get("DB_NAME", "postgres")
DB_USER = os.environ.get("DB_USER", "postgres")
DB_PASSWORD = os.environ.get("DB_PASSWORD", "<PASSWORD>")

# Connect to PostgreSQL database
try:
    conn = psycopg2.connect(
        host=DB_HOST,
        port=DB_PORT,
        dbname=DB_NAME,
        user=DB_USER,
        password=DB_PASSWORD,
    )
    conn.autocommit = True
    cursor = conn.cursor()
    print("Database connection successful!")
except Exception as e:
    print(f"Error connecting to the database: {e}")
    exit()

# Load the JSON file
data = dict()
for i in range(5):
    with open(f"/home/mohammad/Public/Drug-cell/DrugCell-public/data/cosmic_tier_{i}.json", "r") as file:
        d = json.load(file)
        for k, v in d.items():
            if k not in data:
                data[k] = v
            else:
                data[k].update(v)

# Mutation type mapping
mapping = {
    "complex - deletion inframe": 1,
    "complex - frameshift": 2,
    "deletion - frameshift": 3,
    "deletion - in frame": 4,
    "insertion - frameshift": 5,
    "insertion - in frame": 6,
    "nonstop extension": 7,
    "substitution - coding silent": 8,
    "substitution - missense": 9,
    "substitution - nonsense": 10,
    "unknown": 11,
    "whole gene deletion": 12,
    "frameshift": 13,
    "complex - insertion inframe": 14,
}

# Insert data into the PostgreSQL table
for gene, d in data.items():
    for codon, boxes in d.items():
        for box in boxes:
            box_type = box["box_type"].lower()
            if not box_type:
                continue
            box_type_index = mapping.get(box_type, None)
            if box_type_index is None:
                print(f"Unknown box type: {box_type}")
                continue
            box_name = box["box_name"]
            tier = box["tier"]
            coding_sequence_mutation = box_name["1"]
            amino_acid_mutation = box_name["2"]

            try:
                # query = f"""
                #     INSERT INTO cosmic.cmc (gene, codon_number, mutation_type, coding_sequence_mutation, amino_acid_mutation, tier)
                #     VALUES ("{gene}", {codon}, {box_type_index}, "{coding_sequence_mutation}", "{amino_acid_mutation}", {tier})
                #     RETURNING id;
                #     """
                # if gene == "CTNNB1" and codon in {"13", "21", "22", "32", "33", "34", "37", "40", "41", "45"} and len(amino_acid_mutation) > 255:
                #     print(query)
                cursor.execute(
                    """
                    INSERT INTO cosmic.cmc (gene, codon_number, mutation_type, coding_sequence_mutation, amino_acid_mutation, tier)
                    VALUES (%s, %s, %s, %s, %s, %s)
                    """,
                    (gene, int(codon), box_type_index, coding_sequence_mutation, amino_acid_mutation, tier)
                )
            except Exception as e:
                print(f"Error inserting data for gene {gene}, codon {codon}: {e}")

# Close database connection
cursor.close()
conn.close()
print("Data insertion complete and database connection closed.")
