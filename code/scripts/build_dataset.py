import os

import psycopg2


def build_dataset():
    conn = psycopg2.connect(
        dbname=os.environ.get("POSTGRES_DB", "postgres"),
        user=os.environ.get("POSTGRES_USER", "postgres"),
        password=os.environ.get("POSTGRES_PASSWORD", "<PASSWORD>"),
        host=os.environ.get("POSTGRES_HOST", "localhost"),
        port=os.environ.get("POSTGRES_PORT", 5432),
    )
    cursor = conn.cursor()

    select_query = """
        SELECT base_gene_name as gene, transcript, amino_acid_mutation as aa_mutation 
        FROM cosmic.cell_line_mutations
        GROUP BY base_gene_name, transcript, amino_acid_mutation;
    """

    cursor.execute(select_query)
    rows = cursor.fetchall()
    print(len(rows))
    csv = ""
    for row in rows:
        csv += row[0] + "," + row[1] + "," + row[2] + "\n"
    with open("data/unique_combs.csv", "w") as fp:
        fp.write(csv)
    conn.commit()

    cursor.close()
    conn.close()


build_dataset()
