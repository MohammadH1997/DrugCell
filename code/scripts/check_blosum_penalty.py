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
            if row[4]:
                continue
            else:
                # Assuming the file has columns: gene_name, aa_mutation, transcript, blosum_penalty
                gene_name = row[0]
                aa_mutation = row[2]
                transcript = row[1]
                blosum_penalty = int(float(row[3]))  # Assuming the penalty is in the fourth column
                data.append((blosum_penalty, gene_name, transcript))
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

    # conn = psycopg2.connect(
    #     dbname="postgres",
    #     user="root",
    #     password="Muy1sTSNpnkjXalbQwVUDnBE",
    #     host="etna.liara.cloud",
    #     port="34965"
    # )
    # cursor = conn.cursor()

    def _update(data, i=0):
        update_query = """
            UPDATE cosmic.cell_line_mutations
            SET real_penalty = %s
            WHERE base_gene_name = %s AND base_transcript = %s ;
        """
        c = conn_pool.getconn()
        cursor = c.cursor()
        for row in data:
            cursor.execute(update_query, row)
        c.commit()
        conn_pool.putconn(c)
        print(f"task {i} done")

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


def update_transcript_seq_length():
    redis_conn = Redis(host="localhost", port=6379, db=0)

    # s = set(
    #     """""".split(",")
    # )

    # Define the prefix to remove
    prefix = "cds_sequence:"
    prefix_length = len(prefix)

    # Iterate over the matched keys and extract IDs
    var_list = [
        (redis_conn.strlen(x), x.decode("utf-8")[prefix_length:])
        for x in redis_conn.scan_iter(match="cds_sequence:*")
        if x.decode("utf-8").startswith(prefix)
    ]

    print(len(var_list), var_list[:10])

    conn = psycopg2.connect(
        dbname=os.environ.get("POSTGRES_DB", "postgres"),
        user=os.environ.get("POSTGRES_USER", "postgres"),
        password=os.environ.get("POSTGRES_PASSWORD", "<PASSWORD>"),
        host=os.environ.get("POSTGRES_HOST", "localhost"),
        port=os.environ.get("POSTGRES_PORT", 5432),
    )
    cursor = conn.cursor()

    update_query = """
        UPDATE cosmic.cell_line_mutations
        SET transcript_sequence_length = %s
        WHERE base_transcript = %s ;
    """

    # cursor.executemany(update_query, var_list)
    # conn.commit()
    # cursor.close()

    for count, row in enumerate(var_list):
        cursor.execute(update_query, row)
        if count % 100 == 0:
            conn.commit()
            print("commit:" + str(count))
    conn.commit()

    cursor.close()
    conn.close()


# Call the function with the path to your file
file_path = "/home/mohammad/Public/Drug-cell/DrugCell-public/blosum_real_penalties2.csv"
update_blosum_penalty(file_path)
# update_transcript_seq_length()
