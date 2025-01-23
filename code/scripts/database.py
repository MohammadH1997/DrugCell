import os

import psycopg2
import json

# Connect to PostgreSQL
connection = psycopg2.connect(
    dbname=os.environ.get("DB_NAME"),
    user=os.environ.get("DB_USER"),
    password=os.environ.get("DB_PASSWORD"),
    host=os.environ.get("DB_HOST"),
    port=os.environ.get("DB_PORT"),
)
cursor = connection.cursor()

# Sample data (replace this with your actual JSON)
data = {
    "ABL1": {
        "255": [
            {"box_type": "Substitution - Missense", "box_name": {"0": "COSV107386955", "1": "c.763G>A", "2": "p.E255K"}, "tier": 3},
            {"box_type": "Substitution - coding silent", "box_name": {"0": "COSV100583676", "1": "c.765G>A", "2": "p.E255="}, "tier": 4},
            {"box_type": "Deletion - In frame", "box_name": {"0": "COSV59334390", "1": "c.607_879del", "2": "p.L203_K293del"}, "tier": 4}
        ],
        "263": [
            {"box_type": "Substitution - Missense", "box_name": {"0": "COSV59334206", "1": "c.787A>G", "2": "p.M263V"}, "tier": 3},
            {"box_type": "Deletion - In frame", "box_name": {"0": "COSV59334390", "1": "c.607_879del", "2": "p.L203_K293del"}, "tier": 4}
        ]
    }
}

# Insert data into the database
for gene, positions in data.items():
    for position, entries in positions.items():
        for entry in entries:
            box_type = entry["box_type"]
            box_name = json.dumps(entry["box_name"])  # Convert dict to JSON string
            tier = entry["tier"]
            query = """
            INSERT INTO mutations (gene, position, box_type, box_name, tier)
            VALUES (%s, %s, %s, %s, %s);
            """
            cursor.execute(query, (gene, int(position), box_type, box_name, tier))

# Commit and close
connection.commit()
cursor.close()
connection.close()
