import json
import time

import pandas as pd
import requests

urls = []

with open("data/urls.txt", "r") as file:
    for line in file:
        urls.append(line.strip())

gene_list = []
with open("data/gene_API.txt", "r") as file:
    for line in file:
        gene_list.append(line.strip())

for i in range(len(urls)):
    if i <= 2725:
        print("skip " + str(i))
        continue
    url = urls[i]

    payload = {}
    headers = {
        "accept": "application/json, text/javascript, */*; q=0.01",
        "user-agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/123.0.0.0 Safari/537.36"
    }
    time.sleep(4)
    response = requests.request(method="GET",
                                url=url,
                                headers=headers,
                                data=payload,
                                cookies={
                                    "ncbi_sid": "",
                                    "_ga": "GA1.1.1029951559.1729864783",
                                    "WebEnv": "1dAvCxsyeuPTqv3Kt-4kIT2Y4M_hR4eo2y5cjlfif9B5s%40CE8B07AC71BA4491_0359SID",
                                    "ncbi_pinger": "",
                                })

    print(response.text)

    # Parse JSON-like structure from string content
    try:
        # Assuming the data starts with a large JSON object
        json_data = json.loads(response.text)

        # Extract data based on observed structure
        extracted_data = []
        for item in json_data.get("vars", []):
            # Add extracted fields from each item
            extracted_data.append({
                "id": item.get("id", None),
                "gene": (gene_list[i], None),
                "ci": item.get("ci", None),
                "revstat": item.get("revstat", None),
                "desc": item.get("desc", None),
                "loc_from": item["locs"][0].get("from", None) if item.get("locs") else None,
                "loc_to": item["locs"][0].get("to", None) if item.get("locs") else None,
                "loc_assm": item["locs"][0].get("assm", None) if item.get("locs") else None,
                "loc_chr": item["locs"][0].get("chr", None) if item.get("locs") else None,
            })

        # Convert to DataFrame and save as CSV
        df_extracted = pd.DataFrame(extracted_data)

        out = df_extracted.to_csv(path_or_buf="clinvar_data.csv",
                                  mode="a",
                                  header=not pd.io.common.file_exists("clinvar_data.csv"),
                                  index=False)
        print(i)


    except json.JSONDecodeError as e:
        output_path = f"Error parsing JSON: {e}"

    print(out)
