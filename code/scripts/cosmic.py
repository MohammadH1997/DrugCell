import requests
cid = 905934
url = f"https://cancer.sanger.ac.uk/cell_lines/sample/mutations?id={cid}&export=json"

payload = {}
headers = {}

response = requests.request("GET", url, headers=headers, data=payload)
x = response.text
gene_mutations = list(set([(item.split("\t")[0].split("_")[0]) for item in x.split("\n") if item]))
gene_mutations.remove("Gene")
print(len(gene_mutations),gene_mutations)
with open (f"m{cid}_mutations.txt", "w") as fd:
    fd.write("\n".join(gene_mutations))
