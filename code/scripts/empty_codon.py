import json

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

with open("/home/mohammad/Public/Drug-cell/DrugCell-public/codon_list_filtered.json", "r") as file:
    orig_data = json.load(file)
l1 = []
for g, v in data.items():
    for c, _ in v.items():
        l1.append((g, c))
l2 = []
for g, v in orig_data.items():
    for c, _ in v.items():
        l2.append((g, c))
print(set(l2) - set(l1))

# Function to find genes and codons with no data
def find_empty_genes_codons(data):
    empty_entries = []
    for gene, codons in data.items():
        for codon, details in codons.items():
            if not details:  # Check if details list is empty
                empty_entries.append((gene, codon))
    return empty_entries

# Call the function and print results
# empty_genes_codons = find_empty_genes_codons(data)
# if empty_genes_codons:
#     print("Genes and codons with no data:")
#     for gene, codon in empty_genes_codons:
#         print(f"Gene: {gene}, Codon: {codon}")
# else:
#     print("No empty data found.")

# print(data)
