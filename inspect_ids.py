import mygene
import json

# Initialize MyGene
mg = mygene.MyGeneInfo()

# IDs from the file
ids_to_check = ['3', '7', '15', '135', '145', '172', '192', '236']

results_entrez = mg.querymany(ids_to_check, scopes='entrezgene', fields='symbol,name,taxid', species='human')

with open('inspection_results.txt', 'w') as f:
    json.dump(results_entrez, f, indent=2)
