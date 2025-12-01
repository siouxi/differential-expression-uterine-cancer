import mygene
import pandas as pd
import json

output = []

# 1. Check local file
output.append("--- Checking local file ---")
try:
    df = pd.read_csv('results/GSE285498/gene_annotation_ensembl_human.tsv', sep='\t')
    df['ENTREZID'] = df['ENTREZID'].astype(str)
    
    found_3 = df[df['ENTREZID'] == '3']
    found_7 = df[df['ENTREZID'] == '7']
    
    output.append(f"Found ID 3 in local file: {not found_3.empty}")
    if not found_3.empty:
        output.append(str(found_3.to_dict(orient='records')))
        
    output.append(f"Found ID 7 in local file: {not found_7.empty}")
    if not found_7.empty:
        output.append(str(found_7.to_dict(orient='records')))
        
except Exception as e:
    output.append(f"Error reading local file: {e}")

# 2. Query MyGene for symbols
output.append("\n--- Querying MyGene for symbols A2M, AATK ---")
mg = mygene.MyGeneInfo()
symbols = ['A2M', 'AATK']
results = mg.querymany(symbols, scopes='symbol', fields='entrezgene,name,taxid', species='human')
output.append(json.dumps(results, indent=2))

with open('check_results.txt', 'w') as f:
    f.write('\n'.join(output))
