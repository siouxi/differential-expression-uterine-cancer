import mygene
import pandas as pd
import os

# Paths
base_path = r"c:\Users\nicol\Desktop\Proyectos\differential-expression-uterine-cancer"
input_file = os.path.join(base_path, "standard", "GSE102787", "GSE102787_DEGS_Significativos.csv")
output_file = os.path.join(base_path, "standard", "GSE102787", "GSE102787_DEGS_Significativos_standard.csv")

print("Reading input file...")
df_input = pd.read_csv(input_file)
print(f"Input file has {len(df_input)} rows")
print(f"Columns: {list(df_input.columns)}")

# Get gene symbols
gene_symbols = df_input['Gene'].tolist()
print(f"\nQuerying MyGene.info for {len(gene_symbols)} gene symbols...")

# Initialize MyGene
mg = mygene.MyGeneInfo()

# Query MyGene.info
results = mg.querymany(
    gene_symbols,
    scopes='symbol',
    fields='ensembl.gene,entrezgene',
    species='human',
    returnall=True
)

# Process results
print(f"\nProcessing results...")
mapped_data = []

for gene_symbol in gene_symbols:
    # Find matching result
    match = None
    for r in results['out']:
        if r.get('query') == gene_symbol:
            match = r
            break
    
    if match and 'ensembl' in match:
        # Handle cases where ensembl is a list or dict
        ensembl_data = match['ensembl']
        if isinstance(ensembl_data, list):
            ensembl_id = ensembl_data[0].get('gene') if ensembl_data else None
        else:
            ensembl_id = ensembl_data.get('gene')
    else:
        ensembl_id = None
    
    entrez_id = match.get('entrezgene') if match else None
    
    mapped_data.append({
        'Ensembl_ID': ensembl_id,
        'Entrez_ID': int(entrez_id) if entrez_id and not pd.isna(entrez_id) else None,
        'Gene_Symbol': gene_symbol
    })

# Create output dataframe
df_mapped = pd.DataFrame(mapped_data)

# Merge with original data to get logFC
df_output = df_mapped.merge(df_input, left_on='Gene_Symbol', right_on='Gene', how='left')
df_output = df_output[['Ensembl_ID', 'Entrez_ID', 'Gene_Symbol', 'logFC']]

# Report mapping statistics
total_genes = len(df_output)
mapped_ensembl = df_output['Ensembl_ID'].notna().sum()
mapped_entrez = df_output['Entrez_ID'].notna().sum()

print(f"\nMapping Statistics:")
print(f"Total genes: {total_genes}")
print(f"Mapped to Ensembl ID: {mapped_ensembl} ({mapped_ensembl/total_genes*100:.1f}%)")
print(f"Mapped to Entrez ID: {mapped_entrez} ({mapped_entrez/total_genes*100:.1f}%)")
print(f"Unmapped genes (Ensembl): {total_genes - mapped_ensembl}")
print(f"Unmapped genes (Entrez): {total_genes - mapped_entrez}")

# Show unmapped genes
unmapped_ensembl = df_output[df_output['Ensembl_ID'].isna()]['Gene_Symbol'].tolist()
if unmapped_ensembl:
    print(f"\nUnmapped Gene Symbols (Ensembl) ({len(unmapped_ensembl)}):")
    print(", ".join(unmapped_ensembl[:15]))  # Show first 15
    if len(unmapped_ensembl) > 15:
        print(f"... and {len(unmapped_ensembl) - 15} more")

# Save output
df_output.to_csv(output_file, index=False)
print(f"\nSaved standardized file to: {output_file}")
