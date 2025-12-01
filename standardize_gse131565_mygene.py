import mygene
import pandas as pd
import os

# Paths
base_path = r"c:\Users\nicol\Desktop\Proyectos\differential-expression-uterine-cancer"
input_file = os.path.join(base_path, "standard", "GSE131565", "DESeq2_CBD_vs_Control.csv")
output_file = os.path.join(base_path, "standard", "GSE131565", "DESeq2_CBD_vs_Control_standard.csv")

print("Reading input file...")
df_input = pd.read_csv(input_file)
print(f"Input file has {len(df_input)} rows")
print(f"Columns: {list(df_input.columns)}")

# Get entrez IDs (they're already integers)
entrez_ids = df_input['gene_id'].tolist()
print(f"\nQuerying MyGene.info for {len(entrez_ids)} Entrez IDs...")

# Initialize MyGene
mg = mygene.MyGeneInfo()

# Query MyGene.info - use entrezgene as the scope
results = mg.querymany(
    entrez_ids,
    scopes='entrezgene',
    fields='ensembl.gene,symbol',
    species='human',
    returnall=True
)

# Process results
print(f"\nProcessing results...")
mapped_data = []

for entrez_id in entrez_ids:
    # Find matching result
    match = None
    for r in results['out']:
        if r.get('query') == str(entrez_id):
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
    
    gene_symbol = match.get('symbol') if match else None
    
    mapped_data.append({
        'Ensembl_ID': ensembl_id,
        'Entrez_ID': entrez_id,
        'Gene_Symbol': gene_symbol
    })

# Create output dataframe
df_mapped = pd.DataFrame(mapped_data)

# Merge with original data to get log2FoldChange
df_output = df_mapped.merge(df_input, left_on='Entrez_ID', right_on='gene_id', how='left')
df_output = df_output[['Ensembl_ID', 'Entrez_ID', 'Gene_Symbol', 'log2FoldChange']]

# Report mapping statistics
total_genes = len(df_output)
mapped_ensembl = df_output['Ensembl_ID'].notna().sum()
mapped_symbol = df_output['Gene_Symbol'].notna().sum()

print(f"\nMapping Statistics:")
print(f"Total genes: {total_genes}")
print(f"Mapped to Ensembl ID: {mapped_ensembl} ({mapped_ensembl/total_genes*100:.1f}%)")
print(f"Mapped to Gene Symbol: {mapped_symbol} ({mapped_symbol/total_genes*100:.1f}%)")
print(f"Unmapped genes (Ensembl): {total_genes - mapped_ensembl}")
print(f"Unmapped genes (Symbol): {total_genes - mapped_symbol}")

# Show unmapped genes
unmapped_ensembl = df_output[df_output['Ensembl_ID'].isna()]['Entrez_ID'].tolist()
if unmapped_ensembl:
    print(f"\nUnmapped Entrez IDs (Ensembl) ({len(unmapped_ensembl)}):")
    print(", ".join(map(str, unmapped_ensembl[:15])))  # Show first 15
    if len(unmapped_ensembl) > 15:
        print(f"... and {len(unmapped_ensembl) - 15} more")

# Save output
df_output.to_csv(output_file, index=False)
print(f"\nSaved standardized file to: {output_file}")
