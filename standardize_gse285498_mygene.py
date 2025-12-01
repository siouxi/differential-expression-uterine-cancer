import mygene
import pandas as pd
import os

# Files to process
files = [
    "CBD_vs_Eto_DEGs_FDR0.05_logFC1.csv",
    "Eto_vs_Mock_DEGs_FDR0.05_logFC1.csv",
    "Merge_vs_Mock_DEGs_FDR0.05_logFC1.csv"
]

base_path = r"c:\Users\nicol\Desktop\Proyectos\differential-expression-uterine-cancer"
input_dir = os.path.join(base_path, "standard", "GSE285498")
annotation_file = os.path.join(base_path, "results", "GSE285498", "gene_annotation_ensembl_human.tsv")

# Initialize MyGene
mg = mygene.MyGeneInfo()

# Load local annotation file
print(f"Loading local annotation file: {annotation_file}")
local_map = {}
try:
    df_local = pd.read_csv(annotation_file, sep='\t')
    # Create mapping: EnsemblID -> {EntrezID, Symbol}
    for _, row in df_local.iterrows():
        ensembl = row['ENSEMBL_ID']
        if pd.notna(ensembl):
            local_map[ensembl] = {
                'entrez': str(int(row['ENTREZID'])) if pd.notna(row['ENTREZID']) else None,
                'symbol': row['SYMBOL']
            }
    print(f"Loaded {len(local_map)} entries from local file.\n")
except Exception as e:
    print(f"Error loading local annotation file: {e}")
    print("Proceeding with MyGene only.\n")

for filename in files:
    input_file = os.path.join(input_dir, filename)
    output_file = os.path.join(input_dir, filename.replace(".csv", "_standard.csv"))
    
    print(f"{'='*60}")
    print(f"Processing {filename}...")
    print(f"{'='*60}")
    
    # Read input file
    df_input = pd.read_csv(input_file)
    print(f"Input file has {len(df_input)} rows")
    
    # Get Ensembl IDs from gene_id column
    ensembl_ids = df_input['gene_id'].tolist()
    
    print(f"Querying MyGene.info for {len(ensembl_ids)} Ensembl IDs...")
    
    # Query MyGene.info
    try:
        results = mg.querymany(
            ensembl_ids,
            scopes='ensembl.gene',
            fields='entrezgene,symbol',
            species='human',
            returnall=True
        )
        mygene_results = results['out']
    except Exception as e:
        print(f"MyGene query failed: {e}")
        mygene_results = []
    
    # Process results
    print(f"Processing results...")
    mapped_data = []
    
    stats = {
        'total': len(ensembl_ids),
        'mapped_mygene': 0,
        'mapped_local': 0,
        'unmapped': 0
    }
    
    # Create a lookup for MyGene results
    mg_lookup = {}
    for r in mygene_results:
        q = r.get('query')
        if q:
            mg_lookup[q] = r

    for ensembl_id in ensembl_ids:
        entrez_id = None
        gene_symbol = None
        source = None
        
        # 1. Try MyGene
        match = mg_lookup.get(ensembl_id)
        if match and 'notfound' not in match:
            # Check for Entrez ID
            if 'entrezgene' in match:
                entrez_id = str(int(match['entrezgene']))
            
            # Check for Symbol
            if 'symbol' in match:
                gene_symbol = match['symbol']
                
            if entrez_id or gene_symbol:
                source = 'mygene'
        
        # 2. Try Local File if missing data
        if not entrez_id or not gene_symbol:
            local_entry = local_map.get(ensembl_id)
            if local_entry:
                if not entrez_id:
                    entrez_id = local_entry['entrez']
                if not gene_symbol:
                    gene_symbol = local_entry['symbol']
                if entrez_id or gene_symbol:
                    source = 'local' if source != 'mygene' else 'hybrid'
        
        # Update stats
        if source == 'mygene':
            stats['mapped_mygene'] += 1
        elif source in ['local', 'hybrid']:
            stats['mapped_local'] += 1
        else:
            stats['unmapped'] += 1

        mapped_data.append({
            'Ensembl_ID': ensembl_id,
            'Entrez_ID': entrez_id,
            'Gene_Symbol': gene_symbol
        })
    
    # Create output dataframe
    df_mapped = pd.DataFrame(mapped_data)
    
    # Merge with original data to get logFC
    df_output = df_mapped.copy()
    df_output['logFC'] = df_input['logFC'].values
    
    # Reorder columns
    df_output = df_output[['Ensembl_ID', 'Entrez_ID', 'Gene_Symbol', 'logFC']]
    
    # Report mapping statistics
    print(f"\nMapping Statistics:")
    print(f"Total genes: {stats['total']}")
    print(f"Mapped by MyGene: {stats['mapped_mygene']} ({stats['mapped_mygene']/stats['total']*100:.1f}%)")
    print(f"Mapped by Local File: {stats['mapped_local']} ({stats['mapped_local']/stats['total']*100:.1f}%)")
    total_mapped = stats['mapped_mygene'] + stats['mapped_local']
    print(f"Total Mapped: {total_mapped} ({total_mapped/stats['total']*100:.1f}%)")
    print(f"Unmapped: {stats['unmapped']} ({stats['unmapped']/stats['total']*100:.1f}%)")
    
    # Count how many have both IDs
    both_ids = df_output[(df_output['Ensembl_ID'].notna()) & (df_output['Entrez_ID'].notna())]
    print(f"Genes with both Ensembl and Entrez: {len(both_ids)} ({len(both_ids)/stats['total']*100:.1f}%)")
    
    # Save output
    df_output.to_csv(output_file, index=False)
    print(f"Saved to: {output_file}\n")

print(f"{'='*60}")
print("All files processed successfully!")
print(f"{'='*60}")
