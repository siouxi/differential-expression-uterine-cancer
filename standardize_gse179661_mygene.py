import mygene
import pandas as pd
import os
import sys

# Files to process
files = [
    "DESeq2_DEG_FDR0.05_log2FC2_HepG2.csv",
    "DESeq2_DEG_FDR0.05_log2FC2_MHCC97H.csv"
]

base_path = r"c:\Users\nicol\Desktop\Proyectos\differential-expression-uterine-cancer"
input_dir = os.path.join(base_path, "standard", "GSE179661")
annotation_file = os.path.join(base_path, "results", "GSE285498", "gene_annotation_ensembl_human.tsv")

# Initialize MyGene
mg = mygene.MyGeneInfo()

# Load local annotation file
print(f"Loading local annotation file: {annotation_file}")
local_map = {}
try:
    df_local = pd.read_csv(annotation_file, sep='\t')
    # Create mapping: EntrezID -> {EnsemblID, Symbol}
    # Ensure EntrezID is string for matching
    df_local['ENTREZID'] = df_local['ENTREZID'].astype(str)
    
    for _, row in df_local.iterrows():
        entrez = row['ENTREZID']
        if pd.notna(entrez) and entrez != 'nan':
            local_map[entrez] = {
                'ensembl': row['ENSEMBL_ID'],
                'symbol': row['SYMBOL']
            }
    print(f"Loaded {len(local_map)} entries from local file.")
except Exception as e:
    print(f"Error loading local annotation file: {e}")
    print("Proceeding with MyGene only.")

for filename in files:
    input_file = os.path.join(input_dir, filename)
    output_file = os.path.join(input_dir, filename.replace(".csv", "_standard.csv"))
    
    print(f"\n{'='*60}")
    print(f"Processing {filename}...")
    print(f"{'='*60}")
    
    if not os.path.exists(input_file):
        print(f"File not found: {input_file}")
        continue

    # Read input file
    df_input = pd.read_csv(input_file)
    print(f"Input file has {len(df_input)} rows")
    
    # Get entrez IDs
    # Ensure they are strings
    df_input['gene_id'] = df_input['gene_id'].astype(str)
    entrez_ids = df_input['gene_id'].tolist()
    
    print(f"Querying MyGene.info for {len(entrez_ids)} Entrez IDs...")
    
    # Query MyGene.info
    try:
        results = mg.querymany(
            entrez_ids,
            scopes='entrezgene',
            fields='ensembl.gene,symbol',
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
        'total': len(entrez_ids),
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

    for entrez_id in entrez_ids:
        ensembl_id = None
        gene_symbol = None
        source = None
        
        # 1. Try MyGene
        match = mg_lookup.get(entrez_id)
        if match:
            # Check for Ensembl ID
            if 'ensembl' in match:
                ensembl_data = match['ensembl']
                if isinstance(ensembl_data, list):
                    ensembl_id = ensembl_data[0].get('gene') if ensembl_data else None
                else:
                    ensembl_id = ensembl_data.get('gene')
            
            # Check for Symbol
            if 'symbol' in match:
                gene_symbol = match['symbol']
                
            if ensembl_id:
                source = 'mygene'
        
        # 2. Try Local File if Ensembl ID missing
        if not ensembl_id:
            local_entry = local_map.get(entrez_id)
            if local_entry:
                ensembl_id = local_entry['ensembl']
                # Only use local symbol if we don't have one yet
                if not gene_symbol:
                    gene_symbol = local_entry['symbol']
                source = 'local'
        
        # Update stats
        if source == 'mygene':
            stats['mapped_mygene'] += 1
        elif source == 'local':
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
    
    # Merge with original data to get log2FoldChange
    # Ensure types match for merge
    df_mapped['Entrez_ID'] = df_mapped['Entrez_ID'].astype(str)
    df_input['gene_id'] = df_input['gene_id'].astype(str)
    
    df_output = df_mapped.merge(df_input, left_on='Entrez_ID', right_on='gene_id', how='left')
    df_output = df_output[['Ensembl_ID', 'Entrez_ID', 'Gene_Symbol', 'log2FoldChange']]
    
    # Report mapping statistics
    print(f"\nMapping Statistics:")
    print(f"Total genes: {stats['total']}")
    print(f"Mapped by MyGene: {stats['mapped_mygene']} ({stats['mapped_mygene']/stats['total']*100:.1f}%)")
    print(f"Mapped by Local File: {stats['mapped_local']} ({stats['mapped_local']/stats['total']*100:.1f}%)")
    print(f"Total Mapped (Ensembl): {stats['mapped_mygene'] + stats['mapped_local']} ({(stats['mapped_mygene'] + stats['mapped_local'])/stats['total']*100:.1f}%)")
    print(f"Unmapped: {stats['unmapped']} ({stats['unmapped']/stats['total']*100:.1f}%)")
    
    # Save output
    df_output.to_csv(output_file, index=False)
    print(f"Saved to: {output_file}")

print(f"\n{'='*60}")
print("All files processed successfully!")
print(f"{'='*60}")
