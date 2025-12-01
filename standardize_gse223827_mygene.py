import mygene
import pandas as pd
import os

# Files to process
files = [
    "DESeq2_CBD_vs_DMSO.csv",
    "DESeq2_CBD-THC_vs_DMSO.csv",
    "DESeq2_DMSO-EtOH_vs_DMSO.csv",
    "DESeq2_EtOH_vs_DMSO.csv",
    "DESeq2_THC_vs_DMSO.csv"
]

base_path = r"c:\Users\nicol\Desktop\Proyectos\differential-expression-uterine-cancer"
input_dir = os.path.join(base_path, "standard", "GSE223827")
annotation_file = os.path.join(base_path, "results", "GSE285498", "gene_annotation_ensembl_human.tsv")

# Initialize MyGene
mg = mygene.MyGeneInfo()

# Load local annotation file
print(f"Loading local annotation file: {annotation_file}")
local_map = {}
try:
    df_local = pd.read_csv(annotation_file, sep='\t')
    # Create mapping: Symbol -> {EnsemblID, EntrezID}
    for _, row in df_local.iterrows():
        symbol = row['SYMBOL']
        if pd.notna(symbol):
            local_map[symbol] = {
                'ensembl': row['ENSEMBL_ID'],
                'entrez': str(int(row['ENTREZID'])) if pd.notna(row['ENTREZID']) else None
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
    
    # Get gene symbols from gene_id column
    gene_symbols = df_input['gene_id'].tolist()
    
    print(f"Querying MyGene.info for {len(gene_symbols)} Gene Symbols...")
    
    # Query MyGene.info
    try:
        results = mg.querymany(
            gene_symbols,
            scopes='symbol',
            fields='ensembl.gene,entrezgene,symbol',
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
        'total': len(gene_symbols),
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

    for gene_symbol in gene_symbols:
        ensembl_id = None
        entrez_id = None
        final_symbol = gene_symbol  # Keep original symbol
        source = None
        
        # 1. Try MyGene
        match = mg_lookup.get(gene_symbol)
        if match and 'notfound' not in match:
            # Check for Ensembl ID
            if 'ensembl' in match:
                ensembl_data = match['ensembl']
                if isinstance(ensembl_data, list):
                    ensembl_id = ensembl_data[0].get('gene') if ensembl_data else None
                else:
                    ensembl_id = ensembl_data.get('gene')
            
            # Check for Entrez ID
            if 'entrezgene' in match:
                entrez_id = str(int(match['entrezgene']))
            
            # Update symbol if provided
            if 'symbol' in match:
                final_symbol = match['symbol']
                
            if ensembl_id or entrez_id:
                source = 'mygene'
        
        # 2. Try Local File if missing data
        if not ensembl_id or not entrez_id:
            local_entry = local_map.get(gene_symbol)
            if local_entry:
                if not ensembl_id:
                    ensembl_id = local_entry['ensembl']
                if not entrez_id:
                    entrez_id = local_entry['entrez']
                if ensembl_id or entrez_id:
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
            'Gene_Symbol': final_symbol
        })
    
    # Create output dataframe
    df_mapped = pd.DataFrame(mapped_data)
    
    # Merge with original data to get log2FoldChange
    df_output = df_mapped.copy()
    df_output['log2FoldChange'] = df_input['log2FoldChange'].values
    
    # Reorder columns
    df_output = df_output[['Ensembl_ID', 'Entrez_ID', 'Gene_Symbol', 'log2FoldChange']]
    
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
