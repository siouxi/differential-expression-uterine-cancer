import pandas as pd
import os

files = [
    r"results/GSE285498/degs/CBD_vs_Eto_DEGs_FDR0.05_logFC1.tsv",
    r"results/GSE285498/degs/CBD_vs_Mock_DEGs_FDR0.05_logFC1.tsv",
    r"results/GSE285498/degs/Eto_vs_Mock_DEGs_FDR0.05_logFC1.tsv",
    r"results/GSE285498/degs/Merge_vs_Mock_DEGs_FDR0.05_logFC1.tsv"
]

base_path = r"c:\Users\nicol\Desktop\Proyectos\differential-expression-uterine-cancer"
output_base = os.path.join(base_path, "standard")
gse_id = "GSE285498"
output_dir = os.path.join(output_base, gse_id)
os.makedirs(output_dir, exist_ok=True)

for relative_path in files:
    full_path = os.path.join(base_path, relative_path)
    filename = os.path.basename(relative_path)
    filename_no_ext = os.path.splitext(filename)[0]
    
    output_xlsx = os.path.join(output_dir, filename_no_ext + ".xlsx")
    output_csv = os.path.join(output_dir, filename_no_ext + ".csv")
    
    print(f"Processing {filename}...")
    
    try:
        # Check if file is empty or too small
        if os.path.getsize(full_path) == 0:
             print(f"  Warning: File {filename} is empty. Skipping.")
             continue

        df = pd.read_csv(full_path, sep='\t')
        
        if df.empty:
            print(f"  Warning: DataFrame from {filename} is empty. Skipping.")
            continue

        # Rename ENSEMBL to gene_id
        if 'ENSEMBL' in df.columns:
            df.rename(columns={'ENSEMBL': 'gene_id'}, inplace=True)
        
        # Check columns
        cols_to_keep = ['gene_id', 'logFC']
        if all(col in df.columns for col in cols_to_keep):
            # Reorder columns explicitly: gene_id first, then logFC
            df_subset = df[cols_to_keep].copy()
            
            # Save as XLSX
            df_subset.to_excel(output_xlsx, index=False)
            print(f"  Saved to: {output_xlsx}")
            
            # Save as CSV
            df_subset.to_csv(output_csv, index=False)
            print(f"  Saved to: {output_csv}")
        else:
            print(f"  Error: Columns {cols_to_keep} not found (after renaming). Available: {list(df.columns)}")
            
    except Exception as e:
        print(f"  Error processing {filename}: {e}")
