import pandas as pd
import os

files = [
    r"results/GSE173201/deseq2/DEGs_A2780_Cisplatino_vs_Control_FDR0.05_log2FC1.csv",
    r"results/GSE173201/deseq2/DEGs_A2780cis_Cisplatino_vs_A2780cis_Control_FDR0.05_log2FC1.csv",
    r"results/GSE173201/deseq2/DEGs_A2780cis_Control_vs_A2780_Control_FDR0.05_log2FC1.csv"
]

base_path = r"c:\Users\nicol\Desktop\Proyectos\differential-expression-uterine-cancer"
output_base = os.path.join(base_path, "standard")
gse_id = "GSE173201"
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
        df = pd.read_csv(full_path)
        
        # Rename Unnamed: 0 to gene_id
        if 'Unnamed: 0' in df.columns:
            df.rename(columns={'Unnamed: 0': 'gene_id'}, inplace=True)
            
        # Check columns
        cols_to_keep = ['gene_id', 'log2FoldChange']
        if all(col in df.columns for col in cols_to_keep):
            df_subset = df[cols_to_keep].copy()
            
            # Clean gene_id (remove .0 if present)
            # Convert to numeric first to handle potential mixed types, then to string/int
            df_subset['gene_id'] = pd.to_numeric(df_subset['gene_id'], errors='coerce')
            df_subset = df_subset.dropna(subset=['gene_id']) # Drop if gene_id became NaN
            df_subset['gene_id'] = df_subset['gene_id'].astype(int).astype(str)
            
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
