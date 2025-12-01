import pandas as pd
import os

files = [
    r"results/GSE223827/deg/DESeq2_CBD_vs_DMSO.tsv",
    r"results/GSE223827/deg/DESeq2_CBD-THC_vs_DMSO.tsv",
    r"results/GSE223827/deg/DESeq2_DMSO-EtOH_vs_DMSO.tsv",
    r"results/GSE223827/deg/DESeq2_EtOH_vs_DMSO.tsv",
    r"results/GSE223827/deg/DESeq2_THC_vs_DMSO.tsv"
]

base_path = r"c:\Users\nicol\Desktop\Proyectos\differential-expression-uterine-cancer"
output_base = os.path.join(base_path, "standard")
gse_id = "GSE223827"
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
        df = pd.read_csv(full_path, sep='\t')
        
        cols_to_keep = ['gene_id', 'log2FoldChange']
        if all(col in df.columns for col in cols_to_keep):
            df_subset = df[cols_to_keep].copy()
            
            # Save as XLSX
            df_subset.to_excel(output_xlsx, index=False)
            print(f"  Saved to: {output_xlsx}")
            
            # Save as CSV
            df_subset.to_csv(output_csv, index=False)
            print(f"  Saved to: {output_csv}")
        else:
            print(f"  Error: Columns {cols_to_keep} not found in {filename}")
            
    except Exception as e:
        print(f"  Error processing {filename}: {e}")
