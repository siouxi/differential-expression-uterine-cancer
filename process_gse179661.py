import pandas as pd
import os

files_info = [
    {
        "path": r"results/GSE179661/HepG2/deseq2/DESeq2_DEG_FDR0.05_log2FC2.tsv",
        "cell_line": "HepG2"
    },
    {
        "path": r"results/GSE179661/MHCC97H/deseq2/DESeq2_DEG_FDR0.05_log2FC2.tsv",
        "cell_line": "MHCC97H"
    }
]

base_path = r"c:\Users\nicol\Desktop\Proyectos\differential-expression-uterine-cancer"
output_base = os.path.join(base_path, "standard")
gse_id = "GSE179661"
output_dir = os.path.join(output_base, gse_id)
os.makedirs(output_dir, exist_ok=True)

for info in files_info:
    relative_path = info["path"]
    cell_line = info["cell_line"]
    full_path = os.path.join(base_path, relative_path)
    
    filename = os.path.basename(relative_path)
    filename_no_ext = os.path.splitext(filename)[0]
    
    # Append cell line to filename
    new_filename = f"{filename_no_ext}_{cell_line}"
    
    output_xlsx = os.path.join(output_dir, new_filename + ".xlsx")
    output_csv = os.path.join(output_dir, new_filename + ".csv")
    
    print(f"Processing {relative_path} ({cell_line})...")
    
    try:
        df = pd.read_csv(full_path, sep='\t')
        
        cols_to_keep = ['gene_id', 'log2FoldChange']
        if all(col in df.columns for col in cols_to_keep):
            df_subset = df[cols_to_keep].copy()
            
            # Clean gene_id (remove .0 if present)
            df_subset['gene_id'] = pd.to_numeric(df_subset['gene_id'], errors='coerce')
            df_subset = df_subset.dropna(subset=['gene_id'])
            df_subset['gene_id'] = df_subset['gene_id'].astype(int).astype(str)
            
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
