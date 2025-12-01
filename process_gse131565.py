import pandas as pd
import os

# Input file
input_file = r"results/GSE131565/deg/DESeq2_CBD_vs_Control.tsv"
base_path = r"c:\Users\nicol\Desktop\Proyectos\differential-expression-uterine-cancer"
full_input_path = os.path.join(base_path, input_file)

# Output setup
output_base = os.path.join(base_path, "standard")
gse_id = "GSE131565"
output_dir = os.path.join(output_base, gse_id)
os.makedirs(output_dir, exist_ok=True)

filename = os.path.basename(input_file)
filename_no_ext = os.path.splitext(filename)[0]
output_xlsx = os.path.join(output_dir, filename_no_ext + ".xlsx")
output_csv = os.path.join(output_dir, filename_no_ext + ".csv")

print(f"Processing {filename}...")

try:
    # Read TSV
    df = pd.read_csv(full_input_path, sep='\t')
    
    # Check columns
    cols_to_keep = ['gene_id', 'log2FoldChange']
    if all(col in df.columns for col in cols_to_keep):
        df_subset = df[cols_to_keep]
        
        # Save as XLSX
        df_subset.to_excel(output_xlsx, index=False)
        print(f"  Saved to: {output_xlsx}")
        
        # Save as CSV
        df_subset.to_csv(output_csv, index=False)
        print(f"  Saved to: {output_csv}")
    else:
        print(f"  Error: Columns {cols_to_keep} not found. Available columns: {list(df.columns)}")

except Exception as e:
    print(f"  Error processing {filename}: {e}")
