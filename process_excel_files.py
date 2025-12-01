import pandas as pd
import os
import re

# Input files
files = [
    r"results/DEGS_ANGELICA/GSE21656_DEGS_significativos.xlsx",
    r"results/DEGS_ANGELICA/GSE102787_DEGS_Significativos.xlsx",
    r"results/DEGS_ANGELICA/GSE197561_DEGS_significativos.xlsx"
]

base_path = r"c:\Users\nicol\Desktop\Proyectos\differential-expression-uterine-cancer"
output_base = os.path.join(base_path, "standard")

# Ensure output base exists
os.makedirs(output_base, exist_ok=True)

for relative_path in files:
    full_path = os.path.join(base_path, relative_path)
    filename = os.path.basename(relative_path)
    
    # Extract GSE ID for folder name (assuming it starts with GSE followed by numbers)
    match = re.search(r"(GSE\d+)", filename)
    if match:
        gse_id = match.group(1)
    else:
        # Fallback if no GSE ID found, use filename without extension
        gse_id = os.path.splitext(filename)[0]
    
    output_dir = os.path.join(output_base, gse_id)
    os.makedirs(output_dir, exist_ok=True)
    
    output_file_path = os.path.join(output_dir, filename)
    
    print(f"Processing {filename}...")
    try:
        df = pd.read_excel(full_path)
        
        # Check if columns exist
        if 'Gene' in df.columns and 'logFC' in df.columns:
            df_subset = df[['Gene', 'logFC']]
            df_subset.to_excel(output_file_path, index=False)
            print(f"  Saved to: {output_file_path}")
            
            csv_output_path = os.path.splitext(output_file_path)[0] + ".csv"
            df_subset.to_csv(csv_output_path, index=False)
            print(f"  Saved to: {csv_output_path}")
        else:
            print(f"  Error: Columns 'Gene' and/or 'logFC' not found in {filename}")
            
    except Exception as e:
        print(f"  Error processing {filename}: {e}")

print("Processing complete.")
