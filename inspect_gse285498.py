import pandas as pd
import os
import sys

# Set encoding to utf-8 for output
sys.stdout.reconfigure(encoding='utf-8')

files = [
    r"results/GSE285498/degs/CBD_vs_Eto_DEGs_FDR0.05_logFC1.tsv",
    r"results/GSE285498/degs/CBD_vs_Mock_DEGs_FDR0.05_logFC1.tsv",
    r"results/GSE285498/degs/Eto_vs_Mock_DEGs_FDR0.05_logFC1.tsv",
    r"results/GSE285498/degs/Merge_vs_Mock_DEGs_FDR0.05_logFC1.tsv"
]

base_path = r"c:\Users\nicol\Desktop\Proyectos\differential-expression-uterine-cancer"
output_file = "inspection_results_gse285498.txt"

with open(output_file, "w", encoding="utf-8") as f:
    for relative_path in files:
        full_path = os.path.join(base_path, relative_path)
        f.write(f"--- FILE: {relative_path} ---\n")
        try:
            # Read TSV
            df = pd.read_csv(full_path, sep='\t', nrows=2)
            f.write(f"COLUMNS: {list(df.columns)}\n")
            f.write(f"FIRST ROW: {df.iloc[0].to_dict()}\n")
            f.write("-" * 30 + "\n")
        except Exception as e:
            f.write(f"Error reading {relative_path}: {e}\n")

print(f"Done writing to {output_file}")
