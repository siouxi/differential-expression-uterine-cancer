import pandas as pd
import os
import sys

# Set encoding to utf-8 for output
sys.stdout.reconfigure(encoding='utf-8')

files = [
    r"results/GSE223827/deg/DESeq2_CBD_vs_DMSO.tsv",
    r"results/GSE223827/deg/DESeq2_CBD-THC_vs_DMSO.tsv",
    r"results/GSE223827/deg/DESeq2_DMSO-EtOH_vs_DMSO.tsv",
    r"results/GSE223827/deg/DESeq2_EtOH_vs_DMSO.tsv",
    r"results/GSE223827/deg/DESeq2_THC_vs_DMSO.tsv"
]

base_path = r"c:\Users\nicol\Desktop\Proyectos\differential-expression-uterine-cancer"
output_file = "inspection_results_gse223827.txt"

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
