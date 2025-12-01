import pandas as pd

file_path = r"c:\Users\nicol\Desktop\Proyectos\differential-expression-uterine-cancer\standard\GSE21656\GSE21656_DEGS_significativos_standard.csv"

# Read the file
df = pd.read_csv(file_path)

# Convert Entrez_ID to integer (removing .0), handling NaN values
df['Entrez_ID'] = df['Entrez_ID'].apply(lambda x: str(int(x)) if pd.notna(x) else '')

# Save the file
df.to_csv(file_path, index=False)

print(f"✓ Fixed Entrez_ID column in {file_path}")
print(f"  Total rows: {len(df)}")
print(f"  Sample values:")
print(df[['Entrez_ID', 'Gene_Symbol']].head(10))
