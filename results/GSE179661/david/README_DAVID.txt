DAVID Functional Annotation Tool - Input Files
==============================================

This directory contains gene lists prepared for DAVID analysis.

FILE FORMATS:
-------------

1. *_EntrezID_list.txt:
   - Simple list of Entrez Gene IDs (one per line)
   - Use this in DAVID: Select 'ENTREZ_GENE_ID' as identifier

2. *_GeneSymbol_list.txt:
   - Simple list of Official Gene Symbols (one per line)
   - Use this in DAVID: Select 'OFFICIAL_GENE_SYMBOL' as identifier

3. *_EnsemblID_list.txt:
   - Simple list of Ensembl Gene IDs (one per line)
   - Use this in DAVID: Select 'ENSEMBL_GENE_ID' as identifier

4. *_DAVID_FC_format.txt:
   - Tab-separated file: EntrezID\tFoldChange
   - Use for Functional Annotation Clustering with fold change weighting

5. *_complete_table.tsv:
   - Complete table with all annotations for reference

HOW TO USE DAVID:
-----------------

1. Go to: https://david.ncifcrf.gov/
2. Click 'Start Analysis'
3. Upload your gene list file (choose appropriate format)
4. Select the correct identifier type:
   - For *_EntrezID_list.txt: Select 'ENTREZ_GENE_ID'
   - For *_GeneSymbol_list.txt: Select 'OFFICIAL_GENE_SYMBOL'
   - For *_EnsemblID_list.txt: Select 'ENSEMBL_GENE_ID'
5. Select species: Homo sapiens
6. Click 'Submit List'
7. Select annotation categories:
   - GOTERM_BP_DIRECT (GO Biological Process)
   - GOTERM_CC_DIRECT (GO Cellular Component)
   - GOTERM_MF_DIRECT (GO Molecular Function)
   - KEGG_PATHWAY
   - BIOCARTA
   - REACTOME_PATHWAY
   - etc.
8. Click 'Functional Annotation Clustering' for grouped results
9. Click 'Functional Annotation Chart' for detailed results

RECOMMENDED SETTINGS:
--------------------
- EASE Score (p-value threshold): 0.1 (default)
- Count: 2 (minimum genes per term)
- Classification Stringency: Medium

FILES AVAILABLE:
----------------
- DAVID_all_DEGs_*: All significant DEGs (FDR ≤ 0.05, |log2FC| ≥ 1)
- DAVID_upregulated_*: Up-regulated DEGs (CBD > Control)
- DAVID_downregulated_*: Down-regulated DEGs (CBD < Control)

NOTE:
-----
DAVID provides complementary analysis to clusterProfiler.
Compare results from both tools for comprehensive interpretation.
