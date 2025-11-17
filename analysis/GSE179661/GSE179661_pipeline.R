required_pkgs <- c(
  "tidyverse", "here", "jsonlite", "glue", "fs",
  "BiocManager", "GEOquery", "DESeq2", "edgeR",
  "SummarizedExperiment", "apeglm", "EnhancedVolcano",
  "pheatmap", "ComplexHeatmap", "RColorBrewer", "plotly",
  "clusterProfiler", "enrichplot", "msigdbr", "fgsea",
  "org.Hs.eg.db", "AnnotationDbi", "biomaRt", "GSVA",
  "viridis", "VennDiagram", "ggrepel"
)

install_and_load <- function(pkgs) {
  new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(new_pkgs)) {
    message("Installing missing packages: ", paste(new_pkgs, collapse = ", "))
    BiocManager::install(new_pkgs, ask = FALSE, update = FALSE)
  }
  suppressPackageStartupMessages(
    invisible(lapply(pkgs, library, character.only = TRUE))
  )
}

install_and_load(required_pkgs)


geo_id <- "GSE179661"

script_dir <- tryCatch({
  if (require("rstudioapi", quietly = TRUE)) {
    path <- rstudioapi::getActiveDocumentContext()$path
    if (length(path) > 0 && path != "") {
      dirname(path)
    } else {
      NULL
    }
  } else {
    NULL
  }
}, error = function(e) NULL)


base_dir <- dirname(dirname(script_dir))

paths <- list(
  data_dir = file.path(base_dir, "DATA", geo_id),
  results_dir = file.path(base_dir, "results", geo_id),
  qc_dir = file.path(base_dir, "results", geo_id, "qc"),
  deseq_dir = file.path(base_dir, "results", geo_id, "deseq2"),
  enrich_dir = file.path(base_dir, "results", geo_id, "enrichment")
)

fs::dir_create(unlist(paths))

cat("Project paths configured:\n")
cat("  Base directory:", base_dir, "\n")
cat("  Data directory:", paths$data_dir, "\n")
cat("  Results directory:", paths$results_dir, "\n\n")


# 3. Load and process phenotype data ----------------------------------------------------------------
cat("Loading phenotype data...\n")
pheno_file <- file.path(paths$data_dir, paste0(geo_id, "_pheno_data.csv"))

if (!file.exists(pheno_file)) {
  stop(paste("Phenotype file not found:", pheno_file))
}

pheno_raw <- readr::read_csv(pheno_file, show_col_types = FALSE)

# Extract and process sample metadata
# From the pheno_data, we have: title, geo_accession, treatment:ch1, source_name_ch1
sample_sheet <- pheno_raw %>%
  dplyr::select(
    sample_id = title,
    geo_accession,
    treatment = `treatment:ch1`,
    cell_line = source_name_ch1
  ) %>%
  dplyr::mutate(
    # Extract condition from treatment column
    condition = dplyr::case_when(
      grepl("none treated", treatment, ignore.case = TRUE) ~ "Control",
      grepl("CBD", treatment, ignore.case = TRUE) ~ "CBD",
      TRUE ~ "Unknown"
    ),
    # Extract cell line type
    cell_line = dplyr::case_when(
      grepl("MHCC97H", cell_line, ignore.case = TRUE) ~ "MHCC97H",
      grepl("HepG2", cell_line, ignore.case = TRUE) ~ "HepG2",
      TRUE ~ cell_line
    ),
    # Extract replicate number
    replicate = stringr::str_extract(sample_id, "rep\\d+"),
    replicate = as.numeric(stringr::str_replace(replicate, "rep", "")),
    # Create sample IDs that match count matrix column names
    # The count matrix likely uses GSM IDs or sample names from title
    sample_id_clean = geo_accession
  ) %>%
  dplyr::mutate(
    condition = factor(condition, levels = c("Control", "CBD")),
    cell_line = factor(cell_line)
  )

cat("  Number of samples:", nrow(sample_sheet), "\n")
cat("  Conditions:", paste(levels(sample_sheet$condition), collapse = ", "), "\n")
cat("  Cell lines:", paste(levels(sample_sheet$cell_line), collapse = ", "), "\n\n")


# 4. Load raw gene counts matrix ----------------------------------------------------------------
cat("Loading raw gene counts matrix...\n")
counts_file <- file.path(paths$data_dir, paste0(geo_id, "_Raw_gene_counts_matrix.txt.gz"))

if (!file.exists(counts_file)) {
  stop(paste("Counts file not found:", counts_file))
}

# Read gzipped counts file
# The file format may vary - try different approaches
counts_df <- tryCatch({
  readr::read_tsv(counts_file, comment = "#", show_col_types = FALSE)
}, error = function(e) {
  # Alternative: read as table
  read.table(gzfile(counts_file), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
})

# Identify gene ID column (could be Geneid, gene_id, or first column)
gene_col <- NULL
if ("gene_id" %in% colnames(counts_df)) {
  gene_col <- "gene_id"
} else if ("gene_id" %in% colnames(counts_df)) {
  gene_col <- "gene_id"
} else {
  gene_col <- colnames(counts_df)[1]
}

# Prepare count matrix
# Remove annotation columns if present
annotation_cols <- c("Chr", "Start", "End", "Length", "Strand", "Geneid", "gene_id")
cols_to_remove <- intersect(annotation_cols, colnames(counts_df))
cols_to_remove <- setdiff(cols_to_remove, gene_col)

if (length(cols_to_remove) > 0) {
  counts_mat <- counts_df %>%
    dplyr::select(-dplyr::all_of(cols_to_remove)) %>%
    tibble::column_to_rownames(gene_col) %>%
    as.matrix()
} else {
  counts_mat <- counts_df %>%
    tibble::column_to_rownames(gene_col) %>%
    as.matrix()
}

# Ensure counts are numeric
counts_mat <- apply(counts_mat, 2, as.numeric)
rownames(counts_mat) <- rownames(counts_df)

cat("  Count matrix dimensions:", dim(counts_mat), "\n")
cat("  Sample names in matrix:", paste(head(colnames(counts_mat), 5), collapse = ", "), "...\n\n")

# Match sample names between count matrix and phenotype data
cat("Matching sample names...\n")
sample_cols <- colnames(counts_mat)

# Diagnostic: show what we have
cat("  Count matrix columns (", length(sample_cols), "):", paste(sample_cols, collapse = ", "), "\n")
cat("  Phenotype sample_id (", length(sample_sheet$sample_id), "):", paste(sample_sheet$sample_id, collapse = ", "), "\n")
cat("  Phenotype geo_accession (", length(sample_sheet$geo_accession), "):", paste(sample_sheet$geo_accession, collapse = ", "), "\n")

# Helper function to normalize sample names
# Converts "H_0 rep1" -> "h_0_1", "HE_40 rep2" -> "he_40_2", etc.
normalize_sample_name <- function(name) {
  # Convert to lowercase
  name_lower <- tolower(name)
  # Replace " rep" or "rep" with "_"
  name_lower <- gsub("\\s*rep\\s*", "_", name_lower)
  # Remove extra spaces and normalize underscores
  name_lower <- gsub("\\s+", "_", name_lower)
  name_lower <- gsub("_+", "_", name_lower)
  # Remove trailing/leading underscores
  name_lower <- gsub("^_|_$", "", name_lower)
  return(name_lower)
}

# Create normalized versions for matching
sample_sheet <- sample_sheet %>%
  dplyr::mutate(
    sample_id_normalized = normalize_sample_name(sample_id),
    geo_accession_normalized = tolower(geo_accession)
  )

sample_cols_normalized <- tolower(sample_cols)

# Strategy 1: Direct match with normalized sample_id
cat("  Strategy 1: Matching by normalized sample_id...\n")
matched_indices <- match(sample_cols_normalized, sample_sheet$sample_id_normalized)
matched_samples <- sample_cols[!is.na(matched_indices)]

if (length(matched_samples) > 0) {
  cat("    ✓ Matched", length(matched_samples), "samples by normalized sample_id\n")
  sample_sheet$sample_id_clean[matched_indices[!is.na(matched_indices)]] <- sample_cols[!is.na(matched_indices)]
} else {
  # Strategy 2: Direct match with geo_accession
  cat("  Strategy 2: Matching by geo_accession...\n")
  matched_indices <- match(sample_cols_normalized, sample_sheet$geo_accession_normalized)
  matched_samples <- sample_cols[!is.na(matched_indices)]
  
  if (length(matched_samples) > 0) {
    cat("    ✓ Matched", length(matched_samples), "samples by geo_accession\n")
    sample_sheet$sample_id_clean[matched_indices[!is.na(matched_indices)]] <- sample_cols[!is.na(matched_indices)]
  } else {
    # Strategy 3: Partial matching - check if count matrix names contain parts of sample_id
    cat("  Strategy 3: Partial matching by pattern...\n")
    matched_samples <- character(0)
    matched_pheno_indices <- integer(0)
    
    for (i in seq_along(sample_cols)) {
      col_name <- sample_cols_normalized[i]
      # Try to match with normalized sample_id
      for (j in seq_along(sample_sheet$sample_id_normalized)) {
        pheno_name <- sample_sheet$sample_id_normalized[j]
        # Check if one contains the other (flexible matching)
        if (grepl(gsub("_", ".", col_name), pheno_name, ignore.case = TRUE) ||
            grepl(gsub("_", ".", pheno_name), col_name, ignore.case = TRUE) ||
            # Check for pattern like "h_0" matching "H_0 rep1"
            grepl(paste0("^", gsub("_\\d+$", "", col_name)), pheno_name, ignore.case = TRUE)) {
          if (!(j %in% matched_pheno_indices)) {
            matched_samples <- c(matched_samples, sample_cols[i])
            matched_pheno_indices <- c(matched_pheno_indices, j)
            sample_sheet$sample_id_clean[j] <- sample_cols[i]
            break
          }
        }
      }
    }
    
    if (length(matched_samples) > 0) {
      cat("    ✓ Matched", length(matched_samples), "samples by partial pattern\n")
    } else {
      # Strategy 4: Try matching by extracting key patterns (H_0, H_40, HE_0, HE_40)
      cat("  Strategy 4: Matching by extracted patterns...\n")
      matched_samples <- character(0)
      matched_pheno_indices <- integer(0)
      
      # Extract pattern from count matrix names (e.g., "h_0_1" -> "h_0")
      extract_pattern <- function(name) {
        # Remove trailing numbers and underscores
        pattern <- gsub("_\\d+$", "", tolower(name))
        return(pattern)
      }
      
      # Extract pattern from phenotype names (e.g., "H_0 rep1" -> "h_0")
      extract_pheno_pattern <- function(name) {
        pattern <- tolower(name)
        pattern <- gsub("\\s*rep\\d+", "", pattern)
        pattern <- gsub("\\s+", "_", pattern)
        return(pattern)
      }
      
      count_patterns <- sapply(sample_cols, extract_pattern)
      pheno_patterns <- sapply(sample_sheet$sample_id, extract_pheno_pattern)
      
      for (i in seq_along(count_patterns)) {
        match_idx <- which(pheno_patterns == count_patterns[i])
        if (length(match_idx) > 0 && !(match_idx[1] %in% matched_pheno_indices)) {
          matched_samples <- c(matched_samples, sample_cols[i])
          matched_pheno_indices <- c(matched_pheno_indices, match_idx[1])
          sample_sheet$sample_id_clean[match_idx[1]] <- sample_cols[i]
        }
      }
      
      if (length(matched_samples) > 0) {
        cat("    ✓ Matched", length(matched_samples), "samples by extracted patterns\n")
      }
    }
  }
}

# Final check: if still no matches, use sequential assignment as fallback
if (length(matched_samples) == 0 || length(matched_samples) < nrow(sample_sheet)) {
  cat("  Strategy 5: Fallback - sequential assignment\n")
  unmatched_pheno <- which(is.na(sample_sheet$sample_id_clean) | sample_sheet$sample_id_clean == sample_sheet$geo_accession)
  unmatched_counts <- setdiff(sample_cols, matched_samples)
  
  if (length(unmatched_pheno) > 0 && length(unmatched_counts) > 0) {
    n_assign <- min(length(unmatched_pheno), length(unmatched_counts))
    for (i in 1:n_assign) {
      sample_sheet$sample_id_clean[unmatched_pheno[i]] <- unmatched_counts[i]
      matched_samples <- c(matched_samples, unmatched_counts[i])
    }
    cat("    ⚠ Assigned", n_assign, "samples sequentially (verify manually!)\n")
  }
}

# Ensure all samples have a clean ID
if (any(is.na(sample_sheet$sample_id_clean))) {
  sample_sheet$sample_id_clean[is.na(sample_sheet$sample_id_clean)] <- sample_sheet$geo_accession[is.na(sample_sheet$sample_id_clean)]
}

# Show matching results
cat("\n  Matching results:\n")
matching_df <- data.frame(
  Phenotype_ID = sample_sheet$sample_id,
  Count_Matrix_ID = sample_sheet$sample_id_clean,
  Condition = sample_sheet$condition,
  Cell_Line = sample_sheet$cell_line,
  Matched = sample_sheet$sample_id_clean %in% sample_cols
)
print(matching_df)

# Filter and reorder
matched_samples <- unique(matched_samples[matched_samples %in% sample_cols])
sample_sheet <- sample_sheet %>%
  dplyr::filter(sample_id_clean %in% matched_samples | sample_id_clean %in% sample_cols) %>%
  dplyr::arrange(match(sample_id_clean, sample_cols))

# Ensure counts matrix has the matched samples in correct order
available_samples <- intersect(sample_sheet$sample_id_clean, sample_cols)
counts_mat <- counts_mat[, available_samples, drop = FALSE]
sample_sheet <- sample_sheet %>%
  dplyr::filter(sample_id_clean %in% available_samples) %>%
  dplyr::arrange(match(sample_id_clean, colnames(counts_mat)))

# Final verification
if (ncol(counts_mat) != nrow(sample_sheet)) {
  warning("Mismatch: count matrix has ", ncol(counts_mat), " columns but sample_sheet has ", nrow(sample_sheet), " rows")
}

cat("\n  ✓ Final matched samples:", length(available_samples), "\n")
cat("  ✓ Count matrix columns:", ncol(counts_mat), "\n")
cat("  ✓ Sample sheet rows:", nrow(sample_sheet), "\n\n")


# 5. Create DESeq2 dataset and pre-filtering ------------------------------------------------
cat("Creating DESeq2 dataset...\n")

# Create sample metadata for DESeq2
col_data <- sample_sheet %>%
  dplyr::select(sample_id = sample_id_clean, condition, cell_line, replicate) %>%
  as.data.frame()
rownames(col_data) <- col_data$sample_id

# Create DESeq2 dataset
# Design: account for cell line and condition
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = col_data,
  design = ~ cell_line + condition
)

# Pre-filtering: keep genes with sufficient expression
# IMPORTANT FOR CANCER STUDIES: Use less aggressive filtering to preserve 
# sub-expressed genes that may be biologically relevant (e.g., tumor suppressors)
# Strategy: at least 5 counts in at least 25% of samples (less strict than typical)
# This preserves genes that might be downregulated in cancer while removing 
# genes with essentially no signal (technical zeros)
n_samples <- ncol(dds)
min_samples_with_counts <- max(2, ceiling(n_samples * 0.25))  # At least 25% or minimum 2 (less strict)
count_threshold <- 5  # Lower threshold (5 instead of 10) to preserve sub-expressed genes

keep <- rowSums(DESeq2::counts(dds) >= count_threshold) >= min_samples_with_counts
dds <- dds[keep, ]

cat("  Pre-filtering criteria (adjusted for cancer study - preserves sub-expressed genes):\n")
cat("    - Count threshold:", count_threshold, "counts per sample (lowered from 10 to preserve downregulated genes)\n")
cat("    - Minimum samples:", min_samples_with_counts, "of", n_samples, "samples (", 
    round(min_samples_with_counts/n_samples*100, 1), "% - less strict to capture sub-expressed genes)\n")
cat("    - Rationale: In cancer, sub-expressed genes (e.g., tumor suppressors) are biologically relevant\n")
cat("  Genes after filtering:", nrow(dds), "\n")
cat("  Samples:", ncol(dds), "\n\n")


# 6. Quality control and normalization --------------------------------------------------------
cat("Performing quality control and normalization...\n")
cat("  IMPORTANT: Multiple normalizations are calculated, but only DESeq2 size factors\n")
cat("  are used in the differential expression analysis. Others are for QC only.\n\n")

# Size factors (DESeq2)
# Estimates: normalization factors for each sample (accounts for library size differences)
# 
# MEDIDAS ESTADÍSTICAS QUE USA:
# - Median ratio method: Calcula la mediana de las razones (ratios) entre cada muestra y una 
#   muestra de referencia (pseudo-reference) construida como la mediana geométrica de todos los genes
# - Fórmula: sizeFactor_i = median( counts_ij / geometric_mean(counts_j) ) para todos los genes j
# - Excluye genes con conteos extremos (outliers) usando el método de mediana absoluta (MAD)
# - Resultado: Factor de normalización por muestra (típicamente entre 0.5 y 2.0)
# 
# Access with: sizeFactors(dds) or colData(dds)$sizeFactor
dds <- DESeq2::estimateSizeFactors(dds)
cat("  DESeq2 size factors:\n")
print(sizeFactors(dds))
cat("\n")

# Dispersion estimation (DESeq2)
# Estimates: gene-wise dispersion, fitted dispersion, and prior dispersion
# 
# MEDIDAS ESTADÍSTICAS QUE USA:
# 1. Gene-wise dispersion (dispGeneEst):
#    - Método de momentos: Var = μ + α*μ² donde α es la dispersión
#    - Calcula dispersión inicial para cada gen usando la fórmula: α = (Var - μ) / μ²
#    - Usa estimación de máxima verosimilitud (MLE) para genes con suficientes conteos
# 
# 2. Fitted dispersion (dispFit):
#    - Regresión local (loess): Ajusta una curva suave de dispersión vs. expresión media
#    - Modelo: log(α) ~ f(log(μ)) donde μ es la expresión media del gen
#    - Predice dispersión esperada basada en el nivel de expresión
# 
# 3. Prior dispersion (dispersion):
#    - Shrinkage hacia el valor ajustado: Combina dispersión gen-específica y ajustada
#    - Fórmula: α_final = w*α_gene + (1-w)*α_fitted, donde w es un peso basado en la precisión
#    - Usa distribución a priori (prior) para estabilizar estimaciones de genes con pocos datos
# 
# Access with: dispersions(dds), mcols(dds)$dispGeneEst, mcols(dds)$dispFit, mcols(dds)$dispersion
dds <- DESeq2::estimateDispersions(dds)
cat("  Dispersion estimation completed\n")
cat("  Gene-wise dispersion range:", 
    round(range(dispersions(dds), na.rm = TRUE), 4), "\n\n")

# Variance stabilizing transformation (VST)
# Transforms counts to log2-like scale with stabilized variance
# 
# MEDIDAS ESTADÍSTICAS QUE USA:
# - Transformación basada en la dispersión estimada: Usa los valores de dispersión (α) calculados
#   previamente para estabilizar la varianza
# - Fórmula aproximada: VST(x) = ∫[1/sqrt(μ + α*μ²)] dμ desde 0 hasta x
#   donde μ es la expresión media y α es la dispersión
# - Para valores grandes: Se aproxima a log2(x) (comportamiento logarítmico)
# - Para valores pequeños: Mantiene la escala lineal (evita problemas con ceros)
# - Resultado: Datos transformados donde la varianza es aproximadamente constante (homocedasticidad)
# - Ventaja sobre log2: No necesita pseudocounts y la varianza es más estable
# 
# Parámetro 'blind = FALSE': Usa la información de las condiciones experimentales en la 
# transformación (mejor para análisis exploratorio cuando ya conocemos los grupos)
# 
# Access with: assay(vsd) - returns transformed counts
vsd <- DESeq2::vst(dds, blind = FALSE)
cat("  VST transformation completed\n")
cat("  VST matrix dimensions:", dim(assay(vsd)), "\n\n")

# TMM normalization (edgeR) for comparison
# Estimates: normalization factors using Trimmed Mean of M-values method
# 
# MEDIDAS ESTADÍSTICAS QUE USA:
# 1. M-values (log fold changes): M = log2(counts_sample / counts_reference)
#    - Calcula el log2 de la razón entre cada muestra y una muestra de referencia
# 
# 2. A-values (mean expression): A = (log2(counts_sample) + log2(counts_reference)) / 2
#    - Calcula el promedio de expresión en escala logarítmica
# 
# 3. Trimmed Mean (media recortada):
#    - Elimina el 30% superior e inferior de los M-values (por defecto trim=0.3)
#    - Calcula la media de los M-values restantes (40% central)
#    - Esto elimina genes diferencialmente expresados y genes con baja expresión
# 
# 4. Normalization factor: 
#    - Factor = 2^(trimmed_mean_M)
#    - Ajusta el tamaño de biblioteca de cada muestra relativo a la referencia
# 
# 5. Reference sample:
#    - Usa la muestra con el percentil 75 más cercano a la mediana de los conteos totales
#    - O puede usar una muestra específica como referencia
# 
# Ventaja: Robusto a genes diferencialmente expresados (al eliminar extremos)
# 
# Access with: dds_edgeR$samples$norm.factors
dds_edgeR <- edgeR::DGEList(counts = counts_mat[rownames(dds), colnames(dds)])
dds_edgeR <- edgeR::calcNormFactors(dds_edgeR, method = "TMM")
cat("  edgeR TMM normalization factors:\n")
print(dds_edgeR$samples$norm.factors)
cat("\n")

# Calculate CPM for QC
# Calculates: Counts Per Million (normalized by library size and TMM factors)
# 
# MEDIDAS ESTADÍSTICAS QUE USA:
# - Fórmula: CPM = (counts / library_size) * 1,000,000
#   donde library_size = suma de conteos de todos los genes en la muestra
# 
# - Con factores TMM: CPM = (counts / (library_size * TMM_factor)) * 1,000,000
#   - Ajusta por diferencias en tamaño de biblioteca Y por sesgos de composición (TMM)
# 
# - Interpretación: Número de conteos por millón de reads mapeados
#   - Permite comparar expresión entre muestras con diferentes tamaños de biblioteca
#   - Ejemplo: 10 CPM significa 10 conteos por cada millón de reads totales
# 
# - Uso: Control de calidad, visualización, filtrado (ej: genes con CPM < 1 en todas las muestras)
# 
# Parámetro 'log = FALSE': Devuelve CPM en escala lineal (no logarítmica)
# 
# Access with: cpm_mat - matrix of CPM values
cpm_mat <- edgeR::cpm(dds_edgeR, log = FALSE)
cat("  CPM matrix calculated\n")
cat("  CPM matrix dimensions:", dim(cpm_mat), "\n")
cat("  CPM range:", round(range(cpm_mat, na.rm = TRUE), 2), "\n\n")

# Save normalization summary for inspection
normalization_summary <- data.frame(
  Sample = colnames(dds),
  DESeq2_SizeFactor = sizeFactors(dds),
  edgeR_TMM_Factor = dds_edgeR$samples$norm.factors,
  Library_Size = colSums(counts(dds)),
  Mean_CPM = colMeans(cpm_mat)
)
readr::write_tsv(normalization_summary, 
                 file.path(paths$qc_dir, "normalization_factors_summary.tsv"))
cat("  Normalization summary saved to:", 
    file.path(paths$qc_dir, "normalization_factors_summary.tsv"), "\n\n")

# QC plots directory
fs::dir_create(paths$qc_dir)

# 6.1. Sample count distribution (boxplot) with colors
# Prepare data for boxplot with colors based on condition and cell line
# Ensure counts_mat columns match col_data rows
counts_log2 <- log2(counts_mat[, rownames(col_data), drop = FALSE] + 1)

# Create color scheme based on condition and cell line
# Colors: Control=Blue shades, CBD=Purple/Pink shades
# Intensity: MHCC97H=Lighter, HepG2=Darker
boxplot_colors <- character(ncol(counts_log2))
sample_labels <- character(ncol(counts_log2))

for (i in 1:ncol(counts_log2)) {
  sample_name <- colnames(counts_log2)[i]
  sample_idx <- which(rownames(col_data) == sample_name)
  
  if (length(sample_idx) > 0) {
    condition <- as.character(col_data$condition[sample_idx])
    cell_line <- as.character(col_data$cell_line[sample_idx])
    
    # Assign colors based on condition and cell line combination
    if (condition == "Control" && cell_line == "MHCC97H") {
      boxplot_colors[i] <- "#6BAED6"  # Light blue
      sample_labels[i] <- paste0(sample_name, "\n(Control-MHCC97H)")
    } else if (condition == "Control" && cell_line == "HepG2") {
      boxplot_colors[i] <- "#2E86AB"  # Dark blue
      sample_labels[i] <- paste0(sample_name, "\n(Control-HepG2)")
    } else if (condition == "CBD" && cell_line == "MHCC97H") {
      boxplot_colors[i] <- "#D67BA8"  # Light purple/pink
      sample_labels[i] <- paste0(sample_name, "\n(CBD-MHCC97H)")
    } else if (condition == "CBD" && cell_line == "HepG2") {
      boxplot_colors[i] <- "#A23B72"  # Dark purple
      sample_labels[i] <- paste0(sample_name, "\n(CBD-HepG2)")
    } else {
      boxplot_colors[i] <- "gray"
      sample_labels[i] <- sample_name
    }
  } else {
    boxplot_colors[i] <- "gray"
    sample_labels[i] <- sample_name
  }
}

pdf(file.path(paths$qc_dir, "count_distribution_boxplot.pdf"), width = 14, height = 7)
par(mar = c(10, 4, 4, 6))  # Increase margins for sample names and legend
boxplot(counts_log2, 
        main = "Raw counts distribution (log2) by Condition and Cell Line",
        xlab = "", 
        ylab = "Log2(counts + 1)",
        las = 2, 
        cex.axis = 0.65,
        col = boxplot_colors,
        border = "black",
        names = colnames(counts_log2))

# Add legend
legend("topright", 
       legend = c("Control - MHCC97H", "Control - HepG2", "CBD - MHCC97H", "CBD - HepG2"),
       fill = c("#6BAED6", "#2E86AB", "#D67BA8", "#A23B72"),
       cex = 0.8,
       title = "Condition - Cell Line")

# Add grid for better readability
grid(nx = NA, ny = NULL, col = "gray90", lty = "dotted")

dev.off()

# Also create a ggplot version for better aesthetics
# Match col_data to counts_log2 column order
col_data_ordered <- col_data[colnames(counts_log2), , drop = FALSE]
boxplot_df <- data.frame(
  Sample = rep(colnames(counts_log2), each = nrow(counts_log2)),
  Log2Counts = as.vector(counts_log2),
  Condition = rep(col_data_ordered$condition, each = nrow(counts_log2)),
  CellLine = rep(col_data_ordered$cell_line, each = nrow(counts_log2))
)

boxplot_gg <- ggplot2::ggplot(boxplot_df, ggplot2::aes(x = Sample, y = Log2Counts, fill = Condition, alpha = CellLine)) +
  ggplot2::geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  ggplot2::scale_fill_manual(values = c("Control" = "#2E86AB", "CBD" = "#A23B72"),
                            name = "Condition") +
  ggplot2::scale_alpha_manual(values = c("MHCC97H" = 0.7, "HepG2" = 1.0),
                             name = "Cell Line") +
  ggplot2::labs(
    title = "Raw counts distribution (log2)",
    subtitle = "Colors: Condition (Control=Blue, CBD=Purple), Transparency: Cell Line (MHCC97H=Lighter, HepG2=Darker)",
    x = "Samples",
    y = "Log2(counts + 1)"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
    plot.subtitle = ggplot2::element_text(size = 9),
    legend.position = "right"
  )

ggplot2::ggsave(file.path(paths$qc_dir, "count_distribution_boxplot_ggplot.pdf"), 
                boxplot_gg, width = 12, height = 7)

# Calculate and report mean values for boxplot data
cat("Boxplot data statistics (log2(counts + 1)):\n")
cat("  ===========================================\n")
boxplot_stats <- data.frame(
  Sample = colnames(counts_log2),
  Mean = colMeans(counts_log2),
  Median = apply(counts_log2, 2, median),
  SD = apply(counts_log2, 2, sd),
  Min = apply(counts_log2, 2, min),
  Max = apply(counts_log2, 2, max),
  Q1 = apply(counts_log2, 2, quantile, probs = 0.25),
  Q3 = apply(counts_log2, 2, quantile, probs = 0.75)
)

# Add condition and cell line info
boxplot_stats <- boxplot_stats %>%
  dplyr::left_join(
    col_data %>% 
      tibble::rownames_to_column("Sample") %>%
      dplyr::select(Sample, Condition = condition, CellLine = cell_line),
    by = "Sample"
  ) %>%
  dplyr::select(Sample, Condition, CellLine, Mean, Median, SD, Min, Max, Q1, Q3)

cat("  Summary statistics by sample:\n")
print(boxplot_stats, row.names = FALSE)
cat("\n")

# Overall statistics
cat("  Overall statistics (all samples combined):\n")
cat("    Overall mean:", round(mean(counts_log2), 4), "\n")
cat("    Overall median:", round(median(counts_log2), 4), "\n")
cat("    Overall SD:", round(sd(counts_log2), 4), "\n")
cat("\n")

# Statistics by condition
cat("  Statistics by condition:\n")
for (cond in levels(col_data$condition)) {
  cond_samples <- rownames(col_data)[col_data$condition == cond]
  cond_data <- counts_log2[, cond_samples, drop = FALSE]
  cat("    ", cond, ":\n")
  cat("      Mean:", round(mean(cond_data), 4), "\n")
  cat("      Median:", round(median(cond_data), 4), "\n")
  cat("      SD:", round(sd(cond_data), 4), "\n")
}
cat("\n")

# Statistics by cell line
cat("  Statistics by cell line:\n")
for (cl in levels(col_data$cell_line)) {
  cl_samples <- rownames(col_data)[col_data$cell_line == cl]
  cl_data <- counts_log2[, cl_samples, drop = FALSE]
  cat("    ", cl, ":\n")
  cat("      Mean:", round(mean(cl_data), 4), "\n")
  cat("      Median:", round(median(cl_data), 4), "\n")
  cat("      SD:", round(sd(cl_data), 4), "\n")
}
cat("\n")

# Save statistics to file
readr::write_tsv(boxplot_stats, 
                 file.path(paths$qc_dir, "boxplot_statistics.tsv"))
cat("  Boxplot statistics saved to:", 
    file.path(paths$qc_dir, "boxplot_statistics.tsv"), "\n\n")

# 6.2. PCA plot
pca_data <- DESeq2::plotPCA(vsd, intgroup = c("condition", "cell_line"), returnData = TRUE)
pca_plot <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2, 
                                                     color = condition, shape = cell_line)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::labs(title = "PCA: Condition and Cell Line",
                x = paste0("PC1: ", round(attr(pca_data, "percentVar")[1] * 100, 1), "% variance"),
                y = paste0("PC2: ", round(attr(pca_data, "percentVar")[2] * 100, 1), "% variance")) +
  ggplot2::theme_bw()

ggplot2::ggsave(file.path(paths$qc_dir, "PCA_condition_cell_line.pdf"), 
                pca_plot, width = 8, height = 6)

# 6.3. Sample distance heatmap
sample_dists <- dist(t(SummarizedExperiment::assay(vsd)))
sample_dist_mat <- as.matrix(sample_dists)
rownames(sample_dist_mat) <- colnames(vsd)
colnames(sample_dist_mat) <- colnames(vsd)

pdf(file.path(paths$qc_dir, "sample_distance_heatmap.pdf"), width = 8, height = 7)
pheatmap::pheatmap(
  mat = sample_dist_mat,
  color = RColorBrewer::brewer.pal(9, "Blues"),
  clustering_method = "complete",
  main = "Sample-to-sample distances (VST)",
  annotation_col = col_data[, c("condition", "cell_line"), drop = FALSE]
)
dev.off()

# 6.4. Density plot of normalized counts
pdf(file.path(paths$qc_dir, "normalized_counts_density.pdf"), width = 8, height = 6)
plotDensity <- function(mat, main) {
  plot(density(log2(mat[, 1] + 1)), main = main, xlab = "Log2(counts + 1)", 
       ylab = "Density", col = 1, lwd = 2)
  for (i in 2:ncol(mat)) {
    lines(density(log2(mat[, i] + 1)), col = i, lwd = 2)
  }
  legend("topright", legend = colnames(mat), col = 1:ncol(mat), lty = 1, cex = 0.6)
}
plotDensity(SummarizedExperiment::assay(vsd), "VST normalized counts")
dev.off()

cat("  QC plots saved to:", paths$qc_dir, "\n\n")

# 6.5. Generate QC plots AFTER pre-filtering for comparison
cat("Generating QC plots after pre-filtering (for comparison)...\n")

# Boxplot after filtering
counts_filtered_log2 <- log2(DESeq2::counts(dds) + 1)
col_data_filtered <- col_data[colnames(dds), , drop = FALSE]

# Create colors for filtered data
boxplot_colors_filtered <- character(ncol(counts_filtered_log2))
for (i in 1:ncol(counts_filtered_log2)) {
  sample_name <- colnames(counts_filtered_log2)[i]
  sample_idx <- which(rownames(col_data_filtered) == sample_name)
  
  if (length(sample_idx) > 0) {
    condition <- as.character(col_data_filtered$condition[sample_idx])
    cell_line <- as.character(col_data_filtered$cell_line[sample_idx])
    
    if (condition == "Control" && cell_line == "MHCC97H") {
      boxplot_colors_filtered[i] <- "#6BAED6"
    } else if (condition == "Control" && cell_line == "HepG2") {
      boxplot_colors_filtered[i] <- "#2E86AB"
    } else if (condition == "CBD" && cell_line == "MHCC97H") {
      boxplot_colors_filtered[i] <- "#D67BA8"
    } else if (condition == "CBD" && cell_line == "HepG2") {
      boxplot_colors_filtered[i] <- "#A23B72"
    } else {
      boxplot_colors_filtered[i] <- "gray"
    }
  }
}

pdf(file.path(paths$qc_dir, "count_distribution_boxplot_AFTER_filtering.pdf"), width = 14, height = 7)
par(mar = c(10, 4, 4, 6))
boxplot(counts_filtered_log2, 
        main = "Counts distribution AFTER pre-filtering (log2)",
        xlab = "", 
        ylab = "Log2(counts + 1)",
        las = 2, 
        cex.axis = 0.65,
        col = boxplot_colors_filtered,
        border = "black",
        names = colnames(counts_filtered_log2))

legend("topright", 
       legend = c("Control - MHCC97H", "Control - HepG2", "CBD - MHCC97H", "CBD - HepG2"),
       fill = c("#6BAED6", "#2E86AB", "#D67BA8", "#A23B72"),
       cex = 0.8,
       title = "Condition - Cell Line")
grid(nx = NA, ny = NULL, col = "gray90", lty = "dotted")
dev.off()

# Statistics after filtering
cat("  Statistics after pre-filtering:\n")
boxplot_stats_filtered <- data.frame(
  Sample = colnames(counts_filtered_log2),
  Mean = colMeans(counts_filtered_log2),
  Median = apply(counts_filtered_log2, 2, median),
  SD = apply(counts_filtered_log2, 2, sd),
  Min = apply(counts_filtered_log2, 2, min),
  Max = apply(counts_filtered_log2, 2, max)
) %>%
  dplyr::left_join(
    col_data_filtered %>% 
      tibble::rownames_to_column("Sample") %>%
      dplyr::select(Sample, Condition = condition, CellLine = cell_line),
    by = "Sample"
  )

cat("    Overall mean (after filtering):", round(mean(counts_filtered_log2), 4), "\n")
cat("    Overall median (after filtering):", round(median(counts_filtered_log2), 4), "\n")
cat("    Overall SD (after filtering):", round(sd(counts_filtered_log2), 4), "\n")
cat("    Genes included:", nrow(counts_filtered_log2), "\n\n")

readr::write_tsv(boxplot_stats_filtered, 
                 file.path(paths$qc_dir, "boxplot_statistics_AFTER_filtering.tsv"))

# Comparison plot: Before vs After filtering
comparison_df <- data.frame(
  Stage = rep(c("Before Filtering", "After Filtering"), 
              each = length(as.vector(counts_log2))),
  Log2Counts = c(as.vector(counts_log2), as.vector(counts_filtered_log2)),
  Sample = rep(rep(colnames(counts_log2), each = nrow(counts_log2)), 2)
)

# Only include samples that exist in both
common_samples <- intersect(colnames(counts_log2), colnames(counts_filtered_log2))
comparison_df <- comparison_df %>%
  dplyr::filter(Sample %in% common_samples)

comparison_plot <- ggplot2::ggplot(comparison_df, 
                                    ggplot2::aes(x = Stage, y = Log2Counts, fill = Stage)) +
  ggplot2::geom_boxplot(alpha = 0.7, outlier.size = 0.3) +
  ggplot2::scale_fill_manual(values = c("Before Filtering" = "#E8E8E8", 
                                        "After Filtering" = "#4A90E2")) +
  ggplot2::labs(
    title = "Comparison: Before vs After Pre-filtering",
    subtitle = paste0("Before: ", nrow(counts_log2), " genes | After: ", 
                     nrow(counts_filtered_log2), " genes"),
    x = "",
    y = "Log2(counts + 1)"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none")

ggplot2::ggsave(file.path(paths$qc_dir, "comparison_before_after_filtering.pdf"), 
                comparison_plot, width = 8, height = 6)

# Summary comparison table
comparison_summary <- data.frame(
  Metric = c("Number of genes", "Overall mean", "Overall median", "Overall SD", "Min", "Max"),
  Before_Filtering = c(
    nrow(counts_log2),
    round(mean(counts_log2), 4),
    round(median(counts_log2), 4),
    round(sd(counts_log2), 4),
    round(min(counts_log2), 4),
    round(max(counts_log2), 4)
  ),
  After_Filtering = c(
    nrow(counts_filtered_log2),
    round(mean(counts_filtered_log2), 4),
    round(median(counts_filtered_log2), 4),
    round(sd(counts_filtered_log2), 4),
    round(min(counts_filtered_log2), 4),
    round(max(counts_filtered_log2), 4)
  ),
  Change = c(
    paste0(round((nrow(counts_filtered_log2) - nrow(counts_log2)) / nrow(counts_log2) * 100, 1), "%"),
    paste0(round((mean(counts_filtered_log2) - mean(counts_log2)) / mean(counts_log2) * 100, 1), "%"),
    paste0(round((median(counts_filtered_log2) - median(counts_log2)) / median(counts_log2) * 100, 1), "%"),
    paste0(round((sd(counts_filtered_log2) - sd(counts_log2)) / sd(counts_log2) * 100, 1), "%"),
    paste0(round((min(counts_filtered_log2) - min(counts_log2)) / abs(min(counts_log2)) * 100, 1), "%"),
    paste0(round((max(counts_filtered_log2) - max(counts_log2)) / max(counts_log2) * 100, 1), "%")
  )
)

readr::write_tsv(comparison_summary, 
                 file.path(paths$qc_dir, "filtering_comparison_summary.tsv"))
cat("  Comparison summary:\n")
print(comparison_summary, row.names = FALSE)
cat("\n")
cat("  Comparison plots saved to:", paths$qc_dir, "\n")
cat("  Note: PCA and other VST-based plots already use filtered data\n")
cat("        (VST was calculated after pre-filtering)\n\n")


# 7. Differential expression analysis (CBD vs Control) ----------------------------------------
cat("Performing differential expression analysis...\n")
cat("  NOTE: DESeq2 uses the size factors estimated earlier (line 373)\n")
cat("  NOTE: TMM and CPM are only for QC/comparison, NOT used in DE analysis\n\n")

# Run DESeq2
# This uses the size factors already estimated (line 373)
# DESeq() internally applies normalization using those size factors
dds <- DESeq2::DESeq(dds)

# Extract results for CBD vs Control
res <- DESeq2::results(dds, contrast = c("condition", "CBD", "Control"))

# Apply LFC shrinkage
res_shrink <- DESeq2::lfcShrink(dds, coef = "condition_CBD_vs_Control", type = "apeglm")

# Add gene annotation using clusterProfiler::bitr()
# Try to get gene symbols from rownames (could be Ensembl IDs or gene symbols)
gene_ids <- rownames(res_shrink)

cat("  Annotating genes using clusterProfiler...\n")
cat("  Sample gene IDs:", paste(head(gene_ids, 5), collapse = ", "), "...\n")

# Determine ID type more carefully
# Check for Ensembl IDs (ENSG...)
is_ensembl <- any(grepl("^ENSG", gene_ids, ignore.case = TRUE))
# Check if they look like gene symbols (letters, numbers, dashes, but not all numeric)
looks_like_symbol <- any(grepl("^[A-Za-z][A-Za-z0-9-]*$", gene_ids) & !grepl("^ENSG", gene_ids, ignore.case = TRUE))
# Check if they're numeric (might be Entrez IDs)
is_numeric <- all(grepl("^[0-9]+$", gene_ids))

cat("  ID type detection:\n")
cat("    Looks like Ensembl:", is_ensembl, "\n")
cat("    Looks like symbols:", looks_like_symbol, "\n")
cat("    Looks like numeric:", is_numeric, "\n")

# Convert to data frame first
res_tbl <- res_shrink %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id")

# Map genes using clusterProfiler::bitr()
# First, try to determine the actual ID type by testing a sample
test_ids <- head(gene_ids, min(100, length(gene_ids)))
id_type_detected <- NULL

# Try different ID types to see which one works
if (is_ensembl) {
  cat("  Attempting mapping as Ensembl IDs...\n")
  gene_map <- tryCatch({
    result <- clusterProfiler::bitr(
      geneID = gene_ids,
      fromType = "ENSEMBL",
      toType = c("SYMBOL", "ENTREZID", "GENENAME"),
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    )
    id_type_detected <<- "ENSEMBL"
    result %>%
      dplyr::rename(ensgene = ENSEMBL, symbol = SYMBOL, entrez = ENTREZID, description = GENENAME)
  }, error = function(e) {
    cat("    Failed as Ensembl IDs:", e$message, "\n")
    return(NULL)
  })
} else if (is_numeric) {
  cat("  Attempting mapping as Entrez IDs...\n")
  gene_map <- tryCatch({
    result <- clusterProfiler::bitr(
      geneID = as.character(gene_ids),
      fromType = "ENTREZID",
      toType = c("SYMBOL", "ENSEMBL", "GENENAME"),
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    )
    id_type_detected <<- "ENTREZID"
    result %>%
      dplyr::rename(entrez = ENTREZID, symbol = SYMBOL, ensgene = ENSEMBL, description = GENENAME)
  }, error = function(e) {
    cat("    Failed as Entrez IDs:", e$message, "\n")
    return(NULL)
  })
} else {
  # Try as gene symbols
  cat("  Attempting mapping as gene symbols...\n")
  gene_map <- tryCatch({
    # First verify that at least some IDs are valid symbols
    valid_symbols <- tryCatch({
      AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL")
    }, error = function(e) NULL)
    
    if (!is.null(valid_symbols)) {
      # Check if any of our IDs match valid symbols
      matching <- intersect(gene_ids, valid_symbols)
      if (length(matching) > 0) {
        cat("    Found", length(matching), "matching symbols out of", length(gene_ids), "IDs\n")
        result <- clusterProfiler::bitr(
          geneID = gene_ids,
          fromType = "SYMBOL",
          toType = c("ENSEMBL", "ENTREZID", "GENENAME"),
          OrgDb = org.Hs.eg.db::org.Hs.eg.db
        )
        id_type_detected <<- "SYMBOL"
        result %>%
          dplyr::rename(symbol = SYMBOL, ensgene = ENSEMBL, entrez = ENTREZID, description = GENENAME)
      } else {
        cat("    No matching symbols found. IDs may not be standard gene symbols.\n")
        return(NULL)
      }
    } else {
      return(NULL)
    }
  }, error = function(e) {
    cat("    Failed as gene symbols:", e$message, "\n")
    return(NULL)
  })
}

# Process the mapping results
# IMPORTANT: Some genes may map to multiple symbols/IDs (one-to-many relationship)
# We need to handle duplicates to avoid increasing the number of rows
if (!is.null(gene_map) && nrow(gene_map) > 0) {
  cat("  Successfully mapped", nrow(gene_map), "genes using", id_type_detected, "as source\n")
  
  # Check for duplicates in gene_map (one gene mapping to multiple symbols/IDs)
  if (id_type_detected == "ENSEMBL") {
    # Count duplicates
    dup_count <- sum(duplicated(gene_map$ensgene))
    if (dup_count > 0) {
      cat("  Warning: Found", dup_count, "genes with multiple mappings (one-to-many)\n")
      cat("  Taking first mapping for each gene to preserve original gene count\n")
      # Keep only first mapping for each gene
      gene_map <- gene_map %>%
        dplyr::group_by(ensgene) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()
    }
    
    res_tbl <- res_tbl %>%
      dplyr::rename(ensgene = gene_id) %>%
      dplyr::left_join(gene_map, by = "ensgene") %>%
      dplyr::arrange(padj)
  } else if (id_type_detected == "ENTREZID") {
    # Count duplicates
    dup_count <- sum(duplicated(gene_map$entrez))
    if (dup_count > 0) {
      cat("  Warning: Found", dup_count, "genes with multiple mappings (one-to-many)\n")
      cat("  Taking first mapping for each gene to preserve original gene count\n")
      gene_map <- gene_map %>%
        dplyr::group_by(entrez) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()
    }
    
    res_tbl <- res_tbl %>%
      dplyr::mutate(entrez = as.character(gene_id)) %>%
      dplyr::left_join(gene_map, by = "entrez") %>%
      dplyr::arrange(padj)
  } else if (id_type_detected == "SYMBOL") {
    # Count duplicates
    dup_count <- sum(duplicated(gene_map$symbol))
    if (dup_count > 0) {
      cat("  Warning: Found", dup_count, "genes with multiple mappings (one-to-many)\n")
      cat("  Taking first mapping for each gene to preserve original gene count\n")
      gene_map <- gene_map %>%
        dplyr::group_by(symbol) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()
    }
    
    res_tbl <- res_tbl %>%
      dplyr::rename(symbol = gene_id) %>%
      dplyr::left_join(gene_map, by = "symbol") %>%
      dplyr::arrange(padj)
  }
  
  # Verify that we didn't create duplicates
  original_count <- length(gene_ids)
  final_count <- nrow(res_tbl)
  if (final_count > original_count) {
    cat("  Warning: Row count increased from", original_count, "to", final_count, "\n")
    cat("  Removing duplicates to preserve original gene count...\n")
    # Keep only first occurrence of each gene_id
    res_tbl <- res_tbl %>%
      dplyr::distinct(gene_id, .keep_all = TRUE) %>%
      dplyr::arrange(padj)
    cat("  Final count after deduplication:", nrow(res_tbl), "\n")
  }
} else {
  # Fallback: use original IDs as symbols and try alternative mapping methods
  cat("  clusterProfiler::bitr failed. Trying alternative annotation methods...\n")
  
  # Try using biomaRt as fallback
  tryCatch({
    cat("  Attempting annotation with biomaRt...\n")
    mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    # Try to detect ID type from sample
    sample_id <- head(gene_ids, 1)
    if (grepl("^ENSG", sample_id, ignore.case = TRUE)) {
      attrs <- c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id", "description")
      filters <- "ensembl_gene_id"
    } else if (grepl("^[0-9]+$", sample_id)) {
      attrs <- c("entrezgene_id", "hgnc_symbol", "ensembl_gene_id", "description")
      filters <- "entrezgene_id"
    } else {
      attrs <- c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id", "description")
      filters <- "hgnc_symbol"
    }
    
    gene_map_bm <- biomaRt::getBM(
      attributes = attrs,
      filters = filters,
      values = gene_ids,
      mart = mart
    )
    
    if (nrow(gene_map_bm) > 0) {
      cat("  Successfully mapped", nrow(gene_map_bm), "genes using biomaRt\n")
      
      # Handle duplicates in biomaRt results (one gene may map to multiple symbols)
      if (filters == "ensembl_gene_id") {
        # Check for duplicates
        dup_count <- sum(duplicated(gene_map_bm$ensembl_gene_id))
        if (dup_count > 0) {
          cat("  Warning: Found", dup_count, "genes with multiple mappings\n")
          gene_map_bm <- gene_map_bm %>%
            dplyr::group_by(ensembl_gene_id) %>%
            dplyr::slice(1) %>%
            dplyr::ungroup()
        }
        
        gene_map_bm <- gene_map_bm %>%
          dplyr::rename(ensgene = ensembl_gene_id, symbol = hgnc_symbol, 
                       entrez = entrezgene_id, description = description)
        res_tbl <- res_tbl %>%
          dplyr::rename(ensgene = gene_id) %>%
          dplyr::left_join(gene_map_bm, by = "ensgene") %>%
          dplyr::arrange(padj)
      } else if (filters == "entrezgene_id") {
        # Check for duplicates
        dup_count <- sum(duplicated(gene_map_bm$entrezgene_id))
        if (dup_count > 0) {
          cat("  Warning: Found", dup_count, "genes with multiple mappings\n")
          gene_map_bm <- gene_map_bm %>%
            dplyr::group_by(entrezgene_id) %>%
            dplyr::slice(1) %>%
            dplyr::ungroup()
        }
        
        gene_map_bm <- gene_map_bm %>%
          dplyr::rename(entrez = entrezgene_id, symbol = hgnc_symbol, 
                       ensgene = ensembl_gene_id, description = description)
        res_tbl <- res_tbl %>%
          dplyr::mutate(entrez = as.character(gene_id)) %>%
          dplyr::left_join(gene_map_bm, by = "entrez") %>%
          dplyr::arrange(padj)
      } else {
        # Check for duplicates
        dup_count <- sum(duplicated(gene_map_bm$hgnc_symbol))
        if (dup_count > 0) {
          cat("  Warning: Found", dup_count, "genes with multiple mappings\n")
          gene_map_bm <- gene_map_bm %>%
            dplyr::group_by(hgnc_symbol) %>%
            dplyr::slice(1) %>%
            dplyr::ungroup()
        }
        
        gene_map_bm <- gene_map_bm %>%
          dplyr::rename(symbol = hgnc_symbol, ensgene = ensembl_gene_id, 
                       entrez = entrezgene_id, description = description)
        res_tbl <- res_tbl %>%
          dplyr::rename(symbol = gene_id) %>%
          dplyr::left_join(gene_map_bm, by = "symbol") %>%
          dplyr::arrange(padj)
      }
      
      # Final check for duplicates after join
      original_count <- length(gene_ids)
      final_count <- nrow(res_tbl)
      if (final_count > original_count) {
        cat("  Warning: Row count increased from", original_count, "to", final_count, "\n")
        cat("  Removing duplicates...\n")
        res_tbl <- res_tbl %>%
          dplyr::distinct(gene_id, .keep_all = TRUE) %>%
          dplyr::arrange(padj)
        cat("  Final count after deduplication:", nrow(res_tbl), "\n")
      }
    } else {
      cat("  biomaRt returned no results. Using original IDs as symbols.\n")
      res_tbl <- res_tbl %>%
        dplyr::rename(symbol = gene_id) %>%
        dplyr::mutate(ensgene = NA, entrez = NA, description = NA) %>%
        dplyr::arrange(padj)
    }
  }, error = function(e) {
    cat("  biomaRt also failed:", e$message, "\n")
    cat("  Using original IDs as symbols (no annotation available)\n")
    # Final fallback: use original IDs
    res_tbl <- res_tbl %>%
      dplyr::rename(symbol = gene_id) %>%
      dplyr::mutate(ensgene = NA, entrez = NA, description = NA) %>%
      dplyr::arrange(padj)
  })
}

cat("  Annotated", sum(!is.na(res_tbl$symbol)), "out of", nrow(res_tbl), "genes\n")

# Export full results
fs::dir_create(paths$deseq_dir)
readr::write_tsv(res_tbl, file.path(paths$deseq_dir, "DESeq2_results_full.tsv"))

# Significant DEGs (FDR ≤ 0.05, |log2FC| ≥ 1)
# IMPORTANT: In cancer studies, both upregulated and downregulated genes are relevant
# - Upregulated: oncogenes, proliferation genes
# - Downregulated: tumor suppressors, apoptosis genes, cell cycle inhibitors
deg_tbl <- res_tbl %>%
  dplyr::filter(!is.na(padj), padj <= 0.05, abs(log2FoldChange) >= 1) %>%
  dplyr::arrange(padj)

readr::write_tsv(deg_tbl, file.path(paths$deseq_dir, "DESeq2_DEG_FDR0.05_log2FC1.tsv"))

# Separate up and down regulated
# In cancer context:
# - Upregulated: Genes increased by CBD (potential oncogenes suppressed, or protective genes activated)
# - Downregulated: Genes decreased by CBD (potential tumor suppressors activated, or oncogenes suppressed)
deg_up <- deg_tbl %>% dplyr::filter(log2FoldChange > 0)
deg_down <- deg_tbl %>% dplyr::filter(log2FoldChange < 0)

readr::write_tsv(deg_up, file.path(paths$deseq_dir, "DEG_upregulated.tsv"))
readr::write_tsv(deg_down, file.path(paths$deseq_dir, "DEG_downregulated.tsv"))

cat("  Total DEGs (FDR≤0.05, |log2FC|≥1):", nrow(deg_tbl), "\n")
cat("  Upregulated (CBD > Control):", nrow(deg_up), "\n")
cat("  Downregulated (CBD < Control):", nrow(deg_down), "\n")
cat("  NOTE: Both up and down regulated genes are important in cancer studies\n")
cat("        - Downregulated may include tumor suppressors or apoptosis genes\n")
cat("        - Upregulated may include protective genes or oncogene suppression\n\n")

# 7.1. Volcano plot
cat("  Creating volcano plot...\n")

# Prepare data for volcano plot
# Ensure we have valid data and handle missing symbols
volcano_data <- res_tbl %>%
  dplyr::mutate(
    # Use symbol if available, otherwise use gene_id
    label = dplyr::if_else(!is.na(symbol) & symbol != "", symbol, 
                           dplyr::if_else(!is.na(gene_id), gene_id, "Unknown")),
    # Ensure log2FoldChange and padj are numeric and finite
    log2FC = as.numeric(log2FoldChange),
    p_adj = as.numeric(padj),
    # Create significance categories
    significant = !is.na(p_adj) & p_adj <= 0.05 & abs(log2FC) >= 1,
    direction = dplyr::case_when(
      significant & log2FC > 0 ~ "Up",
      significant & log2FC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  dplyr::filter(is.finite(log2FC) & is.finite(p_adj))

cat("    Data points for volcano plot:", nrow(volcano_data), "\n")
cat("    Significant genes:", sum(volcano_data$significant, na.rm = TRUE), "\n")

# Try EnhancedVolcano with error handling
volcano <- tryCatch({
  # Suppress warnings about deprecated size argument
  suppressWarnings({
    EnhancedVolcano::EnhancedVolcano(
      volcano_data,
      lab = volcano_data$label,
      x = "log2FC",
      y = "p_adj",
      pCutoff = 0.05,
      FCcutoff = 1,
      pointSize = 2.0,
      labSize = 3.5,
      title = "CBD vs Control (GSE179661)",
      subtitle = paste0("Total DEGs: ", nrow(deg_tbl), " (Up: ", nrow(deg_up), ", Down: ", nrow(deg_down), ")"),
      caption = paste0("Total genes: ", nrow(volcano_data)),
      xlab = bquote(~Log[2]~ "fold change"),
      ylab = bquote(~-Log[10]~adjusted~italic(P)),
      legendPosition = "right",
      legendLabSize = 12,
      legendIconSize = 5.0,
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      max.overlaps = 20
    )
  })
}, error = function(e) {
  cat("    EnhancedVolcano failed:", e$message, "\n")
  cat("    Creating alternative volcano plot with ggplot2...\n")
  
  # Alternative volcano plot using ggplot2
  volcano_data <- volcano_data %>%
    dplyr::mutate(
      neg_log10_padj = -log10(p_adj + 1e-300)  # Avoid log(0)
    )
  
  ggplot2::ggplot(volcano_data, ggplot2::aes(x = log2FC, y = neg_log10_padj, 
                                              color = direction, label = label)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::scale_color_manual(
      values = c("Up" = "#E31A1C", "Down" = "#1F78B4", "NS" = "gray60"),
      name = "Regulation"
    ) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
    ggplot2::labs(
      title = "Volcano Plot: CBD vs Control (GSE179661)",
      subtitle = paste0("Total DEGs: ", nrow(deg_tbl), " (Up: ", nrow(deg_up), ", Down: ", nrow(deg_down), ")"),
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~Adjusted~P-value),
      caption = paste0("Total genes: ", nrow(volcano_data))
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 11),
      legend.position = "right"
    ) +
    ggrepel::geom_text_repel(
      data = volcano_data %>% dplyr::filter(significant) %>% head(20),
      max.overlaps = 20,
      size = 3,
      box.padding = 0.5
    )
})

# Save volcano plot
tryCatch({
  ggplot2::ggsave(
    file.path(paths$deseq_dir, "volcano_CBD_vs_control.pdf"), 
    volcano, 
    width = 10, 
    height = 8,
    device = "pdf"
  )
  cat("    Volcano plot saved successfully\n")
}, error = function(e) {
  cat("    Error saving volcano plot:", e$message, "\n")
  # Try saving as PNG instead
  tryCatch({
    ggplot2::ggsave(
      file.path(paths$deseq_dir, "volcano_CBD_vs_control.png"), 
      volcano, 
      width = 10, 
      height = 8,
      dpi = 300
    )
    cat("    Volcano plot saved as PNG instead\n")
  }, error = function(e2) {
    warning("Could not save volcano plot: ", e2$message)
  })
})

# 7.2. MA plot
pdf(file.path(paths$deseq_dir, "MA_plot.pdf"), width = 8, height = 6)
DESeq2::plotMA(res_shrink, ylim = c(-5, 5), main = "MA plot: CBD vs Control")
dev.off()

# 7.3. Heatmap of top DEGs
if (nrow(deg_tbl) > 0) {
  top_genes <- deg_tbl %>%
    dplyr::arrange(padj) %>%
    head(50) %>%
    dplyr::pull(symbol)
  
  # Get expression data for top genes
  top_genes_expr <- SummarizedExperiment::assay(vsd)[top_genes[top_genes %in% rownames(vsd)], ]
  
  pdf(file.path(paths$deseq_dir, "heatmap_top50_DEGs.pdf"), width = 10, height = 12)
  pheatmap::pheatmap(
    top_genes_expr,
    scale = "row",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    annotation_col = col_data[, c("condition", "cell_line"), drop = FALSE],
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 7,
    main = "Top 50 DEGs (VST normalized)"
  )
  dev.off()
}


# 8. Pathway and enrichment analysis -----------------------------------------------------------
cat("Performing pathway and enrichment analysis...\n")

# Prepare gene lists for enrichment
if (!is.null(res_tbl$entrez) && sum(!is.na(res_tbl$entrez)) > 0) {
  gene_universe <- unique(res_tbl$entrez[!is.na(res_tbl$entrez)])
  sig_entrez <- unique(deg_tbl$entrez[!is.na(deg_tbl$entrez)])
  sig_entrez_up <- unique(deg_up$entrez[!is.na(deg_up$entrez)])
  sig_entrez_down <- unique(deg_down$entrez[!is.na(deg_down$entrez)])
  
  fs::dir_create(paths$enrich_dir)
  
  # 8.1. GO enrichment (SEA) - Biological Process
  if (length(sig_entrez) > 0) {
    cat("  Running GO enrichment...\n")
    ego_bp <- clusterProfiler::enrichGO(
      gene = sig_entrez,
      universe = gene_universe,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      readable = TRUE
    )
    
    if (nrow(ego_bp) > 0) {
      readr::write_tsv(clusterProfiler::as_tibble(ego_bp), 
                      file.path(paths$enrich_dir, "GO_BP_enrich.tsv"))
      
      pdf(file.path(paths$enrich_dir, "GO_BP_dotplot.pdf"), width = 10, height = 8)
      print(enrichplot::dotplot(ego_bp, showCategory = 20))
      dev.off()
    }
  }
  
  # 8.2. KEGG enrichment
  if (length(sig_entrez) > 0) {
    cat("  Running KEGG enrichment...\n")
    ekegg <- clusterProfiler::enrichKEGG(
      gene = sig_entrez,
      universe = gene_universe,
      organism = "hsa",
      pvalueCutoff = 0.05
    )
    
    if (!is.null(ekegg) && nrow(ekegg) > 0) {
      ekegg <- clusterProfiler::setReadable(ekegg, OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
                                            keyType = "ENTREZID")
      readr::write_tsv(clusterProfiler::as_tibble(ekegg), 
                      file.path(paths$enrich_dir, "KEGG_enrich.tsv"))
      
      pdf(file.path(paths$enrich_dir, "KEGG_dotplot.pdf"), width = 10, height = 8)
      print(enrichplot::dotplot(ekegg, showCategory = 20))
      dev.off()
      
      pdf(file.path(paths$enrich_dir, "KEGG_cnetplot.pdf"), width = 12, height = 10)
      print(enrichplot::cnetplot(ekegg, showCategory = 10, 
                                 foldChange = res_tbl$log2FoldChange))
      dev.off()
    }
  }
  
  # 8.3. GSEA with MSigDB Hallmark
  cat("  Running GSEA (Hallmark pathways)...\n")
  msig_h <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  hallmark_list <- split(msig_h$entrez_gene, msig_h$gs_name)
  
  # Create ranked gene list
  ranked_genes <- res_tbl %>%
    dplyr::filter(!is.na(stat) & !is.na(entrez)) %>%
    dplyr::arrange(desc(stat)) %>%
    dplyr::select(entrez, stat) %>%
    dplyr::distinct(entrez, .keep_all = TRUE) %>%
    tibble::deframe()
  
  if (length(ranked_genes) > 0) {
    fgsea_res <- fgsea::fgsea(
      pathways = hallmark_list,
      stats = ranked_genes,
      minSize = 15,
      maxSize = 500,
      nperm = 10000
    ) %>%
      dplyr::arrange(padj)
    
    readr::write_tsv(fgsea_res, file.path(paths$enrich_dir, "FGSEA_hallmark.tsv"))
    
    # Plot top pathways
    topPathwaysUp <- fgsea_res %>% 
      dplyr::filter(ES > 0, padj < 0.05) %>% 
      head(5) %>% 
      dplyr::pull(pathway)
    topPathwaysDown <- fgsea_res %>% 
      dplyr::filter(ES < 0, padj < 0.05) %>% 
      head(5) %>% 
      dplyr::pull(pathway)
    
    if (length(topPathwaysUp) > 0) {
      pdf(file.path(paths$enrich_dir, "FGSEA_topPathways_up.pdf"), width = 12, height = 8)
      fgsea::plotGseaTable(hallmark_list[topPathwaysUp], ranked_genes, fgsea_res, 
                          gseaParam = 1)
      dev.off()
    }
    
    if (length(topPathwaysDown) > 0) {
      pdf(file.path(paths$enrich_dir, "FGSEA_topPathways_down.pdf"), width = 12, height = 8)
      fgsea::plotGseaTable(hallmark_list[topPathwaysDown], ranked_genes, fgsea_res, 
                          gseaParam = 1)
      dev.off()
    }
  }
} else {
  warning("Could not perform enrichment analysis: missing Entrez IDs")
}

cat("  Enrichment results saved to:", paths$enrich_dir, "\n\n")


# 9. Analysis of key genes and signatures -----------------------------------------------------
cat("Analyzing key genes and pathway signatures...\n")

# 9.1. Key genes from the paper
key_genes <- c("GSDME", "CASP3", "ATF4", "DDIT3", "IGFBP1", "TRPV3")
# DDIT3 is CHOP

# Find these genes in the results
key_genes_data <- res_tbl %>%
  dplyr::filter(symbol %in% key_genes) %>%
  dplyr::select(symbol, log2FoldChange, padj, baseMean)

if (nrow(key_genes_data) > 0) {
  readr::write_tsv(key_genes_data, 
                  file.path(paths$deseq_dir, "key_genes_expression.tsv"))
  
  # Plot key genes expression
  key_genes_expr <- SummarizedExperiment::assay(vsd)[
    key_genes_data$symbol[key_genes_data$symbol %in% rownames(vsd)], 
  ]
  
  if (nrow(key_genes_expr) > 0) {
    pdf(file.path(paths$deseq_dir, "key_genes_heatmap.pdf"), width = 8, height = 6)
    pheatmap::pheatmap(
      key_genes_expr,
      scale = "row",
      annotation_col = col_data[, c("condition", "cell_line"), drop = FALSE],
      main = "Key genes expression (GSDME, CASP3, ATF4, CHOP, IGFBP1, TRPV3)"
    )
    dev.off()
  }
}

# 9.2. Pathway signatures (Pyroptosis, ISR, Glycolysis)
signature_list <- list(
  Pyroptosis = c("GSDME", "GSDMD", "CASP1", "CASP3", "CASP4", "CASP5", 
                 "NLRP3", "PYCARD", "IL1B", "IL18"),
  ISR = c("EIF2AK3", "EIF2AK4", "ATF4", "ATF3", "DDIT3", "PPP1R15A", 
          "PPP1R15B", "XBP1", "ERN1"),
  Glycolysis = c("HK1", "HK2", "PFKM", "PFKP", "ALDOA", "ALDOB", "ENO1", 
                 "ENO2", "PKM", "LDHA", "LDHB", "GAPDH")
)

# Calculate signature scores using ssGSEA
signature_gs <- lapply(signature_list, function(g) {
  intersect(g, rownames(vsd))
})

# Filter out empty gene sets
signature_gs <- signature_gs[sapply(signature_gs, length) > 0]

if (length(signature_gs) > 0) {
  sig_scores <- GSVA::gsva(
    expr = SummarizedExperiment::assay(vsd),
    gset.idx.list = signature_gs,
    method = "ssgsea",
    mx.diff = TRUE,
    verbose = FALSE
  )
  
  signature_df <- sig_scores %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::left_join(
      col_data %>% tibble::rownames_to_column("sample_id"),
      by = "sample_id"
    )
  
  readr::write_tsv(signature_df, 
                  file.path(paths$deseq_dir, "signature_scores.tsv"))
  
  # Plot signature scores
  sig_plot <- signature_df %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(names(signature_gs)), 
      names_to = "signature", 
      values_to = "score"
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = condition, y = score, fill = condition)) +
    ggplot2::geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.1, size = 2) +
    ggplot2::facet_wrap(~signature, scales = "free_y") +
    ggplot2::labs(
      title = "Pathway signature scores (ssGSEA, VST counts)",
      x = "Condition",
      y = "Signature score"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
  
  ggplot2::ggsave(
    file.path(paths$deseq_dir, "signature_scores.pdf"), 
    sig_plot, width = 10, height = 6
  )
  
  # Statistical test for signature differences
  sig_stats <- signature_df %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(names(signature_gs)), 
      names_to = "signature", 
      values_to = "score"
    ) %>%
    dplyr::group_by(signature) %>%
    dplyr::summarise(
      control_mean = mean(score[condition == "Control"]),
      cbd_mean = mean(score[condition == "CBD"]),
      p_value = tryCatch({
        t.test(score ~ condition)$p.value
      }, error = function(e) NA)
    )
  
  readr::write_tsv(sig_stats, 
                  file.path(paths$deseq_dir, "signature_statistics.tsv"))
}

cat("  Key genes and signature analysis completed.\n\n")


# 10. Reproducibility and session info --------------------------------------------------------
cat("Saving session information...\n")
writeLines(capture.output(sessionInfo()), 
          file.path(paths$results_dir, "sessionInfo.txt"))

cat("\n=== Pipeline completed successfully! ===\n")
cat("Results saved to:", paths$results_dir, "\n")
cat("\nSummary:\n")
cat("  - QC plots:", paths$qc_dir, "\n")
cat("  - DE results:", paths$deseq_dir, "\n")
cat("  - Enrichment results:", paths$enrich_dir, "\n")

# End of pipeline --------------------------------------------------------------------------------
