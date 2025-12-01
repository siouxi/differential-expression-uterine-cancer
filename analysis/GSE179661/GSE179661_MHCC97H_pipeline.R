# Pipeline de Análisis GSE179661 - Línea Celular MHCC97H
# Basado en el pipeline original, filtrado para MHCC97H

# 1. Configuración y Carga de Paquetes ------------------------------------------------------------------
required_pkgs <- c(
  "tidyverse", "here", "jsonlite", "glue", "fs",
  "BiocManager", "GEOquery", "DESeq2", "edgeR",
  "SummarizedExperiment", "apeglm", "EnhancedVolcano",
  "pheatmap", "ComplexHeatmap", "RColorBrewer", "plotly",
  "clusterProfiler", "enrichplot", "msigdbr", "fgsea",
  "org.Hs.eg.db", "AnnotationDbi", "biomaRt", "GSVA",
  "corrplot", "vsn",
  "viridis", "VennDiagram", "ggrepel"
)

install_and_load <- function(pkgs) {
  # Verificar qué paquetes ya están instalados (chequeo rápido)
  installed <- rownames(installed.packages())
  new_pkgs <- pkgs[!pkgs %in% installed]
  
  if (length(new_pkgs) > 0) {
    cat("⚠️  Instalando paquetes faltantes (esto puede tardar un poco):", paste(new_pkgs, collapse = ", "), "\n")
    BiocManager::install(new_pkgs, ask = FALSE, update = FALSE)
    cat("✅ Instalación de paquetes completada\n\n")
  } else {
    cat("✅ Todos los paquetes requeridos ya están instalados\n\n")
  }
  
  # Cargar paquetes con mensajes de progreso para los grandes
  cat("Cargando paquetes...\n")
  for (pkg in pkgs) {
    if (pkg %in% c("org.Hs.eg.db", "clusterProfiler", "biomaRt", "msigdbr")) {
      cat("  Cargando", pkg, "(esto puede tardar un momento)...\n")
    }
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE, quietly = TRUE)
    )
  }
  cat("✅ Todos los paquetes cargados\n\n")
}

install_and_load(required_pkgs)


geo_id <- "GSE179661"
cell_line_id <- "MHCC97H"

# Prefijo común para títulos de gráficos
plot_title_prefix <- paste(geo_id, "/", cell_line_id)

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

# Rutas modificadas para resultados específicos de MHCC97H
paths <- list(
  data_dir = file.path(base_dir, "DATA", geo_id),
  results_dir = file.path(base_dir, "results", geo_id, cell_line_id),
  qc_dir = file.path(base_dir, "results", geo_id, cell_line_id, "qc"),
  deseq_dir = file.path(base_dir, "results", geo_id, cell_line_id, "deseq2"),
  enrich_dir = file.path(base_dir, "results", geo_id, cell_line_id, "enrichment")
)

fs::dir_create(unlist(paths))

cat("Rutas del proyecto configuradas para", cell_line_id, ":\n")
cat("  Directorio base:", base_dir, "\n")
cat("  Directorio de datos:", paths$data_dir, "\n")
cat("  Directorio de resultados:", paths$results_dir, "\n\n")


# 3. Cargar y procesar datos de fenotipo ----------------------------------------------------------------
cat("Cargando datos de fenotipo...\n")
pheno_file <- file.path(paths$data_dir, paste0(geo_id, "_pheno_data.csv"))

if (!file.exists(pheno_file)) {
  stop(paste("Archivo de fenotipo no encontrado:", pheno_file))
}

pheno_raw <- readr::read_csv(pheno_file, show_col_types = FALSE)

# Extraer y procesar metadatos de muestras
sample_sheet <- pheno_raw %>%
  dplyr::select(
    sample_id = title,
    geo_accession,
    treatment = `treatment:ch1`,
    cell_line = source_name_ch1
  ) %>%
  dplyr::mutate(
    # Extraer condición de la columna treatment
    condition = dplyr::case_when(
      grepl("none treated", treatment, ignore.case = TRUE) ~ "Control",
      grepl("CBD", treatment, ignore.case = TRUE) ~ "CBD",
      TRUE ~ "Unknown"
    ),
    # Extraer tipo de línea celular
    cell_line = dplyr::case_when(
      grepl("MHCC97H", cell_line, ignore.case = TRUE) ~ "MHCC97H",
      grepl("HepG2", cell_line, ignore.case = TRUE) ~ "HepG2",
      TRUE ~ cell_line
    ),
    # Extraer número de réplica
    replicate = stringr::str_extract(sample_id, "rep\\d+"),
    replicate = as.numeric(stringr::str_replace(replicate, "rep", "")),
    # Crear IDs de muestra que coincidan con los nombres de columna de la matriz de conteo
    sample_id_clean = geo_accession
  ) %>%
  dplyr::mutate(
    condition = factor(condition, levels = c("Control", "CBD")),
    cell_line = factor(cell_line)
  ) %>%
  # FILTRAR SOLO PARA MHCC97H
  dplyr::filter(cell_line == "MHCC97H")

cat("  Número de muestras (MHCC97H):", nrow(sample_sheet), "\n")
cat("  Condiciones:", paste(levels(sample_sheet$condition), collapse = ", "), "\n\n")


# 4. Cargar matriz de conteos génicos crudos ----------------------------------------------------------------
cat("Cargando matriz de conteos génicos crudos...\n")
counts_file <- file.path(paths$data_dir, paste0(geo_id, "_Raw_gene_counts_matrix.txt.gz"))

if (!file.exists(counts_file)) {
  stop(paste("Archivo de conteos no encontrado:", counts_file))
}

# Leer archivo de conteos comprimido
counts_df <- tryCatch({
  readr::read_tsv(counts_file, comment = "#", show_col_types = FALSE)
}, error = function(e) {
  read.table(gzfile(counts_file), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
})

# Identificar columna de ID de gen
gene_col <- NULL
if ("gene_id" %in% colnames(counts_df)) {
  gene_col <- "gene_id"
} else if ("gene_id" %in% colnames(counts_df)) {
  gene_col <- "gene_id"
} else {
  gene_col <- colnames(counts_df)[1]
}

# Preparar matriz de conteo
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

# Asegurar que los conteos sean numéricos
counts_mat <- apply(counts_mat, 2, as.numeric)
rownames(counts_mat) <- rownames(counts_df)

cat("  Dimensiones de matriz de conteo (todas las muestras):", dim(counts_mat), "\n")

# Coincidir nombres de muestras
cat("Coincidiendo nombres de muestras...\n")
sample_cols <- colnames(counts_mat)

# Función auxiliar para normalizar nombres de muestras
normalize_sample_name <- function(name) {
  name_lower <- tolower(name)
  name_lower <- gsub("\\s*rep\\s*", "_", name_lower)
  name_lower <- gsub("\\s+", "_", name_lower)
  name_lower <- gsub("_+", "_", name_lower)
  name_lower <- gsub("^_|_$", "", name_lower)
  return(name_lower)
}

sample_sheet <- sample_sheet %>%
  dplyr::mutate(
    sample_id_normalized = normalize_sample_name(sample_id),
    geo_accession_normalized = tolower(geo_accession)
  )

sample_cols_normalized <- tolower(sample_cols)

# Estrategia 1: Coincidencia directa con sample_id normalizado
matched_indices <- match(sample_cols_normalized, sample_sheet$sample_id_normalized)
matched_samples <- sample_cols[!is.na(matched_indices)]

if (length(matched_samples) > 0) {
  sample_sheet$sample_id_clean[matched_indices[!is.na(matched_indices)]] <- sample_cols[!is.na(matched_indices)]
} else {
  # Estrategia 2: Coincidencia directa con geo_accession
  matched_indices <- match(sample_cols_normalized, sample_sheet$geo_accession_normalized)
  matched_samples <- sample_cols[!is.na(matched_indices)]
  
  if (length(matched_samples) > 0) {
    sample_sheet$sample_id_clean[matched_indices[!is.na(matched_indices)]] <- sample_cols[!is.na(matched_indices)]
  } else {
    # Estrategia 3: Coincidencia parcial
    matched_samples <- character(0)
    matched_pheno_indices <- integer(0)
    
    for (i in seq_along(sample_cols)) {
      col_name <- sample_cols_normalized[i]
      for (j in seq_along(sample_sheet$sample_id_normalized)) {
        pheno_name <- sample_sheet$sample_id_normalized[j]
        if (grepl(gsub("_", ".", col_name), pheno_name, ignore.case = TRUE) ||
            grepl(gsub("_", ".", pheno_name), col_name, ignore.case = TRUE) ||
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
  }
}

# Asegurar que todas las muestras tengan un ID limpio
if (any(is.na(sample_sheet$sample_id_clean))) {
  sample_sheet$sample_id_clean[is.na(sample_sheet$sample_id_clean)] <- sample_sheet$geo_accession[is.na(sample_sheet$sample_id_clean)]
}

# Filtrar y reordenar
matched_samples <- unique(matched_samples[matched_samples %in% sample_cols])
sample_sheet <- sample_sheet %>%
  dplyr::filter(sample_id_clean %in% matched_samples | sample_id_clean %in% sample_cols) %>%
  dplyr::arrange(match(sample_id_clean, sample_cols))

# Asegurar que la matriz de conteos tenga las muestras coincidentes en el orden correcto
available_samples <- intersect(sample_sheet$sample_id_clean, sample_cols)
counts_mat <- counts_mat[, available_samples, drop = FALSE]
sample_sheet <- sample_sheet %>%
  dplyr::filter(sample_id_clean %in% available_samples) %>%
  dplyr::arrange(match(sample_id_clean, colnames(counts_mat)))

cat("\n  ✓ Muestras coincidentes finales (MHCC97H):", length(available_samples), "\n")
cat("  ✓ Columnas de matriz de conteo:", ncol(counts_mat), "\n")
cat("  ✓ Filas de hoja de muestras:", nrow(sample_sheet), "\n\n")


# 5. Crear dataset DESeq2 y pre-filtrado ------------------------------------------------
cat("Creando dataset DESeq2...\n")

col_data <- sample_sheet %>%
  dplyr::select(sample_id = sample_id_clean, condition, replicate) %>%
  as.data.frame()
rownames(col_data) <- col_data$sample_id

# Crear dataset DESeq2 - DISEÑO CAMBIADO A ~ condition
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = col_data,
  design = ~ condition
)

# Pre-filtrado
n_samples <- ncol(dds)
min_samples_with_counts <- max(2, ceiling(n_samples * 0.25))
count_threshold <- 5

keep <- rowSums(DESeq2::counts(dds) >= count_threshold) >= min_samples_with_counts
dds <- dds[keep, ]

cat("  Genes después del filtrado:", nrow(dds), "\n")
cat("  Muestras:", ncol(dds), "\n\n")


# 6. Control de calidad y normalización --------------------------------------------------------
cat("Realizando control de calidad y normalización...\n")

dds <- DESeq2::estimateSizeFactors(dds)

# Colores QC para MHCC97H (Control vs CBD)
qc_colors <- c("#6BAED6", "#D67BA8") # Azul para Control, Rosa para CBD
sample_colors_qc <- ifelse(col_data$condition == "Control", qc_colors[1], qc_colors[2])
names(sample_colors_qc) <- colnames(dds)

dds <- DESeq2::estimateDispersions(dds)
vsd <- DESeq2::vst(dds, blind = FALSE)

# TMM y CPM para QC
dds_edgeR <- edgeR::DGEList(counts = counts_mat[rownames(dds), colnames(dds)])
dds_edgeR <- edgeR::calcNormFactors(dds_edgeR, method = "TMM")
cpm_mat <- edgeR::cpm(dds_edgeR, log = FALSE)

# Guardar resumen de normalización
normalization_summary <- data.frame(
  Sample = colnames(dds),
  Condition = col_data$condition,
  DESeq2_SizeFactor = sizeFactors(dds),
  edgeR_TMM_Factor = dds_edgeR$samples$norm.factors,
  Library_Size = colSums(counts(dds)),
  Mean_CPM = colMeans(cpm_mat)
)
readr::write_tsv(normalization_summary, file.path(paths$qc_dir, "normalization_factors_summary.tsv"))

# Directorio de gráficos QC
fs::dir_create(paths$qc_dir)

# 6.1. Boxplot
counts_log2 <- log2(counts_mat[, rownames(col_data), drop = FALSE] + 1)
pdf(file.path(paths$qc_dir, "count_distribution_boxplot.pdf"), width = 10, height = 7)
par(mar = c(10, 4, 4, 6))
boxplot(counts_log2, 
        main = paste("MHCC97H - Distribución de Conteos"),
        ylab = "Log2(conteos + 1)",
        las = 2, 
        col = sample_colors_qc[colnames(counts_log2)],
        names = colnames(counts_log2))
legend("topright", legend = c("Control", "CBD"), fill = qc_colors, title = "Condición")
dev.off()

# 6.2. PCA plot
pca_data <- DESeq2::plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
pca_plot <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2, color = condition)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::scale_color_manual(values = c("Control" = qc_colors[1], "CBD" = qc_colors[2])) +
  ggplot2::labs(title = "PCA - MHCC97H") +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(paths$qc_dir, "PCA_condition.pdf"), pca_plot, width = 8, height = 6)

# 6.3. Heatmap de distancia de muestras
sample_dists <- dist(t(SummarizedExperiment::assay(vsd)))
sample_dist_mat <- as.matrix(sample_dists)
pdf(file.path(paths$qc_dir, "sample_distance_heatmap.pdf"), width = 8, height = 7)
pheatmap::pheatmap(
  mat = sample_dist_mat,
  color = RColorBrewer::brewer.pal(9, "Blues"),
  clustering_method = "complete",
  main = "Distancias entre Muestras - MHCC97H",
  annotation_col = col_data[, "condition", drop = FALSE]
)
dev.off()


# 7. Análisis de expresión diferencial (CBD vs Control) ----------------------------------------
cat("Realizando análisis de expresión diferencial...\n")

dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds, contrast = c("condition", "CBD", "Control"))
res_shrink <- DESeq2::lfcShrink(dds, coef = "condition_CBD_vs_Control", type = "apeglm")

# Anotación
gene_ids <- rownames(res_shrink)
cat("  Anotando genes usando clusterProfiler...\n")

# Intentar Entrez primero
gene_map <- tryCatch({
  clusterProfiler::bitr(gene_ids, fromType = "ENTREZID", toType = c("SYMBOL", "ENSEMBL", "GENENAME"), OrgDb = org.Hs.eg.db::org.Hs.eg.db, drop = TRUE)
}, error = function(e) NULL)

res_tbl <- res_shrink %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id")

if (!is.null(gene_map)) {
  # Deduplicar mapeo
  gene_map <- gene_map %>% dplyr::distinct(ENTREZID, .keep_all = TRUE)
  
  res_tbl <- res_tbl %>%
    dplyr::mutate(gene_id_char = as.character(gene_id)) %>%
    dplyr::left_join(gene_map, by = c("gene_id_char" = "ENTREZID")) %>%
    dplyr::rename(symbol = SYMBOL, ensgene = ENSEMBL, description = GENENAME, entrez = gene_id_char)
} else {
  # Respaldo si los IDs no son Entrez (ej. Ensembl o Símbolo) - simplificado para este script
  # Asumiendo Entrez basado en el éxito del script original
  res_tbl$symbol <- res_tbl$gene_id
}

# Exportar resultados completos
fs::dir_create(paths$deseq_dir)
readr::write_tsv(res_tbl, file.path(paths$deseq_dir, "DESeq2_results_full.tsv"))

# DEGs Significativos (FDR <= 0.05, |log2FC| >= 2)
deg_tbl <- res_tbl %>%
  dplyr::filter(!is.na(padj), padj <= 0.05, abs(log2FoldChange) >= 2) %>%
  dplyr::arrange(padj)

readr::write_tsv(deg_tbl, file.path(paths$deseq_dir, "DESeq2_DEG_FDR0.05_log2FC2.tsv"))

deg_up <- deg_tbl %>% dplyr::filter(log2FoldChange > 0)
deg_down <- deg_tbl %>% dplyr::filter(log2FoldChange < 0)

readr::write_tsv(deg_up, file.path(paths$deseq_dir, "DEG_upregulated.tsv"))
readr::write_tsv(deg_down, file.path(paths$deseq_dir, "DEG_downregulated.tsv"))

cat("  Total DEGs:", nrow(deg_tbl), "\n")
cat("  Regulados al alza:", nrow(deg_up), "\n")
cat("  Regulados a la baja:", nrow(deg_down), "\n\n")

# 7.1. Volcano plot
volcano_data <- res_tbl %>%
  dplyr::mutate(
    label = dplyr::coalesce(symbol, gene_id),
    significant = !is.na(padj) & padj <= 0.05 & abs(log2FoldChange) >= 2,
    direction = dplyr::case_when(
      significant & log2FoldChange > 0 ~ "Up",
      significant & log2FoldChange < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )

volcano_plot <- ggplot2::ggplot(volcano_data, ggplot2::aes(x = log2FoldChange, y = -log10(padj), color = direction)) +
  ggplot2::geom_point(alpha = 0.6, size = 1.5) +
  ggplot2::scale_color_manual(values = c("Up" = "#E31A1C", "Down" = "#1F78B4", "NS" = "gray60")) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggplot2::geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  ggplot2::labs(title = "Volcano Plot - MHCC97H (CBD vs Control)") +
  ggplot2::theme_bw()

ggplot2::ggsave(file.path(paths$deseq_dir, "volcano_CBD_vs_control.pdf"), volcano_plot, width = 8, height = 6)

# 7.3. Heatmap de top DEGs
if (nrow(deg_tbl) > 0) {
  top_genes <- deg_tbl %>% head(50) %>% dplyr::pull(symbol)
  # Manejar símbolos faltantes
  top_genes <- top_genes[!is.na(top_genes)]
  
  # Mapear símbolos de vuelta a rownames (usualmente IDs Entrez)
  # Esta parte depende de cómo estén configurados los rownames en vsd. Son gene_ids.
  # Necesitamos encontrar los gene_ids correspondientes a estos símbolos
  top_gene_ids <- res_tbl %>% dplyr::filter(symbol %in% top_genes) %>% dplyr::pull(gene_id)
  
  top_genes_expr <- SummarizedExperiment::assay(vsd)[top_gene_ids[top_gene_ids %in% rownames(vsd)], ]
  
  # Usar símbolos como nombres de fila para heatmap si están disponibles
  rownames(top_genes_expr) <- res_tbl$symbol[match(rownames(top_genes_expr), res_tbl$gene_id)]
  
  pdf(file.path(paths$deseq_dir, "heatmap_top50_DEGs.pdf"), width = 8, height = 10)
  pheatmap::pheatmap(
    top_genes_expr,
    scale = "row",
    annotation_col = col_data[, "condition", drop = FALSE],
    main = "Top 50 DEGs - MHCC97H"
  )
  dev.off()
}

# 8. Análisis de Enriquecimiento (Simplificado) -----------------------------------------------------------
cat("Realizando análisis de enriquecimiento...\n")
fs::dir_create(paths$enrich_dir)

# Preparar listas de genes
gene_universe <- as.character(res_tbl$entrez[!is.na(res_tbl$entrez)])
sig_entrez <- as.character(deg_tbl$entrez[!is.na(deg_tbl$entrez)])

if (length(sig_entrez) > 0) {
  # GO BP
  ego_bp <- clusterProfiler::enrichGO(
    gene = sig_entrez,
    universe = gene_universe,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    readable = TRUE
  )
  
  if (!is.null(ego_bp)) {
    readr::write_tsv(as.data.frame(ego_bp), file.path(paths$enrich_dir, "GO_BP_all_DEGs.tsv"))
    pdf(file.path(paths$enrich_dir, "GO_BP_dotplot.pdf"), width = 10, height = 8)
    print(enrichplot::dotplot(ego_bp, showCategory = 20, title = "GO BP - MHCC97H"))
    dev.off()
  }
  
  # KEGG
  ekegg <- clusterProfiler::enrichKEGG(
    gene = sig_entrez,
    universe = gene_universe,
    organism = "hsa",
    pvalueCutoff = 0.05
  )
  
  if (!is.null(ekegg)) {
    readr::write_tsv(as.data.frame(ekegg), file.path(paths$enrich_dir, "KEGG_all_DEGs.tsv"))
  }
}

cat("¡Análisis para MHCC97H completado!\n")
