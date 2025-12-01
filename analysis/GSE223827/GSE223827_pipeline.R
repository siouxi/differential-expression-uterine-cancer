## GSE223827 - Pipeline BRB-seq (unq refseq UMI)
## Carga de datos de expresión y metadatos

## Paquetes básicos
suppressPackageStartupMessages({
  library(tidyverse)
  library(matrixStats)
  library(pheatmap)
  library(RColorBrewer)
  library(uwot)
  library(DESeq2)
  library(limma)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

geo_id <- "GSE223827"

## Detectar base directory igual que en GSE285498 --------------------
## Asumimos que este script vive en analysis/GSE223827 y que
## la raíz del proyecto es dos niveles arriba.

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

if (is.null(script_dir)) {
  # Fallback: usar current working directory
  script_dir <- getwd()
}

base_dir <- dirname(dirname(script_dir))

paths <- list(
  data_dir    = file.path(base_dir, "DATA", geo_id),
  results_dir = file.path(base_dir, "results", geo_id)
)

if (!dir.exists(paths$results_dir)) {
  dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
}

cat("Rutas del proyecto configuradas para", geo_id, ":\n")
cat("  Directorio base:   ", base_dir, "\n")
cat("  Directorio de datos:   ", paths$data_dir, "\n")
cat("  Directorio de resultados:", paths$results_dir, "\n\n")

## Rutas a los archivos dentro de DATA/GSE223827 ---------------------

brb_file   <- file.path(paths$data_dir, paste0(geo_id, "_BRBSeq.unq.refseq.umi.dat.txt.gz"))
pheno_file <- file.path(paths$data_dir, paste0(geo_id, "_pheno_data.csv"))

if (!file.exists(brb_file)) {
  stop("No se encontró el archivo BRB-seq: ", brb_file)
}
if (!file.exists(pheno_file)) {
  stop("No se encontró el archivo de pheno_data: ", pheno_file)
}

## ------------------------------------------------------------------
## 1. Cargar matriz de conteos BRB-seq (UMIs por gen y muestra)
## ------------------------------------------------------------------
## El archivo .dat.txt.gz es una matriz de expresión en formato texto,
## comprimida con gzip. Usamos gzfile() + read.table().

stopifnot(file.exists(brb_file))

expr_mat <- read.table(
  gzfile(brb_file),
  header       = TRUE,   # primera fila = nombres de columnas (muestras)
  sep          = "\t",   # BRB-seq dat suele ser tab-delimited
  row.names    = 1,      # primera columna = ID de gen (RefSeq)
  check.names  = FALSE,  # no forzar nombres válidos de columnas
  comment.char = "",
  quote        = ""
)

## expr_mat: filas = genes (RefSeq), columnas = muestras
dim(expr_mat)

## Confirmar que no existan nombres de muestra duplicados en la matriz
dup_expr_mat <- colnames(expr_mat)[duplicated(colnames(expr_mat))]
if (length(dup_expr_mat)) {
  stop("Se detectaron columnas duplicadas en expr_mat: ", paste(dup_expr_mat, collapse = ", "))
} else {
  cat("Verificación OK: no hay columnas duplicadas en expr_mat (", ncol(expr_mat), " muestras ).\n\n")
}

## Resumen rápido de tratamientos inferidos desde los nombres de columnas
col_tokens <- strsplit(colnames(expr_mat), "_", fixed = TRUE)
treatment_labels <- vapply(
  col_tokens,
  function(tokens) {
    if (length(tokens) >= 2) {
      # último token = ID del explanto (EHxxxx); nos quedamos con el penúltimo
      tokens[length(tokens) - 1]
    } else {
      NA_character_
    }
  },
  FUN.VALUE = character(1)
)
cat("Tratamientos detectados a partir de colnames(expr_mat) (penúltimo bloque):\n")
print(sort(table(treatment_labels, useNA = "ifany"), decreasing = TRUE))
cat("Número total de tratamientos únicos:", dplyr::n_distinct(treatment_labels, na.rm = TRUE), "\n\n")

## ------------------------------------------------------------------
## 2. Cargar metadatos de muestras (pheno data desde GEO)
## ------------------------------------------------------------------

stopifnot(file.exists(pheno_file))

 pheno <- read.csv(pheno_file, stringsAsFactors = FALSE, check.names = FALSE)

## Nos quedamos con columnas clave y construimos variables útiles
pheno_clean <- pheno %>%
  transmute(
    sample_title   = title,
    geo_accession  = geo_accession,
    age_pcw        = `age (pcw):ch1`,
    dev_stage      = `developmental stage:ch1`,
    exposure       = `exposure:ch1`,
    gender         = `gender:ch1`,
    outlier        = `outlier:ch1`
  )

## Confirmar que no existan muestras duplicadas en pheno_clean
dup_pheno <- pheno_clean$sample_title[duplicated(pheno_clean$sample_title)]
if (length(dup_pheno)) {
  stop("Se detectaron sample_title duplicados en pheno_clean: ", paste(dup_pheno, collapse = ", "))
} else {
  cat("Verificación OK: no hay sample_title duplicados en pheno_clean (", nrow(pheno_clean), " registros ).\n\n")
}

## ------------------------------------------------------------------
## 3. Construir sample_sheet alineado con expr_mat
## ------------------------------------------------------------------

common_samples <- intersect(colnames(expr_mat), pheno_clean$sample_title)
if (length(common_samples) == 0) {
  stop("No se encontraron muestras comunes entre expr_mat y pheno_data.")
}

sample_sheet <- pheno_clean %>%
  filter(sample_title %in% common_samples) %>%
  distinct(sample_title, .keep_all = TRUE)

sample_sheet <- sample_sheet[match(colnames(expr_mat), sample_sheet$sample_title), ]
if (any(is.na(sample_sheet$sample_title))) {
  stop("Ocurrió un desajuste al ordenar sample_sheet. Revisa los nombres de columnas.")
}
stopifnot(nrow(sample_sheet) == ncol(expr_mat))

sample_sheet <- as.data.frame(sample_sheet)
rownames(sample_sheet) <- NULL

sample_sheet <- sample_sheet %>%
  mutate(
    testis_id = stringr::str_extract(sample_title, "EH\\d+$"),
    exposure  = factor(exposure, levels = c("DMSO", "EtOH", "DMSO-EtOH", "CBD", "THC", "CBD-THC")),
    dev_stage = factor(dev_stage, levels = sort(unique(dev_stage))),
    age_pcw   = factor(age_pcw, levels = unique(age_pcw)),
    outlier   = factor(outlier, levels = unique(outlier))
  ) %>%
  tibble::column_to_rownames("sample_title")

expr_mat <- expr_mat[, rownames(sample_sheet)]

cat("Resumen de condiciones (exposure):\n")
print(table(sample_sheet$exposure, useNA = "ifany"))
cat("\nResumen de testículos (nº de explantos por EH):\n")
print(sort(table(sample_sheet$testis_id), decreasing = TRUE))
cat("\n")

qc_dir <- file.path(paths$results_dir, "qc")
if (!dir.exists(qc_dir)) {
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
}

deg_dir <- file.path(paths$results_dir, "deg")
if (!dir.exists(deg_dir)) {
  dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)
}

save_plot <- function(plot_obj, filename, width = 10, height = 6) {
  outfile <- file.path(qc_dir, filename)
  ggplot2::ggsave(outfile, plot_obj, width = width, height = height)
  cat("Plot guardado en:", outfile, "\n")
}

## ------------------------------------------------------------------
## 4. Métricas básicas de QC por muestra
## ------------------------------------------------------------------

library_size <- colSums(expr_mat)
detected_genes <- colSums(expr_mat > 0)
pct_zeros <- (1 - detected_genes / nrow(expr_mat)) * 100

qc_metrics <- tibble::tibble(
  sample        = colnames(expr_mat),
  library_size  = library_size,
  library_size_million = library_size / 1e6,
  detected_genes = detected_genes,
  pct_zeros     = pct_zeros,
  exposure      = sample_sheet$exposure,
  dev_stage     = sample_sheet$dev_stage,
  age_pcw       = sample_sheet$age_pcw,
  testis_id     = sample_sheet$testis_id
)

cat("Resumen numérico global de métricas QC:\n")
print(summary(qc_metrics[, c("library_size", "detected_genes", "pct_zeros")]))
cat("\nResumen por exposición (medianas):\n")
qc_by_exposure <- qc_metrics %>%
  group_by(exposure) %>%
  summarise(
    n_muestras       = dplyr::n(),
    med_library      = stats::median(library_size),
    med_detected     = stats::median(detected_genes),
    med_pct_zeros    = stats::median(pct_zeros),
    .groups = "drop"
  )
print(qc_by_exposure)
cat("\n")

qc_rules <- list(
  min_library   = 8e5,   # 0.8 millones de UMIs
  min_detected  = 3000,  # genes con UMIs > 0
  max_pct_zeros = 85     # porcentaje máximo de ceros aceptable
)

cat("Criterios sugeridos para clasificar muestras (puedes ajustarlos):\n")
print(qc_rules)

qc_metrics <- qc_metrics %>%
  mutate(
    qc_flag = dplyr::case_when(
      library_size < qc_rules$min_library ~ "LOW_LIB",
      detected_genes < qc_rules$min_detected ~ "LOW_GENES",
      pct_zeros > qc_rules$max_pct_zeros ~ "HIGH_ZEROS",
      TRUE ~ "OK"
    )
  )

cat("\nConteo de muestras por categoría QC:\n")
print(table(qc_metrics$qc_flag))

if (any(qc_metrics$qc_flag != "OK")) {
  cat("\nDetalle de muestras a revisar:\n")
  print(
    qc_metrics %>%
      filter(qc_flag != "OK") %>%
      arrange(qc_flag, library_size)
  )
  cat("\n")
}

readr::write_tsv(qc_metrics, file.path(qc_dir, "qc_metrics.tsv"))
cat("Metrics por muestra guardadas en qc_metrics.tsv\n\n")

sample_sheet$library_size_million <- qc_metrics$library_size_million

detected_limits <- range(qc_metrics$detected_genes, na.rm = TRUE)
padding <- 0.02 * diff(detected_limits)
if (!is.finite(padding)) padding <- 0
detected_limits <- detected_limits + c(-padding, padding)
detected_limits[1] <- max(detected_limits[1], 0)

library_plot <- qc_metrics %>%
  mutate(sample = factor(sample, levels = sample[order(library_size)])) %>%
  ggplot2::ggplot(ggplot2::aes(x = sample, y = library_size / 1e6, fill = exposure)) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::labs(x = "Muestra", y = "UMIs (millones)", fill = "Exposición",
                title = "Tamaño de biblioteca por muestra") +
  ggplot2::theme_minimal(base_size = 11)
save_plot(library_plot, "library_size_barplot.pdf", width = 10, height = 14)

zeros_plot <- ggplot2::ggplot(qc_metrics,
                              ggplot2::aes(x = library_size / 1e6, y = pct_zeros,
                                           color = exposure, shape = dev_stage)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::labs(x = "UMIs (millones)", y = "% de ceros",
                title = "Porcentaje de ceros vs tamaño de biblioteca") +
  ggplot2::theme_minimal(base_size = 11)
save_plot(zeros_plot, "pct_zeros_vs_library.pdf")

detected_plot <- ggplot2::ggplot(qc_metrics,
                                 ggplot2::aes(x = exposure, y = detected_genes,
                                              color = dev_stage)) +
  ggplot2::geom_jitter(width = 0.2, height = 0, size = 2) +
  ggplot2::stat_summary(fun = median, geom = "crossbar",
                        width = 0.45, color = "black") +
  ggplot2::labs(x = "Exposición", y = "Genes detectados (>0 UMI)",
                title = "Genes detectados por tratamiento") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::coord_cartesian(ylim = detected_limits)
save_plot(detected_plot, "detected_genes_by_exposure.pdf")

stage_levels_detected <- levels(sample_sheet$dev_stage)
if (length(stage_levels_detected) == 0) {
  stage_levels_detected <- unique(sample_sheet$dev_stage)
}
stage_levels_detected <- stage_levels_detected[!is.na(stage_levels_detected)]

for (stage in stage_levels_detected) {
  stage_df <- qc_metrics %>% filter(dev_stage == stage)
  if (nrow(stage_df) < 3) {
    cat("Saltando plot de genes detectados para etapa", stage, "- menos de 3 muestras.\n")
    next
  }
  stage_plot <- ggplot2::ggplot(stage_df,
                                ggplot2::aes(x = exposure, y = detected_genes, color = exposure)) +
    ggplot2::geom_jitter(width = 0.2, height = 0, size = 2) +
    ggplot2::stat_summary(fun = median, geom = "crossbar",
                          width = 0.45, color = "black") +
    ggplot2::labs(
      title = paste0("Genes detectados por tratamiento (", stage, ")"),
      x = "Exposición",
      y = "Genes detectados (>0 UMI)"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::coord_cartesian(ylim = detected_limits)

  stage_slug <- gsub("[^0-9A-Za-z]+", "_", stage)
  save_plot(stage_plot, paste0("detected_genes_by_exposure_", stage_slug, ".pdf"))
}

## ------------------------------------------------------------------
## 5. Filtrado de genes y normalización (logCPM)
## ------------------------------------------------------------------

min_umi     <- 5
min_samples <- 3
keep_genes <- rowSums(expr_mat >= min_umi) >= min_samples
expr_filtered <- expr_mat[keep_genes, ]
genes_retained <- sum(keep_genes)
cat("Genes retenidos tras el filtro (>= ", min_umi, " UMIs en ", min_samples,
    " muestras): ", genes_retained, " de ", length(keep_genes), "\n\n", sep = "")

if (genes_retained == 0) {
  stop("Ningún gen superó el filtro definido. Ajusta los umbrales.")
}

library_size_filt <- colSums(expr_filtered)
if (any(library_size_filt == 0)) {
  stop("Al menos una muestra tiene 0 UMIs después del filtrado; revisa los parámetros.")
}
cpm <- t(t(expr_filtered) / library_size_filt * 1e6)
log_cpm <- log2(cpm + 1)

## ------------------------------------------------------------------
## 6. PCA y correlación entre muestras
## ------------------------------------------------------------------

top_var <- min(2000, nrow(log_cpm))
top_genes <- order(matrixStats::rowVars(log_cpm), decreasing = TRUE)[seq_len(top_var)]
log_cpm_top <- log_cpm[top_genes, ]

pca_res <- prcomp(t(log_cpm_top), center = TRUE, scale. = FALSE)
pca_scores <- tibble::as_tibble(pca_res$x[, 1:2], rownames = "sample")
names(pca_scores)[names(pca_scores) == "sample"] <- "sample"
pc_cols <- setdiff(names(pca_scores), "sample")
names(pca_scores)[names(pca_scores) == pc_cols[1]] <- "PC1"
names(pca_scores)[names(pca_scores) == pc_cols[2]] <- "PC2"
pca_df <- qc_metrics %>%
  left_join(pca_scores, by = "sample")

pvar <- (pca_res$sdev^2) / sum(pca_res$sdev^2)

pca_plot <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2,
                                                 color = exposure, shape = dev_stage)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::labs(
    title = "PCA (logCPM, top 2000 genes)",
    x = paste0("PC1 (", scales::percent(pvar[1]), ")"),
    y = paste0("PC2 (", scales::percent(pvar[2]), ")")
  ) +
  ggplot2::theme_minimal(base_size = 11)
save_plot(pca_plot, "pca_logcpm_top2000.pdf")

pca_plot_testis <- ggplot2::ggplot(pca_df,
                                   ggplot2::aes(x = PC1, y = PC2,
                                                color = testis_id, shape = exposure)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::labs(
    title = "PCA coloreado por testículo (EH)",
    x = paste0("PC1 (", scales::percent(pvar[1]), ")"),
    y = paste0("PC2 (", scales::percent(pvar[2]), ")"),
    color = "Testis ID",
    shape = "Exposure"
  ) +
  ggplot2::theme_minimal(base_size = 11)
save_plot(pca_plot_testis, "pca_logcpm_top2000_by_testis.pdf")

pca_plot_qc <- ggplot2::ggplot(pca_df,
                               ggplot2::aes(x = PC1, y = PC2,
                                            color = library_size_million,
                                            shape = qc_flag)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::scale_color_distiller(palette = "YlGnBu", direction = 1,
                                 name = "UMIs (M)") +
  ggplot2::labs(
    title = "PCA coloreado por tamaño de biblioteca",
    x = paste0("PC1 (", scales::percent(pvar[1]), ")"),
    y = paste0("PC2 (", scales::percent(pvar[2]), ")"),
    shape = "QC flag"
  ) +
  ggplot2::theme_minimal(base_size = 11)
save_plot(pca_plot_qc, "pca_logcpm_top2000_by_QC.pdf")

## PCA corrigiendo por testis_id (removeBatchEffect) -----------------
cat("El 22%/9% de varianza en el PCA se explica sobre todo por testis_id.\n")
cat("Generando PCA de logCPM corregido por batch (testículo) para aislar exposición...\n\n")

design_exposure <- model.matrix(~ exposure, data = sample_sheet)
log_cpm_batch_corrected <- limma::removeBatchEffect(
  log_cpm_top,
  batch = sample_sheet$testis_id,
  design = design_exposure
)

pca_batch <- prcomp(t(log_cpm_batch_corrected), center = TRUE, scale. = FALSE)
pca_batch_scores <- tibble::as_tibble(pca_batch$x[, 1:2], rownames = "sample")
names(pca_batch_scores)[names(pca_batch_scores) == "sample"] <- "sample"
pc_cols_batch <- setdiff(names(pca_batch_scores), "sample")
names(pca_batch_scores)[names(pca_batch_scores) == pc_cols_batch[1]] <- "PC1"
names(pca_batch_scores)[names(pca_batch_scores) == pc_cols_batch[2]] <- "PC2"

pca_batch_df <- qc_metrics %>%
  left_join(pca_batch_scores, by = "sample")

pvar_batch <- (pca_batch$sdev^2) / sum(pca_batch$sdev^2)

pca_batch_plot <- ggplot2::ggplot(
  pca_batch_df,
  ggplot2::aes(x = PC1, y = PC2, color = exposure, shape = dev_stage)
) +
  ggplot2::geom_point(size = 3) +
  ggplot2::labs(
    title = "PCA (logCPM corregido por testículo)",
    subtitle = "removeBatchEffect(batch = testis_id, design = ~ exposure)",
    x = paste0("PC1 (", scales::percent(pvar_batch[1]), ")"),
    y = paste0("PC2 (", scales::percent(pvar_batch[2]), ")")
  ) +
  ggplot2::theme_minimal(base_size = 11)
save_plot(pca_batch_plot, "pca_logcpm_top2000_batch_corrected.pdf")

## PCA corrigiendo por tamaño de biblioteca (covariate) --------------
cat("Probando corrección continua por tamaño de biblioteca (UMIs)...\n")

covariate_umIs <- matrix(sample_sheet$library_size_million, ncol = 1)
colnames(covariate_umIs) <- "library_size_million"
design_cov <- model.matrix(~ exposure, data = sample_sheet)

log_cpm_cov_corrected <- limma::removeBatchEffect(
  log_cpm_top,
  covariates = covariate_umIs,
  design = design_cov
)

pca_cov <- prcomp(t(log_cpm_cov_corrected), center = TRUE, scale. = FALSE)
pca_cov_scores <- tibble::as_tibble(pca_cov$x[, 1:2], rownames = "sample")
pc_cols_cov <- setdiff(names(pca_cov_scores), "sample")
names(pca_cov_scores)[names(pca_cov_scores) == pc_cols_cov[1]] <- "PC1"
names(pca_cov_scores)[names(pca_cov_scores) == pc_cols_cov[2]] <- "PC2"
pca_cov_df <- qc_metrics %>%
  left_join(pca_cov_scores, by = "sample")

pvar_cov <- (pca_cov$sdev^2) / sum(pca_cov$sdev^2)

pca_cov_plot <- ggplot2::ggplot(
  pca_cov_df,
  ggplot2::aes(x = PC1, y = PC2, color = exposure, shape = dev_stage)
) +
  ggplot2::geom_point(size = 3) +
  ggplot2::labs(
    title = "PCA (logCPM corregido por UMIs)",
    subtitle = "removeBatchEffect(covariates = library_size_million)",
    x = paste0("PC1 (", scales::percent(pvar_cov[1]), ")"),
    y = paste0("PC2 (", scales::percent(pvar_cov[2]), ")")
  ) +
  ggplot2::theme_minimal(base_size = 11)
save_plot(pca_cov_plot, "pca_logcpm_top2000_library_corrected.pdf")

## UMAP global -------------------------------------------------------
set.seed(42)
umap_coords <- uwot::umap(
  t(log_cpm_top),
  n_neighbors = 15,
  min_dist = 0.3,
  metric = "cosine",
  verbose = FALSE
)
umap_df <- as.data.frame(umap_coords)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$sample <- rownames(t(log_cpm_top))
umap_df <- umap_df %>%
  left_join(qc_metrics, by = "sample")

umap_plot <- ggplot2::ggplot(umap_df,
                             ggplot2::aes(x = UMAP1, y = UMAP2,
                                          color = exposure, shape = dev_stage)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::labs(
    title = "UMAP (logCPM, top 2000 genes)",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  ggplot2::theme_minimal(base_size = 11)
save_plot(umap_plot, "umap_logcpm_top2000.pdf")

stage_levels_umap <- stage_levels
if (length(stage_levels_umap) == 0) {
  stage_levels_umap <- unique(sample_sheet$dev_stage)
}
stage_levels_umap <- stage_levels_umap[!is.na(stage_levels_umap)]

for (stage in stage_levels_umap) {
  stage_samples <- rownames(sample_sheet)[sample_sheet$dev_stage == stage]
  if (length(stage_samples) < 5) {
    cat("Saltando UMAP por etapa", stage, "- menos de 5 muestras.\n")
    next
  }
  stage_matrix <- log_cpm_top[, stage_samples, drop = FALSE]
  stage_umap <- uwot::umap(
    t(stage_matrix),
    n_neighbors = min(15, length(stage_samples) - 1),
    min_dist = 0.3,
    metric = "cosine",
    verbose = FALSE
  )
  stage_df <- as.data.frame(stage_umap)
  colnames(stage_df) <- c("UMAP1", "UMAP2")
  stage_df$sample <- rownames(stage_df)
  stage_df <- stage_df %>%
    left_join(qc_metrics, by = "sample")

  stage_plot <- ggplot2::ggplot(stage_df,
                                ggplot2::aes(x = UMAP1, y = UMAP2, color = exposure)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      title = paste0("UMAP (", stage, ", logCPM top 2000)"),
      x = "UMAP1",
      y = "UMAP2"
    ) +
    ggplot2::theme_minimal(base_size = 11)

  stage_slug <- gsub("[^0-9A-Za-z]+", "_", stage)
  save_plot(stage_plot, paste0("umap_logcpm_", stage_slug, ".pdf"))
}

## ------------------------------------------------------------------
## 7. Transformación VST (DESeq2) y PCA
## ------------------------------------------------------------------

cat("Aplicando transformación VST (DESeq2) para estabilizar la varianza...\n")
cat("VST ajusta por tamaño de biblioteca y comprime diferencias, similar a log2 pero con\n")
cat("parámetros derivados del modelo de conteos. Es útil cuando las muestras tienen\n")
cat("profundidades muy distintas porque reduce la dependencia entre media y varianza.\n\n")

deseq_counts <- round(expr_filtered)
dds <- DESeqDataSetFromMatrix(
  countData = deseq_counts,
  colData   = sample_sheet,
  design    = ~ exposure
)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

vsd_top_n <- min(2000, nrow(vsd_mat))
vsd_top_genes <- order(matrixStats::rowVars(vsd_mat), decreasing = TRUE)[seq_len(vsd_top_n)]
vsd_pca <- prcomp(t(vsd_mat[vsd_top_genes, ]), center = TRUE, scale. = FALSE)

vsd_scores <- tibble::as_tibble(vsd_pca$x[, 1:2], rownames = "sample")
vsd_df <- qc_metrics %>%
  left_join(vsd_scores, by = "sample")

vsd_pvar <- (vsd_pca$sdev^2) / sum(vsd_pca$sdev^2)

vsd_pca_plot <- ggplot2::ggplot(vsd_df,
                                ggplot2::aes(x = PC1, y = PC2,
                                             color = exposure, shape = dev_stage)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::labs(
    title = "PCA (DESeq2 VST, top 2000 genes)",
    x = paste0("PC1 (", scales::percent(vsd_pvar[1]), ")"),
    y = paste0("PC2 (", scales::percent(vsd_pvar[2]), ")")
  ) +
  ggplot2::theme_minimal(base_size = 11)
save_plot(vsd_pca_plot, "pca_vst_top2000.pdf")

stage_levels_vst <- stage_levels
if (length(stage_levels_vst) == 0) {
  stage_levels_vst <- unique(sample_sheet$dev_stage)
}
stage_levels_vst <- stage_levels_vst[!is.na(stage_levels_vst)]

for (stage in stage_levels_vst) {
  stage_samples <- rownames(sample_sheet)[sample_sheet$dev_stage == stage]
  if (length(stage_samples) < 3) {
    cat("Saltando PCA VST por etapa", stage, "- menos de 3 muestras.\n")
    next
  }

  stage_vsd <- vsd_mat[vsd_top_genes, stage_samples, drop = FALSE]
  stage_pca <- prcomp(t(stage_vsd), center = TRUE, scale. = FALSE)
  stage_scores <- tibble::as_tibble(stage_pca$x[, 1:2], rownames = "sample")
  stage_df <- qc_metrics %>%
    filter(sample %in% stage_samples) %>%
    left_join(stage_scores, by = "sample")

  stage_pvar <- (stage_pca$sdev^2) / sum(stage_pca$sdev^2)

  stage_plot <- ggplot2::ggplot(stage_df,
                                ggplot2::aes(x = PC1, y = PC2, color = exposure)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      title = paste0("PCA VST (", stage, ", top 2000 genes)"),
      x = paste0("PC1 (", scales::percent(stage_pvar[1]), ")"),
      y = paste0("PC2 (", scales::percent(stage_pvar[2]), ")")
    ) +
    ggplot2::theme_minimal(base_size = 11)

  stage_slug <- gsub("[^0-9A-Za-z]+", "_", stage)
  save_plot(stage_plot, paste0("pca_vst_top2000_", stage_slug, ".pdf"))
}

sample_cor <- stats::cor(log_cpm_top, method = "spearman")
annotation_df <- sample_sheet[, c("exposure", "dev_stage", "age_pcw", "testis_id")]

heatmap_file <- file.path(qc_dir, "sample_correlation_heatmap.pdf")
pdf(heatmap_file, width = 11, height = 10)
pheatmap::pheatmap(
  sample_cor,
  color = colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100),
  annotation_col = annotation_df,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  show_colnames = FALSE,
  show_rownames = FALSE,
  main = "Correlación (logCPM, top 2000 genes)"
)
dev.off()
cat("Heatmap guardado en:", heatmap_file, "\n\n")

## PCA estratificada por etapa de desarrollo -------------------------
stage_levels <- levels(sample_sheet$dev_stage)
if (length(stage_levels) == 0) {
  stage_levels <- unique(sample_sheet$dev_stage)
}
stage_levels <- stage_levels[!is.na(stage_levels)]

for (stage in stage_levels) {
  stage_samples <- rownames(sample_sheet)[sample_sheet$dev_stage == stage]
  if (length(stage_samples) < 3) {
    cat("Saltando PCA por etapa", stage, "- menos de 3 muestras.\n")
    next
  }

  stage_matrix <- log_cpm_top[, stage_samples, drop = FALSE]
  stage_pca <- prcomp(t(stage_matrix), center = TRUE, scale. = FALSE)
  stage_df <- as.data.frame(stage_pca$x[, 1:2])
  stage_df$sample <- rownames(stage_df)
  stage_df <- stage_df %>%
    left_join(qc_metrics, by = "sample")

  stage_pvar <- (stage_pca$sdev^2) / sum(stage_pca$sdev^2)
  stage_plot <- ggplot2::ggplot(
    stage_df,
    ggplot2::aes(x = PC1, y = PC2, color = exposure)
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      title = paste0("PCA (", stage, ", logCPM top 2000)"),
      x = paste0("PC1 (", scales::percent(stage_pvar[1]), ")"),
      y = paste0("PC2 (", scales::percent(stage_pvar[2]), ")")
    ) +
    ggplot2::theme_minimal(base_size = 11)

  stage_slug <- gsub("[^0-9A-Za-z]+", "_", stage)
  save_plot(stage_plot, paste0("pca_logcpm_", stage_slug, ".pdf"))
}

## ------------------------------------------------------------------
## 8. Análisis diferencial (DESeq2 ~ testis_id + exposure)
## ------------------------------------------------------------------

cat("Ejecutando DESeq2 completo (diseño ~ testis_id + exposure)...\n")

sample_sheet$exposure <- droplevels(sample_sheet$exposure)
sample_sheet$exposure <- stats::relevel(sample_sheet$exposure, ref = "DMSO")

deseq_counts <- round(expr_filtered)
dds_full <- DESeqDataSetFromMatrix(
  countData = deseq_counts,
  colData   = sample_sheet,
  design    = ~ testis_id + exposure
)

dds_full <- DESeq(dds_full)
saveRDS(dds_full, file.path(paths$results_dir, paste0(geo_id, "_dds_full.rds")))

contrasts <- setdiff(levels(sample_sheet$exposure), "DMSO")
deg_summary <- list()

for (trt in contrasts) {
  res <- results(dds_full, contrast = c("exposure", trt, "DMSO"), alpha = 0.05)
  res_tbl <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id") %>%
    arrange(padj)

  out_file <- file.path(deg_dir, paste0("DESeq2_", trt, "_vs_DMSO.tsv"))
  readr::write_tsv(res_tbl, out_file)
  cat("Resultados guardados:", out_file, "\n")

  ma_plot <- ggplot2::ggplot(res_tbl, ggplot2::aes(x = baseMean, y = log2FoldChange)) +
    ggplot2::geom_point(
      aes(color = ifelse(!is.na(padj) & padj < 0.05, "sig", "ns")),
      alpha = 0.5, size = 1
    ) +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::scale_color_manual(values = c("ns" = "grey70", "sig" = "firebrick"),
                                name = "padj < 0.05") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(
      title = paste0("MA plot: ", trt, " vs DMSO"),
      x = "baseMean (log10)",
      y = "log2 Fold Change"
    ) +
    ggplot2::theme_minimal(base_size = 11)
  save_plot(ma_plot, paste0("MAplot_", trt, "_vs_DMSO.pdf"), width = 7, height = 5)

  sig_tbl <- res_tbl %>%
    filter(!is.na(padj) & padj < 0.05)

  deg_summary[[trt]] <- tibble::tibble(
    contrast     = paste0(trt, " vs DMSO"),
    total_genes  = nrow(res_tbl),
    significant  = nrow(sig_tbl),
    upregulated   = sum(sig_tbl$log2FoldChange > 0, na.rm = TRUE),
    downregulated = sum(sig_tbl$log2FoldChange < 0, na.rm = TRUE)
  )
}

if (length(deg_summary)) {
  deg_summary_tbl <- dplyr::bind_rows(deg_summary)
  readr::write_tsv(deg_summary_tbl, file.path(deg_dir, "DESeq2_summary.tsv"))
  cat("\nResumen DESeq2 (contrastes vs DMSO):\n")
  print(deg_summary_tbl)
  cat("\n")
} else {
  cat("No se generaron contrastes DESeq2 (¿faltan tratamientos distintos de DMSO?).\n\n")
}

## Regenerar MA plots desde los TSV guardados ------------------------
ma_files <- list.files(deg_dir, pattern = "^DESeq2_.*_vs_DMSO\\.tsv$", full.names = TRUE)
if (length(ma_files)) {
  cat("Generando MA plots a partir de los TSV en results (", length(ma_files), " archivos )...\n", sep = "")
  deg_summary_from_tsv <- list()
  for (file in ma_files) {
    trt_name <- gsub("^DESeq2_|_vs_DMSO\\.tsv$", "", basename(file))
    tbl <- readr::read_tsv(file, show_col_types = FALSE)
    if (!all(c("baseMean", "log2FoldChange", "padj") %in% colnames(tbl))) {
      cat("Saltando", file, "- faltan columnas baseMean/log2FoldChange/padj\n")
      next
    }
    sig_tbl <- tbl %>% filter(!is.na(padj) & padj < 0.05)
    deg_summary_from_tsv[[trt_name]] <- tibble::tibble(
      contrast     = paste0(trt_name, " vs DMSO (TSV)"),
      total_genes  = nrow(tbl),
      significant  = nrow(sig_tbl),
      upregulated   = sum(sig_tbl$log2FoldChange > 0, na.rm = TRUE),
      downregulated = sum(sig_tbl$log2FoldChange < 0, na.rm = TRUE)
    )
    ma_plot <- ggplot2::ggplot(tbl,
                               ggplot2::aes(x = baseMean, y = log2FoldChange)) +
      ggplot2::geom_point(
        ggplot2::aes(color = ifelse(!is.na(padj) & padj < 0.05, "sig", "ns")),
        alpha = 0.5, size = 1
      ) +
      ggplot2::scale_x_continuous(trans = "log10") +
      ggplot2::scale_color_manual(values = c("ns" = "grey70", "sig" = "firebrick"),
                                  name = "padj < 0.05") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::labs(
        title = paste0("MA plot (TSV): ", trt_name, " vs DMSO"),
        x = "baseMean (log10)",
        y = "log2 Fold Change"
      ) +
      ggplot2::theme_minimal(base_size = 11)
    save_plot(ma_plot, paste0("MAplot_", trt_name, "_vs_DMSO_from_tsv.pdf"), width = 7, height = 5)
  }
  if (length(deg_summary_from_tsv)) {
    deg_summary_tsv_tbl <- dplyr::bind_rows(deg_summary_from_tsv)
    readr::write_tsv(deg_summary_tsv_tbl, file.path(deg_dir, "DESeq2_summary_from_tsv.tsv"))
    cat("\nResumen (recalculado desde TSV):\n")
    print(deg_summary_tsv_tbl)
    cat("\n")
  }
  cat("MA plots generados desde TSV.\n\n")
}

## ------------------------------------------------------------------
## 9. Enriquecimiento funcional con clusterProfiler
## ------------------------------------------------------------------

convert_ids <- function(df) {
  # intentamos mapear RefSeq primero, si falla probamos con SYMBOL
  ids_refseq <- tryCatch(
    bitr(df$gene_id, fromType = "REFSEQ",
         toType = c("ENTREZID", "SYMBOL"),
         OrgDb = org.Hs.eg.db),
    error = function(e) NULL
  )
  if (!is.null(ids_refseq) && nrow(ids_refseq) > 0) {
    return(
      df %>%
        inner_join(ids_refseq, by = c("gene_id" = "REFSEQ")) %>%
        distinct(ENTREZID, .keep_all = TRUE)
    )
  }

  ids_symbol <- tryCatch(
    bitr(df$gene_id, fromType = "SYMBOL",
         toType = c("ENTREZID", "SYMBOL"),
         OrgDb = org.Hs.eg.db),
    error = function(e) NULL
  )
  if (!is.null(ids_symbol) && nrow(ids_symbol) > 0) {
    return(
      df %>%
        inner_join(ids_symbol, by = c("gene_id" = "SYMBOL")) %>%
        distinct(ENTREZID, .keep_all = TRUE)
    )
  }

  tibble()
}

ego_dir <- file.path(paths$results_dir, "enrichment")
if (!dir.exists(ego_dir)) dir.create(ego_dir, recursive = TRUE)

for (trt in contrasts) {
  res_file <- file.path(deg_dir, paste0("DESeq2_", trt, "_vs_DMSO.tsv"))
  if (!file.exists(res_file)) next
  res_tbl <- readr::read_tsv(res_file, show_col_types = FALSE)
  if (!all(c("gene_id", "log2FoldChange", "padj") %in% colnames(res_tbl))) next

  converted <- convert_ids(res_tbl)
  universe_genes <- converted$ENTREZID
  sig_genes <- converted %>%
    filter(!is.na(padj) & padj < 0.05) %>%
    pull(ENTREZID)

  if (length(sig_genes) < 5) {
    cat("Saltando GO/KEGG para", trt, "- menos de 5 genes significativos tras mapear.\n")
    next
  }

  ego <- enrichGO(
    gene          = sig_genes,
    universe      = universe_genes,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )

  ekegg <- enrichKEGG(
    gene          = sig_genes,
    universe      = universe_genes,
    organism      = "hsa",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05
  )

  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    go_file <- file.path(ego_dir, paste0("enrichGO_", trt, "_vs_DMSO.tsv"))
    readr::write_tsv(as.data.frame(ego), go_file)
    cat("GO guardado:", go_file, "\n")
  } else {
    cat("Sin términos GO significativos para", trt, "\n")
  }

  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    kegg_file <- file.path(ego_dir, paste0("enrichKEGG_", trt, "_vs_DMSO.tsv"))
    readr::write_tsv(as.data.frame(ekegg), kegg_file)
    cat("KEGG guardado:", kegg_file, "\n")
  } else {
    cat("Sin términos KEGG significativos para", trt, "\n")
  }
}

## ------------------------------------------------------------------
## 10. Guardar objetos para análisis posteriores
## ------------------------------------------------------------------

saveRDS(
  list(
    counts_raw       = expr_mat,
    counts_filtered  = expr_filtered,
    log_cpm          = log_cpm,
    sample_sheet     = sample_sheet,
    qc_metrics       = qc_metrics,
    vsd_mat          = vsd_mat,
    dds_full         = dds_full
  ),
  file.path(paths$results_dir, paste0(geo_id, "_qc_objects.rds"))
)
cat("Objetos de QC guardados en ", file.path(paths$results_dir, paste0(geo_id, "_qc_objects.rds")), "\n")
