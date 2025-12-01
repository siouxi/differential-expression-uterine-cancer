required_pkgs <- c(
  "tidyverse", "fs",
  "BiocManager", "GEOquery", "Biobase",
  "readr", "stringr", "tibble",
  "DESeq2", "matrixStats", "pheatmap", "RColorBrewer"
)

install_and_load <- function(pkgs) {
  new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(new_pkgs)) {
    message("Instalando paquetes faltantes: ", paste(new_pkgs, collapse = ", "))
    BiocManager::install(new_pkgs, ask = FALSE, update = FALSE)
  }
  suppressPackageStartupMessages(
    invisible(lapply(pkgs, library, character.only = TRUE))
  )
}

install_and_load(required_pkgs)


geo_id <- "GSE131565"

## Detectar directorio base asumiendo que este script vive en analysis/GSE131565 ----
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

if (is.null(script_dir) || identical(script_dir, "")) {
  script_dir <- getwd()
}

base_dir <- dirname(dirname(script_dir))

paths <- list(
  data_dir    = file.path(base_dir, "DATA", geo_id),
  results_dir = file.path(base_dir, "results", geo_id)
)

fs::dir_create(unlist(paths))

cat("Rutas del proyecto configuradas para", geo_id, ":\n")
cat("  Directorio base:   ", base_dir, "\n")
cat("  Directorio de datos:   ", paths$data_dir, "\n")
cat("  Directorio de resultados:", paths$results_dir, "\n\n")


## 1. Descargar y cargar datos desde GEO ----------------------------------------
cat("Descargando objeto GEO (ExpressionSet) para", geo_id, "...\n")

gset <- GEOquery::getGEO(
  geo_id,
  GSEMatrix = TRUE,
  AnnotGPL  = FALSE,
  destdir   = paths$data_dir
)

if (length(gset) > 1) {
  cat("  Advertencia: se encontraron", length(gset),
      "ExpressionSets. Usando el primero (gset[[1]]).\n")
}

eset <- gset[[1]]

## Metadatos de muestras (pheno) desde GEO
pheno_raw <- Biobase::pData(eset)
cat("Número de muestras en pheno_data:", nrow(pheno_raw), "\n")
cat("Columnas principales de pheno_data:\n")
print(colnames(pheno_raw))
cat("\n")


## 2. Guardar pheno_data en disco (para reutilizar sin re-descargar) -----------
cat("Guardando pheno_data en DATA/", geo_id, "...\n", sep = "")

pheno_file <- file.path(paths$data_dir, paste0(geo_id, "_pheno_data_from_GEO.csv"))
readr::write_csv(
  pheno_raw %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample_id"),
  pheno_file
)

cat("  pheno_data guardado en:      ", pheno_file, "\n\n")


## 3. Cargar matriz cruda de conteos desde GSE131565_raw_counts.txt ------------
cat("Cargando matriz cruda de conteos desde GSE131565_raw_counts.txt...\n")

counts_file <- file.path(paths$data_dir, "GSE131565_raw_counts.txt")

if (!file.exists(counts_file)) {
  stop("No se encontró el archivo de conteos crudos: ", counts_file)
}

raw_counts_df <- readr::read_tsv(
  counts_file,
  col_types = readr::cols(),
  progress = FALSE
)

cat("  Dimensiones del archivo crudo (filas x columnas): ",
    paste(dim(raw_counts_df), collapse = " x "), "\n")

# Primera columna = Geneid, columnas 7+ son las muestras (Control_1, ..., CBD_3)
expr_mat <- raw_counts_df %>%
  dplyr::select(Geneid, tidyselect::matches("Control_|CBD_")) %>%
  tibble::column_to_rownames("Geneid") %>%
  as.matrix()

cat("  Matriz de conteos construida con dimensiones (genes x muestras): ",
    paste(dim(expr_mat), collapse = " x "), "\n")
cat("  Nombres de muestras en expr_mat:\n")
print(colnames(expr_mat))
cat("\n")

## 4. Construir sample_sheet (Control vs CBD) ----------------------------------
cat("Construyendo tabla de muestras (sample_sheet)...\n")

sample_sheet <- tibble::tibble(
  sample    = colnames(expr_mat),
  condition = dplyr::case_when(
    grepl("^Control", sample, ignore.case = TRUE) ~ "Control",
    grepl("^CBD",     sample, ignore.case = TRUE) ~ "CBD",
    TRUE ~ NA_character_
  )
)

if (any(is.na(sample_sheet$condition))) {
  stop("No se pudo asignar condición (Control/CBD) a todas las muestras:\n",
       paste(sample_sheet$sample[is.na(sample_sheet$condition)], collapse = ", "))
}

sample_sheet <- sample_sheet %>%
  dplyr::mutate(
    condition = factor(condition, levels = c("Control", "CBD"))
  ) %>%
  tibble::column_to_rownames("sample")

stopifnot(identical(rownames(sample_sheet), colnames(expr_mat)))

cat("  Muestras y condiciones:\n")
print(sample_sheet)
cat("\n")


## 5. QC básico por muestra ----------------------------------------------------
qc_dir <- file.path(paths$results_dir, "qc")
deg_dir <- file.path(paths$results_dir, "deg")
fs::dir_create(c(qc_dir, deg_dir))

cat("Calculando métricas básicas de QC...\n")

library_size <- colSums(expr_mat)
detected_genes <- colSums(expr_mat > 0)
pct_zeros <- (1 - detected_genes / nrow(expr_mat)) * 100

qc_metrics <- tibble::tibble(
  sample        = colnames(expr_mat),
  library_size  = library_size,
  library_size_million = library_size / 1e6,
  detected_genes = detected_genes,
  pct_zeros     = pct_zeros,
  condition     = sample_sheet[colnames(expr_mat), "condition", drop = TRUE]
)

cat("Resumen numérico global de métricas QC:\n")
print(summary(qc_metrics[, c("library_size", "detected_genes", "pct_zeros")]))
cat("\n")

readr::write_tsv(qc_metrics, file.path(qc_dir, "qc_metrics.tsv"))
cat("Métricas de QC guardadas en qc_metrics.tsv\n\n")


## 6. DESeq2: Control vs CBD ---------------------------------------------------
cat("Ejecutando DESeq2 (diseño ~ condition, contraste CBD vs Control)...\n")

dds <- DESeqDataSetFromMatrix(
  countData = expr_mat,
  colData   = sample_sheet,
  design    = ~ condition
)

dds$condition <- droplevels(dds$condition)
dds$condition <- relevel(dds$condition, ref = "Control")

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "CBD", "Control"), alpha = 0.05)

res_tbl <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::arrange(padj)

out_file <- file.path(deg_dir, "DESeq2_CBD_vs_Control.tsv")
readr::write_tsv(res_tbl, out_file)
cat("Tabla completa de DESeq2 guardada en:", out_file, "\n")

sig_tbl <- res_tbl %>%
  dplyr::filter(!is.na(padj) & padj < 0.05)

cat("Resumen DESeq2 (CBD vs Control):\n")
cat("  Genes totales:      ", nrow(res_tbl), "\n")
cat("  Genes significativos (padj < 0.05): ", nrow(sig_tbl), "\n")
if (nrow(sig_tbl) > 0) {
  cat("    Up-regulated (CBD > Control):  ",
      sum(sig_tbl$log2FoldChange > 0, na.rm = TRUE), "\n")
  cat("    Down-regulated (CBD < Control):",
      sum(sig_tbl$log2FoldChange < 0, na.rm = TRUE), "\n")
}
cat("\n")

## MA plot simple --------------------------------------------------------------
ma_plot_file <- file.path(deg_dir, "MAplot_CBD_vs_Control_DESeq2.pdf")
pdf(ma_plot_file, width = 6, height = 5)
with(
  res,
  plot(
    x = log10(baseMean + 1),
    y = log2FoldChange,
    pch = 16, cex = 0.5,
    xlab = "log10(baseMean + 1)",
    ylab = "log2 Fold Change (CBD vs Control)",
    main = "MA plot DESeq2"
  )
)
abline(h = 0, col = "red", lty = 2)
dev.off()
cat("MA plot guardado en:", ma_plot_file, "\n\n")

cat("=== Pipeline de RNA-seq (QC básico + DESeq2) para", geo_id, "completada. ===\n")


