required_pkgs <- c(
  "tidyverse", "fs",
  "BiocManager",
  "oligo",                 # lectura y RMA de Affymetrix (GPL6244, Human Gene 1.0 ST)
  "limma",                 # DE en microarrays
  "pheatmap", "RColorBrewer", "ggplot2",
  "matrixStats",           # rowVars para seleccionar genes más variables
  "clusterProfiler", "org.Hs.eg.db",
  "AnnotationDbi",         # anotación vía OrgDb / chip-db
  "hugene10sttranscriptcluster.db"  # anotación específica de GPL6244
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


geo_id <- "GSE57978"

## Diseño biológico de GSE57978 (según GEO / paper)
## - Células madre de glioblastoma (glioma stem cells) de 3 pacientes:
##   0609, CC4121, 4121T3
## - Tratamientos: Veh (etanol vehículo) vs CBD (2 uM, 48h)
## - Objetivo: estimar el efecto de CBD vs Veh ajustando por línea celular

## Detectar directorio base asumiendo que este script vive en analysis/GSE57978 ----
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
  results_dir = file.path(base_dir, "results", geo_id),
  qc_dir      = file.path(base_dir, "results", geo_id, "qc"),
  deg_dir     = file.path(base_dir, "results", geo_id, "deg")
)

fs::dir_create(unlist(paths))

cat("Rutas del proyecto configuradas para", geo_id, ":\n")
cat("  Directorio base:   ", base_dir, "\n")
cat("  Directorio de datos:   ", paths$data_dir, "\n")
cat("  Directorio de resultados:", paths$results_dir, "\n\n")


## 1. Cargar pheno_data ---------------------------------------------------------
cat("Cargando archivo de fenotipos (pheno_data)...\n")

pheno_file <- file.path(paths$data_dir, paste0(geo_id, "_pheno_data.csv"))

if (!file.exists(pheno_file)) {
  stop("No se encontró el archivo de fenotipos: ", pheno_file)
}

pheno_raw <- readr::read_csv(pheno_file, show_col_types = FALSE)

cat("  Número de muestras en pheno_data:", nrow(pheno_raw), "\n")
cat("  Columnas disponibles en pheno_data:\n")
print(colnames(pheno_raw))
cat("\n")

## Construir tabla de muestras con condición (Veh vs CBD) y línea celular -------
sample_sheet <- pheno_raw %>%
  dplyr::transmute(
    sample_title  = title,
    geo_accession = geo_accession,
    condition_raw = title,
    cell_line_raw = title
  ) %>%
  dplyr::mutate(
    condition = dplyr::if_else(
      grepl("CBD", condition_raw, ignore.case = TRUE),
      "CBD",
      "Veh"
    ),
    cell_line = sub("_.*$", "", cell_line_raw)
  ) %>%
  dplyr::select(sample_title, geo_accession, condition, cell_line)

sample_sheet <- sample_sheet %>%
  dplyr::mutate(
    condition = factor(condition, levels = c("Veh", "CBD")),
    cell_line = factor(cell_line)
  )

cat("  Condiciones detectadas:", paste(levels(sample_sheet$condition), collapse = ", "), "\n")
cat("  Líneas celulares:", paste(levels(sample_sheet$cell_line), collapse = ", "), "\n\n")


## 2. Leer archivos CEL y aplicar RMA (Affymetrix) ------------------------------
cat("Leyendo archivos CEL y aplicando RMA...\n")

cel_dir <- file.path(paths$data_dir, "CELS")
if (!dir.exists(cel_dir)) {
  stop("No se encontró el directorio de CELS: ", cel_dir)
}

cel_files <- list.files(cel_dir, pattern = "[.]CEL(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
if (!length(cel_files)) {
  stop("No se encontraron archivos CEL en: ", cel_dir)
}

cat("  Archivos CEL encontrados:\n")
print(basename(cel_files))
cat("\n")

## Leer CELs con oligo
raw_data <- oligo::read.celfiles(cel_files)

## Renombrar muestras usando el nombre de archivo sin extensión
sampleNames(raw_data) <- basename(cel_files) %>%
  sub("\\.CEL\\.gz$", "", ., ignore.case = TRUE) %>%
  sub("\\.CEL$", "", ., ignore.case = TRUE)

cel_sample_df <- tibble::tibble(
  sample_name   = sampleNames(raw_data),
  geo_accession = stringr::str_extract(sample_name, "GSM[0-9]+")
)

if (any(is.na(cel_sample_df$geo_accession))) {
  warning("Algunos archivos CEL no pudieron mapearse a un geo_accession mediante el patrón GSM*. Verifica los nombres.")
}

## Emparejar con pheno_data por geo_accession
sample_sheet <- sample_sheet %>%
  dplyr::left_join(cel_sample_df, by = "geo_accession") %>%
  dplyr::filter(!is.na(sample_name))

if (nrow(sample_sheet) == 0) {
  stop("No se pudieron emparejar muestras entre pheno_data y archivos CEL.")
}

## Ordenar raw_data según sample_sheet
ord <- match(sample_sheet$sample_name, sampleNames(raw_data))
if (any(is.na(ord))) {
  stop("No se pudo ordenar raw_data según sample_sheet; revisa los nombres de muestras.")
}
raw_data <- raw_data[, ord]

## Añadir phenoData al objeto ExpressionSet
pdata <- sample_sheet %>%
  dplyr::select(sample_name, geo_accession, condition, cell_line) %>%
  tibble::column_to_rownames("sample_name")

Biobase::pData(raw_data) <- as.data.frame(pdata)

cat("  Número de muestras en raw_data:", ncol(raw_data), "\n\n")

## Normalización RMA (background + normalización entre arrays) ------------------
eset <- oligo::rma(raw_data, target = "core")
expr_mat <- Biobase::exprs(eset)

cat("  Matriz de expresión normalizada (RMA) con dimensiones:",
    paste(dim(expr_mat), collapse = " x "), "\n\n")


## 3. QC exploratorio (sin MA plots) -------------------------------------------
cat("Realizando QC exploratorio (boxplots, PCA, correlación)...\n")

qc_dir <- paths$qc_dir
fs::dir_create(qc_dir)

## 3.1 Boxplot de intensidades por muestra
pdf(file.path(qc_dir, "boxplot_RMA_intensities_by_sample.pdf"),
    width = 10, height = 6)
par(mar = c(8, 4, 4, 2))
boxplot(expr_mat,
        las = 2,
        main = "Distribución de intensidades (RMA) por muestra",
        ylab = "log2 intensities",
        col = RColorBrewer::brewer.pal(8, "Set2")[Biobase::pData(eset)$condition])
legend("topright",
       legend = levels(Biobase::pData(eset)$condition),
       fill = RColorBrewer::brewer.pal(8, "Set2")[seq_along(levels(Biobase::pData(eset)$condition))],
       cex = 0.8)
dev.off()

## 3.2 PCA global (top 2000 sondas más variables)
top_var <- min(2000, nrow(expr_mat))
row_var <- matrixStats::rowVars(expr_mat)
top_genes <- order(row_var, decreasing = TRUE)[seq_len(top_var)]
expr_top <- expr_mat[top_genes, ]

pca_res <- prcomp(t(expr_top), scale. = TRUE)
pca_df <- data.frame(
  PC1       = pca_res$x[, 1],
  PC2       = pca_res$x[, 2],
  sample    = colnames(expr_top),
  condition = Biobase::pData(eset)$condition,
  cell_line = Biobase::pData(eset)$cell_line
)

p_pca <- ggplot2::ggplot(
  pca_df,
  ggplot2::aes(x = PC1, y = PC2, color = condition, shape = cell_line, label = sample)
) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_text(vjust = -0.7, size = 3) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "PCA (RMA, top 2000 sondas más variables)",
    x = "PC1",
    y = "PC2"
  )

ggplot2::ggsave(
  file.path(qc_dir, "PCA_RMA_top2000.pdf"),
  p_pca, width = 7, height = 5
)

## 3.2b PCA enfatizando tratamiento (residualizando por línea celular)
##      - Quitamos el efecto de "cell_line" con limma::removeBatchEffect
##      - Mantenemos el efecto de "condition" (Veh vs CBD)

expr_mat_treat <- limma::removeBatchEffect(
  expr_mat,
  batch  = Biobase::pData(eset)$cell_line,
  design = model.matrix(~ Biobase::pData(eset)$condition)
)

top_var_treat <- min(2000, nrow(expr_mat_treat))
row_var_treat <- matrixStats::rowVars(expr_mat_treat)
top_genes_treat <- order(row_var_treat, decreasing = TRUE)[seq_len(top_var_treat)]
expr_top_treat <- expr_mat_treat[top_genes_treat, ]

pca_res_treat <- prcomp(t(expr_top_treat), scale. = TRUE)
pca_df_treat <- data.frame(
  PC1       = pca_res_treat$x[, 1],
  PC2       = pca_res_treat$x[, 2],
  sample    = colnames(expr_top_treat),
  condition = Biobase::pData(eset)$condition
)

p_pca_treat <- ggplot2::ggplot(
  pca_df_treat,
  ggplot2::aes(x = PC1, y = PC2, color = condition, label = sample)
) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_text(vjust = -0.7, size = 3) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "PCA (RMA, top 2000, corregido por línea celular)",
    x = "PC1",
    y = "PC2"
  )

ggplot2::ggsave(
  file.path(qc_dir, "PCA_RMA_top2000_condition_residualized_by_cell_line.pdf"),
  p_pca_treat, width = 7, height = 5
)

## 3.3 Heatmap de correlación entre muestras
cor_mat <- stats::cor(expr_mat, method = "pearson")

pdf(file.path(qc_dir, "sample_correlation_heatmap.pdf"),
    width = 7, height = 7)
pheatmap::pheatmap(
  cor_mat,
  color = colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(100),
  main = "Correlación entre muestras (RMA)",
  annotation_col = data.frame(
    condition = Biobase::pData(eset)$condition,
    cell_line = Biobase::pData(eset)$cell_line,
    row.names = colnames(expr_mat)
  )
)
dev.off()

cat("  Archivos de QC guardados en:", qc_dir, "\n\n")


## 4. Análisis de expresión diferencial con limma (CBD vs Veh) -----------------
##    Para parecerse al análisis del paper, usamos un diseño simple ~ condition
##    tratando las 3 líneas como replicados biológicos del efecto CBD vs Veh.
cat("Realizando análisis de DEGs (limma, CBD vs Veh, diseño ~ condition)...\n")

design <- model.matrix(~ condition,
                       data = as.data.frame(Biobase::pData(eset)))

colnames(design)

fit <- limma::lmFit(expr_mat, design)
fit <- limma::eBayes(fit)

## El coeficiente de interés es la diferencia CBD vs Veh
coef_name <- "conditionCBD"
if (!coef_name %in% colnames(fit$coefficients)) {
  stop("No se encontró el coeficiente ", coef_name, " en el diseño limma.")
}

tt <- limma::topTable(
  fit,
  coef = coef_name,
  n = Inf,
  adjust.method = "BH"
)

deg_dir <- paths$deg_dir
fs::dir_create(deg_dir)

## Guardar tabla completa
readr::write_tsv(
  tt %>%
    tibble::rownames_to_column("probe_id"),
  file.path(deg_dir, "CBD_vs_Veh_ALL_probes.tsv")
)

## DEGs significativos (criterio alineado con el paper: FDR <= 0.05, |logFC| >= 1)
lfc_cut <- 1
fdr_cut <- 0.05

degs <- tt %>%
  dplyr::filter(!is.na(adj.P.Val),
                adj.P.Val <= fdr_cut,
                abs(logFC) >= lfc_cut)

degs_up <- degs %>% dplyr::filter(logFC > 0)
degs_dn <- degs %>% dplyr::filter(logFC < 0)

readr::write_tsv(
  degs %>% tibble::rownames_to_column("probe_id"),
  file.path(deg_dir, paste0("CBD_vs_Veh_DEGs_FDR", fdr_cut, "_logFC", lfc_cut, ".tsv"))
)

readr::write_tsv(
  degs_up %>% tibble::rownames_to_column("probe_id"),
  file.path(deg_dir, paste0("CBD_vs_Veh_UP_FDR", fdr_cut, "_logFC", lfc_cut, ".tsv"))
)

readr::write_tsv(
  degs_dn %>% tibble::rownames_to_column("probe_id"),
  file.path(deg_dir, paste0("CBD_vs_Veh_DOWN_FDR", fdr_cut, "_logFC", lfc_cut, ".tsv"))
)

cat("  DEGs totales (FDR <=", fdr_cut, ", |logFC| >=", lfc_cut, "):", nrow(degs), "\n")
cat("    Up-regulated (CBD > Veh):", nrow(degs_up), "\n")
cat("    Down-regulated (CBD < Veh):", nrow(degs_dn), "\n\n")


## 4b. Análisis de expresión diferencial por línea celular ----------------------
cat("Realizando análisis de DEGs por línea celular (modelo ~ condition en cada línea)...\n")

cell_lines <- levels(Biobase::pData(eset)$cell_line)
cl_deg_dir <- file.path(deg_dir, "by_cell_line")
fs::dir_create(cl_deg_dir)

for (cl in cell_lines) {
  cat("  Línea celular:", as.character(cl), "...\n")

  idx_cl <- Biobase::pData(eset)$cell_line == cl
  expr_cl <- expr_mat[, idx_cl, drop = FALSE]
  pdata_cl <- Biobase::pData(eset)[idx_cl, , drop = FALSE]

  ## Necesitamos al menos una muestra Veh y una CBD
  if (length(unique(pdata_cl$condition)) < 2) {
    cat("    [ADVERTENCIA] Solo hay una condición en esta línea; se omite.\n")
    next
  }

  ## En cada línea solo hay 1 Veh y 1 CBD, así que no podemos estimar varianza ni p-valores
  ## con limma (df residuales = 0). Hacemos un análisis descriptivo: logFC = CBD - Veh.

  veh_idx <- which(pdata_cl$condition == "Veh")
  cbd_idx <- which(pdata_cl$condition == "CBD")

  if (length(veh_idx) != 1 || length(cbd_idx) != 1) {
    cat("    [ADVERTENCIA] Se esperaban 1 Veh y 1 CBD para la línea", as.character(cl), "; se omite.\n")
    next
  }

  logFC_cl   <- expr_cl[, cbd_idx, drop = FALSE] - expr_cl[, veh_idx, drop = FALSE]
  aveExpr_cl <- rowMeans(expr_cl)

  tt_cl <- data.frame(
    logFC      = as.numeric(logFC_cl),
    AveExpr    = aveExpr_cl,
    t          = NA_real_,
    P.Value    = NA_real_,
    adj.P.Val  = NA_real_
  )
  rownames(tt_cl) <- rownames(expr_cl)

  ## Guardar tabla completa por línea
  readr::write_tsv(
    tt_cl %>% tibble::rownames_to_column("probe_id"),
    file.path(cl_deg_dir, paste0("CBD_vs_Veh_", cl, "_ALL_probes.tsv"))
  )

  ## Criterio descriptivo por línea: solo por magnitud de logFC (sin FDR)
  degs_cl <- tt_cl %>%
    dplyr::filter(abs(logFC) >= lfc_cut)

  if (nrow(degs_cl) > 0) {
    # Up- y down-regulated dentro de la línea (descriptivo)
    degs_cl_up <- degs_cl %>% dplyr::filter(logFC > 0)
    degs_cl_dn <- degs_cl %>% dplyr::filter(logFC < 0)

    # Todos los genes con |logFC| >= lfc_cut
    readr::write_tsv(
      degs_cl %>% tibble::rownames_to_column("probe_id"),
      file.path(
        cl_deg_dir,
        paste0("CBD_vs_Veh_", cl, "_DEGs_logFC", lfc_cut, "_ALL.tsv")
      )
    )

    # Solo up-regulated (CBD > Veh)
    if (nrow(degs_cl_up) > 0) {
      readr::write_tsv(
        degs_cl_up %>% tibble::rownames_to_column("probe_id"),
        file.path(
          cl_deg_dir,
          paste0("CBD_vs_Veh_", cl, "_UP_logFC", lfc_cut, ".tsv")
        )
      )
    }

    # Solo down-regulated (CBD < Veh)
    if (nrow(degs_cl_dn) > 0) {
      readr::write_tsv(
        degs_cl_dn %>% tibble::rownames_to_column("probe_id"),
        file.path(
          cl_deg_dir,
          paste0("CBD_vs_Veh_", cl, "_DOWN_logFC", lfc_cut, ".tsv")
        )
      )
    }
  }

  cat("    [", as.character(cl), "] genes con |logFC| >=", lfc_cut,
      "(análisis descriptivo, sin FDR):", nrow(degs_cl), "\n\n")
}


## 5. Anotación de sondas a genes (chip GPL6244) ------------------------------
cat("Anotando sondas a genes humanos (hugene10sttranscriptcluster.db)...\n")

probe_ids <- rownames(expr_mat)

annot_tbl <- tryCatch({
  # Para GPL6244 (Human Gene 1.0 ST), las IDs de fila de expr_mat son PROBEID
  at <- AnnotationDbi::select(
    x       = hugene10sttranscriptcluster.db::hugene10sttranscriptcluster.db,
    keys    = probe_ids,
    columns = c("SYMBOL", "ENTREZID", "GENENAME"),
    keytype = "PROBEID"
  )

  # Nos quedamos con una fila por PROBEID para tener un mapping 1:1 sonda→gen
  at <- at %>%
    dplyr::distinct(PROBEID, .keep_all = TRUE) %>%
    dplyr::rename(
      probe_id  = PROBEID,
      symbol    = SYMBOL,
      entrez_id = ENTREZID,
      gene_name = GENENAME
    )

  at
}, error = function(e) {
  warning("Fallo en AnnotationDbi::select para PROBEID; devolviendo tabla vacía: ", e$message)
  tibble::tibble()
})

if (nrow(annot_tbl) > 0) {
  readr::write_tsv(
    annot_tbl,
    file.path(paths$results_dir, "probe_annotation_human_GPL6244.tsv")
  )
  cat("  Archivo de anotación guardado en:",
      file.path(paths$results_dir, "probe_annotation_human_GPL6244.tsv"), "\n\n")
} else {
  cat("  No se pudo generar anotación; revisa manualmente los IDs de sonda.\n\n")
}


cat("=== Pipeline de microarrays para GSE57978 completado (sin MA plots) ===\n")
cat("  - QC plots:", qc_dir, "\n")
cat("  - Tablas de DEGs:", deg_dir, "\n")


