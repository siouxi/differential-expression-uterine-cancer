required_pkgs <- c(
  "tidyverse", "here", "jsonlite", "glue", "fs",
  "BiocManager", "GEOquery", "readr", "stringr",
  "limma", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi", "tibble",
  "pheatmap", "RColorBrewer", "ggplot2"
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


geo_id <- "GSE285498"

# Detectar directorio base (raíz del proyecto) asumiendo que este script vive en analysis/GSE285498
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
  # Fallback: use current working directory
  script_dir <- getwd()
}

base_dir <- dirname(dirname(script_dir))

paths <- list(
  data_dir = file.path(base_dir, "DATA", geo_id),
  results_dir = file.path(base_dir, "results", geo_id)
)

fs::dir_create(unlist(paths))

cat("Rutas del proyecto configuradas:\n")
cat("  Directorio base:", base_dir, "\n")
cat("  Directorio de datos:", paths$data_dir, "\n")
cat("  Directorio de resultados:", paths$results_dir, "\n\n")


## 1. Cargar pheno_data -------------------------------------------------------
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


## 2. Cargar matriz de conteos crudos ----------------------------------------
cat("Cargando matriz de conteos crudos (Raw_gene_counts)...\n")

counts_file <- file.path(paths$data_dir, paste0(geo_id, "_Raw_gene_counts.txt.gz"))

if (!file.exists(counts_file)) {
  stop("No se encontró el archivo de conteos crudos: ", counts_file)
}

## IMPORTANTE:
## El archivo está mal anotado en GEO: el primer renglón contiene SOLO los
## nombres de las muestras (CBD_1, CBD_2, ..., Mock_3), pero NO el nombre
## de la columna de genes. Las filas de datos tienen:
##   gene_id   CBD_1   CBD_2  ...  Mock_3
## Es decir, un campo extra al inicio. Si se lee de forma directa, el
## gene_id queda dentro de la columna CBD_1 y los valores se desplazan.

# Leemos SOLO la primera línea como texto para obtener los nombres de muestras
header_line <- readr::read_lines(counts_file, n_max = 1)
sample_names <- strsplit(header_line, "\t")[[1]]

cat("  Cabecera detectada en el archivo de conteos:\n")
print(sample_names)
cat("\n")

# Ahora leemos el resto del archivo (saltando la primera línea) y
# construimos explícitamente la columna gene_id + nombres de muestras.
counts_df <- readr::read_tsv(
  counts_file,
  skip = 1,
  col_names = c("gene_id", sample_names),
  show_col_types = FALSE
)

cat("  Dimensiones de la matriz de expresión (genes x muestras): ",
    paste(dim(counts_df), collapse = " x "), "\n")
cat("  Nombres de columnas (incluyendo gene_id):\n")
print(colnames(counts_df))
cat("\n")


## 3. Preparar matriz de expresión (log2(FPKM+1)) ----------------------------
cat("Preparando matriz de expresión (log2(FPKM+1))...\n")

# Extraemos ID Ensembl sin versión ni coordenadas:
#   ENSG00000290825.1_1:11868-14409 -> ENSG00000290825
counts_df <- counts_df %>%
  dplyr::mutate(
    ensembl_id = sub("\\..*$", "", sub("_.*$", "", gene_id))
  )

# Agrupamos por ensembl_id y promediamos las filas duplicadas
# (promedio de log2(FPKM+1) por gen y muestra)
counts_df <- counts_df %>%
  dplyr::group_by(ensembl_id) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(sample_names),
                  ~ mean(as.numeric(.x), na.rm = TRUE)),
    .groups = "drop"
  )

# Matriz completa (antes de filtrar)
expr_mat_raw <- counts_df %>%
  dplyr::select(ensembl_id, dplyr::all_of(sample_names)) %>%
  tibble::column_to_rownames("ensembl_id") %>%
  as.matrix()

cat("  Matriz de expresión creada con dimensiones (antes de filtrar):",
    paste(dim(expr_mat_raw), collapse = " x "), "\n")

# Filtro de expresión mínima:
#   - Mantener genes con log2(FPKM+1) >= 1 en al menos 3 muestras
min_expr <- 1       # ~ FPKM >= 1
min_samples <- 3    # al menos 3 de 12 muestras

keep_genes <- rowSums(expr_mat_raw >= min_expr, na.rm = TRUE) >= min_samples
expr_mat <- expr_mat_raw[keep_genes, , drop = FALSE]

cat("  Matriz de expresión después de filtrar genes poco expresados:",
    paste(dim(expr_mat), collapse = " x "), "\n\n")


## 4. Construir tabla de muestras (sample_sheet) -----------------------------
cat("Construyendo tabla de metadatos de muestras...\n")

sample_sheet <- pheno_raw %>%
  dplyr::mutate(
    sample = stringr::str_extract(title, "(Mock|CBD|Eto|Merge)_\\d+"),
    condition = dplyr::case_when(
      grepl("non-treatment", `treatment:ch1`, ignore.case = TRUE) ~ "Mock",
      grepl("Cannabidiol", `treatment:ch1`, ignore.case = TRUE) ~ "CBD",
      grepl("Etoposide treatment", `treatment:ch1`, ignore.case = TRUE) ~ "Eto",
      grepl("Combination treatment", `treatment:ch1`, ignore.case = TRUE) ~ "Merge",
      TRUE ~ "Unknown"
    ),
    replicate = as.integer(stringr::str_extract(sample, "\\d+"))
  ) %>%
  dplyr::select(sample, condition, replicate, geo_accession)

# Filtramos solo las muestras que están en la matriz de expresión y reordenamos
sample_sheet <- sample_sheet %>%
  dplyr::filter(sample %in% sample_names) %>%
  dplyr::arrange(match(sample, sample_names))

# Reordenamos columnas de las matrices de expresión según sample_sheet
expr_mat_raw <- expr_mat_raw[, sample_sheet$sample, drop = FALSE]
expr_mat <- expr_mat[, sample_sheet$sample, drop = FALSE]

# Convertimos condición a factor con Mock como referencia
sample_sheet <- sample_sheet %>%
  dplyr::mutate(
    condition = factor(condition, levels = c("Mock", "CBD", "Eto", "Merge"))
  )

cat("  Número de muestras en sample_sheet:", nrow(sample_sheet), "\n")
cat("  Condiciones:", paste(levels(sample_sheet$condition), collapse = ", "), "\n\n")

# Paleta de colores de QC (misma que en GSE179661)
qc_colors <- c(
  "#90BE6D",  # verde original
  "#A3C4F3",  # azul original
  "#F8961E",  # naranja original
  "#F9C74F",  # amarillo original
  "#577590",  # azul grisáceo profundo
  "#43AA8B",  # verde agua saturado
  "#F94144",  # rojo fuerte
  "#F3722C",  # naranja rojizo
  "#4D908E",  # verde petróleo
  "#277DA1"   # azul turquesa oscuro
)

sample_colors_qc <- rep(qc_colors, length.out = ncol(expr_mat))
names(sample_colors_qc) <- colnames(expr_mat)


## 5. QC exploratorio --------------------------------------------------------
cat("Realizando QC exploratorio...\n")

qc_dir <- file.path(paths$results_dir, "qc")
fs::dir_create(qc_dir)

# 5.1 Boxplots de distribuciones: antes y después del filtrado
pdf(file.path(qc_dir, "boxplot_log2FPKM_by_sample_BEFORE_filtering.pdf"),
    width = 10, height = 6)
par(mar = c(8, 4, 4, 2))
boxplot(expr_mat_raw,
        las = 2,
        main = "",
        ylab = "log2(FPKM + 1)",
        col = sample_colors_qc[colnames(expr_mat_raw)])
legend("topright",
       legend = levels(sample_sheet$condition),
       fill = RColorBrewer::brewer.pal(8, "Set2")[seq_along(levels(sample_sheet$condition))],
       cex = 0.8)
dev.off()

pdf(file.path(qc_dir, "boxplot_log2FPKM_by_sample_AFTER_filtering.pdf"),
    width = 10, height = 6)
par(mar = c(8, 4, 4, 2))
boxplot(expr_mat,
        las = 2,
        main = "",
        ylab = "log2(FPKM + 1)",
        col = sample_colors_qc[colnames(expr_mat)])
legend("topright",
       legend = levels(sample_sheet$condition),
       fill = RColorBrewer::brewer.pal(8, "Set2")[seq_along(levels(sample_sheet$condition))],
       cex = 0.8)
dev.off()

# 5.2 PCA
pca_res <- prcomp(t(expr_mat), scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  sample = rownames(pca_res$x),
  condition = sample_sheet$condition
)

p_pca <- ggplot2::ggplot(pca_df,
                         ggplot2::aes(x = PC1, y = PC2,
                                      color = condition, label = sample)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_text(vjust = -0.7, size = 3) +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "PC1", y = "PC2")

ggplot2::ggsave(file.path(qc_dir, "PCA_log2FPKM.pdf"),
                p_pca, width = 7, height = 5)

# 5.3 Heatmap de correlaciones entre muestras + corrplot.mixed
cor_mat <- stats::cor(expr_mat, method = "pearson", use = "pairwise.complete.obs")

readr::write_tsv(
  as.data.frame(cor_mat) %>% tibble::rownames_to_column("Sample"),
  file.path(qc_dir, "sample_correlation_matrix.tsv")
)

pdf(file.path(qc_dir, "sample_correlation_heatmap_custom_colors.pdf"),
    width = 7, height = 7)
pheatmap::pheatmap(
  cor_mat,
  color = colorRampPalette(qc_colors)(100),
  main = "",
  annotation_col = data.frame(condition = sample_sheet$condition,
                              row.names = sample_sheet$sample)
)
dev.off()

pdf(file.path(qc_dir, "sample_correlation_corrplot_mixed.pdf"),
    width = 7, height = 7)
corrplot::corrplot.mixed(
  cor_mat,
  lower = "number",
  upper = "circle",
  number.cex = 0.6,
  tl.cex = 0.6
)
dev.off()

# 5.4 Densidad de expresión por muestra (log2(FPKM+1))
pdf(file.path(qc_dir, "density_log2FPKM_by_sample.pdf"),
    width = 8, height = 6)
plotDensityExpr <- function(mat, colors) {
  colors_rep <- rep(colors, length.out = ncol(mat))
  plot(
    density(mat[, 1], na.rm = TRUE),
    main = "",
    xlab = "log2(FPKM + 1)",
    ylab = "Densidad",
    col = colors_rep[1],
    lwd = 2
  )
  if (ncol(mat) > 1) {
    for (i in 2:ncol(mat)) {
      lines(density(mat[, i], na.rm = TRUE),
            col = colors_rep[i],
            lwd = 2)
    }
  }
  legend("topright",
         legend = colnames(mat),
         col = colors_rep,
         lty = 1,
         cex = 0.6)
}
plotDensityExpr(expr_mat, sample_colors_qc[colnames(expr_mat)])
dev.off()

# 5.5 Mean–SD plots (ggplot + vsn::meanSdPlot)
gene_means <- rowMeans(expr_mat, na.rm = TRUE)
gene_sds <- apply(expr_mat, 1, sd, na.rm = TRUE)

sdmean_df <- data.frame(
  mean = gene_means,
  sd = gene_sds
)

sdmean_plot <- ggplot2::ggplot(sdmean_df, ggplot2::aes(x = mean, y = sd)) +
  ggplot2::geom_point(color = "#4D908E", alpha = 0.3, size = 0.7) +
  ggplot2::geom_smooth(color = "#F94144", se = FALSE, method = "loess") +
  ggplot2::labs(
    x = "Media de expresión (log2(FPKM+1))",
    y = "Desviación estándar"
  ) +
  ggplot2::theme_bw()

ggplot2::ggsave(
  file.path(qc_dir, "mean_sd_plot_log2FPKM_ggplot.pdf"),
  sdmean_plot,
  width = 7,
  height = 5
)

pdf(file.path(qc_dir, "meanSdPlot_log2FPKM_vsn.pdf"),
    width = 7, height = 5)
vsn::meanSdPlot(expr_mat)
dev.off()

cat("  Archivos de QC guardados en:", qc_dir, "\n\n")


## 6. Análisis de DEGs con limma ---------------------------------------------
cat("Realizando análisis de DEGs (limma, modelo lineal + eBayes)...\n")

design <- model.matrix(~ 0 + condition, data = sample_sheet)
colnames(design) <- levels(sample_sheet$condition)

fit <- limma::lmFit(expr_mat, design)

contrast_matrix <- limma::makeContrasts(
  CBD_vs_Mock   = CBD   - Mock,
  Eto_vs_Mock   = Eto   - Mock,
  Merge_vs_Mock = Merge - Mock,
  CBD_vs_Eto    = CBD   - Eto,
  levels = design
)

fit2 <- limma::contrasts.fit(fit, contrast_matrix)
fit2 <- limma::eBayes(fit2)

deg_dir <- file.path(paths$results_dir, "degs")
fs::dir_create(deg_dir)

save_deg_table <- function(fit_obj, coef_name, file_prefix,
                           lfc_cut = 1, fdr_cut = 0.05) {
  tt <- limma::topTable(
    fit_obj,
    coef = coef_name,
    n = Inf,
    adjust.method = "BH"
  )
  tt$ENSEMBL <- rownames(tt)
  tt$contrast <- coef_name

  # Añadimos dirección básica según logFC
  tt <- tt %>%
    dplyr::mutate(
      direction = dplyr::case_when(
        logFC > 0  ~ "Up",
        logFC < 0  ~ "Down",
        TRUE       ~ "NS"
      )
    )

  # Filtrado por DEGs significativos
  degs <- tt %>%
    dplyr::filter(!is.na(adj.P.Val),
                  adj.P.Val <= fdr_cut,
                  abs(logFC) >= lfc_cut)

  # Separar sobreexpresados (up) e infraexpresados (down)
  degs_up <- degs %>% dplyr::filter(logFC > 0)
  degs_down <- degs %>% dplyr::filter(logFC < 0)

  readr::write_tsv(
    tt,
    file.path(deg_dir,
              paste0(file_prefix, "_ALL_genes.tsv"))
  )
  readr::write_tsv(
    degs,
    file.path(deg_dir,
              paste0(file_prefix, "_DEGs_FDR", fdr_cut,
                     "_logFC", lfc_cut, ".tsv"))
  )
  readr::write_tsv(
    degs_up,
    file.path(deg_dir,
              paste0(file_prefix, "_UP_FDR", fdr_cut,
                     "_logFC", lfc_cut, ".tsv"))
  )
  readr::write_tsv(
    degs_down,
    file.path(deg_dir,
              paste0(file_prefix, "_DOWN_FDR", fdr_cut,
                     "_logFC", lfc_cut, ".tsv"))
  )

  # Anotar genes (solo IDs presentes en este contraste)
  annot_sub <- tryCatch({
    gm <- clusterProfiler::bitr(
      geneID  = tt$ENSEMBL,
      fromType = "ENSEMBL",
      toType   = c("SYMBOL", "ENTREZID", "GENENAME"),
      OrgDb    = org.Hs.eg.db::org.Hs.eg.db
    )

    gm %>%
      dplyr::distinct(ENSEMBL, .keep_all = TRUE)
  }, error = function(e) {
    AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys    = tt$ENSEMBL,
      keytype = "ENSEMBL",
      columns = c("SYMBOL", "ENTREZID", "GENENAME")
    ) %>%
      dplyr::distinct(ENSEMBL, .keep_all = TRUE)
  })

  tt_annot <- tt %>%
    dplyr::left_join(annot_sub, by = c("ENSEMBL" = "ENSEMBL")) %>%
    dplyr::rename(
      ENSEMBL_ID = ENSEMBL,
      SYMBOL     = SYMBOL,
      ENTREZID   = ENTREZID,
      GENENAME   = GENENAME
    )

  degs_annot    <- tt_annot %>% dplyr::filter(adj.P.Val <= fdr_cut,
                                             abs(logFC) >= lfc_cut)
  degs_up_annot <- degs_annot %>% dplyr::filter(direction == "Up")
  degs_dn_annot <- degs_annot %>% dplyr::filter(direction == "Down")

  readr::write_tsv(
    tt_annot,
    file.path(deg_dir,
              paste0(file_prefix, "_ALL_genes_ANNOTATED.tsv"))
  )
  readr::write_tsv(
    degs_annot,
    file.path(deg_dir,
              paste0(file_prefix, "_DEGs_FDR", fdr_cut,
                     "_logFC", lfc_cut, "_ANNOTATED.tsv"))
  )
  readr::write_tsv(
    degs_up_annot,
    file.path(deg_dir,
              paste0(file_prefix, "_UP_FDR", fdr_cut,
                     "_logFC", lfc_cut, "_ANNOTATED.tsv"))
  )
  readr::write_tsv(
    degs_dn_annot,
    file.path(deg_dir,
              paste0(file_prefix, "_DOWN_FDR", fdr_cut,
                     "_logFC", lfc_cut, "_ANNOTATED.tsv"))
  )

  cat("  Contraste", coef_name, "->",
      nrow(degs), "DEGs (FDR <=", fdr_cut,
      ", |logFC| >=", lfc_cut, ")",
      "[Up:", nrow(degs_up), ", Down:", nrow(degs_down), "]\n")

  invisible(tt)
}

res_CBD_vs_Mock   <- save_deg_table(fit2, "CBD_vs_Mock",   "CBD_vs_Mock")
res_Eto_vs_Mock   <- save_deg_table(fit2, "Eto_vs_Mock",   "Eto_vs_Mock")
res_Merge_vs_Mock <- save_deg_table(fit2, "Merge_vs_Mock", "Merge_vs_Mock")
res_CBD_vs_Eto    <- save_deg_table(fit2, "CBD_vs_Eto",    "CBD_vs_Eto")

cat("\n  Tablas de DEGs guardadas en:", deg_dir, "\n\n")


## 6.1. MA plots para cada contraste -----------------------------------------
cat("Generando MA plots para cada contraste...\n")

contrast_names <- colnames(contrast_matrix)

for (cn in contrast_names) {
  # Obtener tabla completa para este contraste
  tt_ma <- limma::topTable(
    fit2,
    coef = cn,
    n = Inf,
    adjust.method = "BH"
  )
  tt_ma$ENSEMBL <- rownames(tt_ma)

  # Clasificar genes en Up/Down/NS según los mismos umbrales de DEGs
  lfc_cut <- 1
  fdr_cut <- 0.05

  tt_ma <- tt_ma %>%
    dplyr::mutate(
      status = dplyr::case_when(
        !is.na(adj.P.Val) & adj.P.Val <= fdr_cut & logFC >= lfc_cut  ~ "Up",
        !is.na(adj.P.Val) & adj.P.Val <= fdr_cut & logFC <= -lfc_cut ~ "Down",
        TRUE                                                         ~ "NS"
      )
    )

  n_up <- sum(tt_ma$status == "Up")
  n_down <- sum(tt_ma$status == "Down")

  outfile <- file.path(deg_dir, paste0("MAplot_", cn, ".pdf"))
  pdf(outfile, width = 7, height = 5)

  p_ma <- ggplot2::ggplot(
    tt_ma,
    ggplot2::aes(x = AveExpr, y = logFC, color = status)
  ) +
    ggplot2::geom_point(alpha = 0.6, size = 1) +
    ggplot2::scale_color_manual(
      values = c(Up = "#E41A1C", Down = "#377EB8", NS = "grey70"),
      name = "Estado",
      breaks = c("Up", "Down", "NS"),
      labels = c(
        paste0("Up (", n_up, ")"),
        paste0("Down (", n_down, ")"),
        "No significativo"
      )
    ) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    ggplot2::labs(
      title = paste("MA plot -", cn),
      x = "AveExpr (nivel medio de expresión)",
      y = "log2 Fold Change"
    ) +
    ggplot2::theme_bw()

  print(p_ma)
  dev.off()

  cat("  MA plot guardado:", outfile,
      " [Up:", n_up, ", Down:", n_down, "]\n")
}

cat("\n")


## 7. Anotación de genes (usando clusterProfiler::bitr) ----------------------
cat("Anotando genes (clusterProfiler::bitr, Ensembl -> símbolo, Entrez, descripción)...\n")

ensembl_ids <- rownames(expr_mat)

annot_tbl <- tryCatch({
  gm <- clusterProfiler::bitr(
    geneID = ensembl_ids,
    fromType = "ENSEMBL",
    toType = c("SYMBOL", "ENTREZID", "GENENAME"),
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )

  gm %>%
    dplyr::distinct(ENSEMBL, .keep_all = TRUE) %>%
    dplyr::rename(
      ENSEMBL_ID = ENSEMBL,
      SYMBOL = SYMBOL,
      ENTREZID = ENTREZID,
      GENENAME = GENENAME
    )
}, error = function(e) {
  cat("  Advertencia: fallo en clusterProfiler::bitr (", e$message, ").\n",
      "  Usando AnnotationDbi::select como alternativa.\n", sep = "")

  AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = ensembl_ids,
    keytype = "ENSEMBL",
    columns = c("SYMBOL", "ENTREZID", "GENENAME")
  ) %>%
    dplyr::distinct(ENSEMBL, .keep_all = TRUE) %>%
    dplyr::rename(
      ENSEMBL_ID = ENSEMBL,
      SYMBOL = SYMBOL,
      ENTREZID = ENTREZID,
      GENENAME = GENENAME
    )
})

readr::write_tsv(
  annot_tbl,
  file.path(paths$results_dir, "gene_annotation_ensembl_human.tsv")
)

cat("  Archivo de anotación guardado en:",
    file.path(paths$results_dir, "gene_annotation_ensembl_human.tsv"), "\n\n")

cat("=== Pipeline de QC + DEGs + anotación para GSE285498 completado ===\n")




