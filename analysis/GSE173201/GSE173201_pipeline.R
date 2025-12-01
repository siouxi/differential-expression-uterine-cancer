############################################################
## GSE173201 – QC y análisis de DEGs con DESeq2
## Este script:
##  1) Descarga la matriz de counts de GEO y la guarda en DATA/GSE173201
##  2) Realiza QC de counts crudos y normalizados
##  3) Ajusta un modelo DESeq2 (~ group)
##  4) Obtiene DEGs para contrastes biológicamente relevantes
## Todas las decisiones metodológicas se explican en los comentarios.
############################################################

############################################################
## 0. Librerías
############################################################

library(dplyr)     # Para manipulación ligera de datos/tablas
library(data.table) # fread rápido para tablas grandes
library(DESeq2)    # Normalización y DEGs
library(vsn)       # QC de media-desviación tras transformación
library(pheatmap)  # Heatmaps de correlación y distancias
library(corrplot)  # Visualización de matrices de correlación
library(fs)        # Manejo de directorios (como en otros pipelines)

############################################################
## 0.1 Directorios y configuración de proyecto (estilo GSE179661)
############################################################

geo_id <- "GSE173201"

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
  # Fallback: directorio de trabajo actual
  script_dir <- getwd()
}

base_dir <- dirname(dirname(script_dir))

paths <- list(
  data_dir    = file.path(base_dir, "DATA",   geo_id),
  results_dir = file.path(base_dir, "results", geo_id),
  qc_dir      = file.path(base_dir, "results", geo_id, "qc"),
  deseq_dir   = file.path(base_dir, "results", geo_id, "deseq2")
)

deseq_dir   <- paths$deseq_dir
qc_dir      <- paths$qc_dir
results_dir <- paths$results_dir

fs::dir_create(unlist(paths))

cat("Rutas del proyecto configuradas (GSE173201):\n")
cat("  Directorio base:", base_dir, "\n")
cat("  Directorio de datos:", paths$data_dir, "\n")
cat("  Directorio de resultados:", paths$results_dir, "\n\n")

############################################################
## 1. Descarga y carga de datos (counts de GEO)
############################################################

# URL de descarga directa de la matriz de counts de GEO
geo_base <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
geo_url  <- paste(
  geo_base,
  "acc=GSE173201",
  "file=GSE173201_raw_counts_GRCh38.p13_NCBI.tsv.gz",
  sep = "&"
)

counts_file <- file.path(paths$data_dir, "GSE173201_raw_counts_GRCh38.p13_NCBI.tsv.gz")

# Descargamos el archivo solo si no existe, para no repetir trabajo
if (!file.exists(counts_file)) {
  download.file(url = geo_url, destfile = counts_file, mode = "wb")
}

# Cargamos la tabla de counts como matriz (genes x muestras)
tbl <- as.matrix(
  fread(counts_file, header = TRUE, colClasses = "integer"),
  rownames = 1
)

############################################################
## 2. Metadatos de muestras
############################################################

# Decisión: definimos explícitamente el diseño biológico de GSE173201
#  - A2780: línea sensible
#  - A2780cis: línea resistente
#  - Control vs Cisplatino (tratamiento)
# Esto nos permitirá luego definir contrastes claros.

rna_names <- c(
  "RNA_S_C1", "RNA_S_C2", "RNA_S_C3",         # A2780 Control
  "RNA_R_C1", "RNA_R_C2", "RNA_R_C3",         # A2780cis Control
  "RNA_S_T1", "RNA_S_T2", "RNA_S_T3",         # A2780 + Cisplatino
  "RNA_R_T1", "RNA_R_T2", "RNA_R_T3"          # A2780cis + Cisplatino
)

colnames(tbl) <- rna_names

sample_info <- data.frame(
  sample_id = rna_names,
  group = factor(
    c(
      "A2780_Control", "A2780_Control", "A2780_Control",
      "A2780cis_Control", "A2780cis_Control", "A2780cis_Control",
      "A2780_Cisplatino", "A2780_Cisplatino", "A2780_Cisplatino",
      "A2780cis_Cisplatino", "A2780cis_Cisplatino", "A2780cis_Cisplatino"
    ),
    # Decisión: fijamos explícitamente el orden de niveles;
    # el baseline de referencia será A2780_Control.
    levels = c(
      "A2780_Control",
      "A2780cis_Control",
      "A2780_Cisplatino",
      "A2780cis_Cisplatino"
    )
  )
)
rownames(sample_info) <- sample_info$sample_id

group_colors <- c(
  "A2780_Control"        = "#A3C4F3",
  "A2780cis_Control"     = "#F9C74F",
  "A2780_Cisplatino"     = "#90BE6D",
  "A2780cis_Cisplatino"  = "#F8961E"
)

############################################################
## 3. Filtrado básico de genes
############################################################

# Decisión: eliminamos genes con muy baja expresión global.
# Criterio usado: al menos 2 muestras con >= 10 counts.
# Esto reduce ruido y acelera DESeq2 sin sesgar los genes realmente expresados.

keep <- rowSums(tbl >= 10) >= 2
tbl  <- tbl[keep, ]

############################################################
## 4. QC con counts crudos
############################################################

# Ordenamos muestras por grupo solo para plots (no afecta análisis)
sample_order <- order(sample_info$group)
tbl           <- tbl[, sample_order]
sample_info   <- sample_info[sample_order, ]

## 4.1 Boxplot log2(counts crudos + 1)
par(mfrow = c(1, 1))
boxplot(
  log2(tbl + 1),
  boxwex   = 0.7,
  notch    = TRUE,
  outline  = FALSE,
  las      = 2,
  main     = "Distribución de counts crudos (log2)",
  col      = group_colors[sample_info$group],
  ylab     = "log2(count + 1)"
)
legend(
  "topright",
  legend = names(group_colors),
  fill   = group_colors,
  border = "gray30",
  cex    = 0.8
)

## 4.2 Densidades de counts crudos (log10)
plot(
  density(log10(tbl[, 1] + 1)),
  type = "n",
  main = "Densidad de counts crudos (log10)",
  xlab = "log10(count + 1)",
  ylab = "Densidad"
)
for (i in seq_len(ncol(tbl))) {
  lines(
    density(log10(tbl[, i] + 1)),
    col = group_colors[sample_info$group[i]],
    lwd = 2
  )
}
legend(
  "topright",
  legend = names(group_colors),
  fill   = group_colors,
  cex    = 0.8
)

## 4.3 Correlación entre muestras (counts crudos log2)
cor_matrix_raw <- cor(log2(tbl + 1))

annotation_df <- data.frame(
  Grupo = sample_info$group
)
rownames(annotation_df) <- rownames(sample_info)

pheatmap(
  cor_matrix_raw,
  annotation_col   = annotation_df,
  annotation_colors = list(Grupo = group_colors),
  show_rownames    = TRUE,
  show_colnames    = TRUE,
  cluster_rows     = TRUE,
  cluster_cols     = TRUE,
  color            = colorRampPalette(c("#90BE6D", "white", "#A3C4F3"))(100),
  main             = "Heatmap de correlación (counts crudos)",
  display_numbers  = FALSE,
  fontsize         = 10,
  fontsize_row     = 8,
  fontsize_col     = 8
)

############################################################
## 5. Construcción del objeto DESeq2 y normalización
############################################################

# Decisión: modelo simple ~ group, suficiente dado el diseño balanceado
# y número de muestras. Para posteriores análisis se podría descomponer
# en ~ cell_line + treatment + cell_line:treatment.

col_data <- data.frame(
  group = sample_info$group
)
rownames(col_data) <- rownames(sample_info)

dds <- DESeqDataSetFromMatrix(
  countData = tbl,
  colData   = col_data,
  design    = ~ group
)

# Prefiltro adicional recomendado por DESeq2:
# genes con suma total de counts < 10 suelen ser ruido técnico.
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Ejecutamos pipeline completo de DESeq2 (incluye estimación de size factors)
dds <- DESeq(dds)

size_factors <- sizeFactors(dds)
round(size_factors, 3)

normalized_counts <- counts(dds, normalized = TRUE)

############################################################
## 6. QC post-normalización
############################################################

## 6.1 Boxplot de counts normalizados
boxplot(
  normalized_counts,
  boxwex  = 0.7,
  notch   = TRUE,
  outline = FALSE,
  las     = 2,
  main    = "Distribución de conteos normalizados",
  col     = group_colors[sample_info$group],
  ylab    = "Counts normalizados (size factors)"
)
legend(
  "topright",
  legend = names(group_colors),
  fill   = group_colors,
  border = "gray30",
  cex    = 0.8
)
## 8. Guardado de resultados
############################################################

# Decisión: guardamos tanto los resultados completos como las listas
# filtradas de DEGs, además de los counts normalizados y los gráficos,
# para usar en visualizaciones o análisis posteriores (en otros scripts).

results_dir <- paths$results_dir
qc_dir <- paths$qc_dir
fs::dir_create(c(results_dir, qc_dir))

############################################################
## 8.1 Guardar gráficos en PDF
############################################################

## Boxplot counts crudos
pdf(file.path(qc_dir, "01_boxplot_counts_crudos_log2.pdf"), width = 8, height = 6)
par(mfrow = c(1, 1))
boxplot(
  log2(tbl + 1),
  boxwex   = 0.7,
  notch    = TRUE,
  outline  = FALSE,
  las      = 2,
  main     = "Distribución de counts crudos (log2)",
  col      = group_colors[sample_info$group],
  ylab     = "log2(count + 1)"
)
legend(
  "topright",
  legend = names(group_colors),
  fill   = group_colors,
  border = "gray30",
  cex    = 0.8
)
dev.off()

## Densidad counts crudos
pdf(file.path(qc_dir, "02_densidad_counts_crudos_log10.pdf"), width = 8, height = 6)
plot(
  density(log10(tbl[, 1] + 1)),
  type = "n",
  main = "Densidad de counts crudos (log10)",
  xlab = "log10(count + 1)",
  ylab = "Densidad"
)
for (i in seq_len(ncol(tbl))) {
  lines(
    density(log10(tbl[, i] + 1)),
    col = group_colors[sample_info$group[i]],
    lwd = 2
  )
}
legend(
  "topright",
  legend = names(group_colors),
  fill   = group_colors,
  cex    = 0.8
)
dev.off()

## Heatmap correlación counts crudos
pheatmap(
  cor_matrix_raw,
  annotation_col    = annotation_df,
  annotation_colors = list(Grupo = group_colors),
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  color             = colorRampPalette(c("#90BE6D", "white", "#A3C4F3"))(100),
  main              = "Heatmap de correlación (counts crudos)",
  display_numbers   = FALSE,
  fontsize          = 10,
  fontsize_row      = 8,
  fontsize_col      = 8,
  filename          = file.path(qc_dir, "03_heatmap_correlacion_counts_crudos.pdf")
)

## Boxplot counts normalizados
pdf(file.path(qc_dir, "04_boxplot_counts_normalizados.pdf"), width = 8, height = 6)
boxplot(
  normalized_counts,
  boxwex  = 0.7,
  notch   = TRUE,
  outline = FALSE,
  las     = 2,
  main    = "Distribución de conteos normalizados",
  col     = group_colors[sample_info$group],
  ylab    = "Counts normalizados (size factors)"
)
legend(
  "topright",
  legend = names(group_colors),
  fill   = group_colors,
  border = "gray30",
  cex    = 0.8
)
dev.off()

## meanSdPlot VST
pdf(file.path(qc_dir, "05_meanSdPlot_VST.pdf"), width = 8, height = 6)
meanSdPlot(vst_mat)
dev.off()

## Heatmap correlación VST
pheatmap(
  cor_matrix_vst,
  annotation_col    = annotation_df,
  annotation_colors = list(Grupo = group_colors),
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  color             = colorRampPalette(c("#F94144", "#F3722C", "#F9C74F",
                                         "#90BE6D", "#43AA8B", "#277DA1",
                                         "#577590", "#A3C4F3"))(100),
  main              = "Heatmap de correlación (VST)",
  display_numbers   = FALSE,
  fontsize          = 10,
  fontsize_row      = 8,
  fontsize_col      = 8,
  filename          = file.path(qc_dir, "06_heatmap_correlacion_VST.pdf")
)

## PCA
pdf(file.path(qc_dir, "07_PCA_VST.pdf"), width = 7, height = 6)
plotPCA(vsd, intgroup = "group")
dev.off()

############################################################
## 8.2 Guardar matrices y tablas de resultados
############################################################

## 6.2 Transformación de varianza estabilizada (VST) para QC
# Decisión: la VST hace que la varianza sea aproximadamente independiente
# de la media, lo que mejora PCA, distancias y heatmaps.

vsd     <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

meanSdPlot(vst_mat)

## 6.3 Correlación con datos transformados (VST)
cor_matrix_vst <- cor(vst_mat)

pheatmap(
  cor_matrix_vst,
  annotation_col    = annotation_df,
  annotation_colors = list(Grupo = group_colors),
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  color             = colorRampPalette(c("#F94144", "#F3722C", "#F9C74F",
                                         "#90BE6D", "#43AA8B", "#277DA1",
                                         "#577590", "#A3C4F3"))(100),
  main              = "Heatmap de correlación (VST)",
  display_numbers   = FALSE,
  fontsize          = 10,
  fontsize_row      = 8,
  fontsize_col      = 8
)

## 6.4 PCA de muestras
plotPCA(vsd, intgroup = "group")

############################################################
## 7. Análisis de DEGs con DESeq2
############################################################

# Decisión: usamos alpha = 0.05 sobre padj (FDR de Benjamini–Hochberg)
# y |log2FC| >= 1 como umbral estándar para definir DEGs "fuertes".

alpha_level <- 0.05
lfc_cutoff  <- 1

# Contraste 1: efecto de cisplatino en línea sensible
res_A2780_Cis_vs_Ctrl <- results(
  dds,
  contrast = c("group", "A2780_Cisplatino", "A2780_Control"),
  alpha    = alpha_level
)

# Contraste 2: efecto de resistencia basal (A2780cis vs A2780, sin droga)
res_A2780cis_Ctrl_vs_A2780_Ctrl <- results(
  dds,
  contrast = c("group", "A2780cis_Control", "A2780_Control"),
  alpha    = alpha_level
)

# Contraste 3: efecto de cisplatino en línea resistente
res_A2780cis_Cis_vs_A2780cis_Ctrl <- results(
  dds,
  contrast = c("group", "A2780cis_Cisplatino", "A2780cis_Control"),
  alpha    = alpha_level
)

# Resúmenes rápidos
summary(res_A2780_Cis_vs_Ctrl)
summary(res_A2780cis_Ctrl_vs_A2780_Ctrl)
summary(res_A2780cis_Cis_vs_A2780cis_Ctrl)

## Contrastes adicionales para MA plots:
##  - Control sensible vs todas las demás condiciones
##  - Control resistente (A2780cis_Control) vs todas las condiciones

# A2780cis + Cisplatino vs A2780 Control
res_A2780cis_Cis_vs_A2780_Ctrl <- results(
  dds,
  contrast = c("group", "A2780cis_Cisplatino", "A2780_Control"),
  alpha    = alpha_level
)

# A2780cis Control vs A2780 Cisplatino
res_A2780cis_Ctrl_vs_A2780_Cis <- results(
  dds,
  contrast = c("group", "A2780cis_Control", "A2780_Cisplatino"),
  alpha    = alpha_level
)

# A2780cis Control vs A2780cis Cisplatino
res_A2780cis_Ctrl_vs_A2780cis_Cis <- results(
  dds,
  contrast = c("group", "A2780cis_Control", "A2780cis_Cisplatino"),
  alpha    = alpha_level
)

# Definimos DEGs aplicando padj y |log2FC|
deg_A2780_Cis_vs_Ctrl <- subset(
  as.data.frame(res_A2780_Cis_vs_Ctrl),
  !is.na(padj) & padj < alpha_level & abs(log2FoldChange) >= lfc_cutoff
)

deg_A2780cis_Ctrl_vs_A2780_Ctrl <- subset(
  as.data.frame(res_A2780cis_Ctrl_vs_A2780_Ctrl),
  !is.na(padj) & padj < alpha_level & abs(log2FoldChange) >= lfc_cutoff
)

deg_A2780cis_Cis_vs_A2780cis_Ctrl <- subset(
  as.data.frame(res_A2780cis_Cis_vs_A2780cis_Ctrl),
  !is.na(padj) & padj < alpha_level & abs(log2FoldChange) >= lfc_cutoff
)

cat("Número de DEGs (padj <", alpha_level, ", |log2FC| >=", lfc_cutoff, "):\n")
cat("  A2780 + Cisplatino vs A2780 Control:", nrow(deg_A2780_Cis_vs_Ctrl), "\n")
cat("  A2780cis Control vs A2780 Control   :", nrow(deg_A2780cis_Ctrl_vs_A2780_Ctrl), "\n")
cat("  A2780cis + Cis vs A2780cis Control :", nrow(deg_A2780cis_Cis_vs_A2780cis_Ctrl), "\n")

############################################################
## 7.1 MA plots y conteos detallados de DEGs (Up / Down / Total)
############################################################

cat("\nResumen detallado de DEGs por contraste (Up / Down / Total):\n")

ma_contrasts <- list(
  A2780_Cis_vs_A2780_Ctrl           = res_A2780_Cis_vs_Ctrl,
  A2780cis_Ctrl_vs_A2780_Ctrl       = res_A2780cis_Ctrl_vs_A2780_Ctrl,
  A2780cis_Cis_vs_A2780cis_Ctrl     = res_A2780cis_Cis_vs_A2780cis_Ctrl,
  A2780cis_Cis_vs_A2780_Ctrl        = res_A2780cis_Cis_vs_A2780_Ctrl,
  A2780cis_Ctrl_vs_A2780_Cis        = res_A2780cis_Ctrl_vs_A2780_Cis,
  A2780cis_Ctrl_vs_A2780cis_Cis     = res_A2780cis_Ctrl_vs_A2780cis_Cis
)

for (cn in names(ma_contrasts)) {
  res_obj <- ma_contrasts[[cn]]
  res_df  <- as.data.frame(res_obj)
  
  # Clasificación Up / Down según umbrales usados en los DEGs
  res_df$status <- NA_character_
  res_df$status[!is.na(res_df$padj) & res_df$padj < alpha_level & res_df$log2FoldChange >= lfc_cutoff]  <- "Up"
  res_df$status[!is.na(res_df$padj) & res_df$padj < alpha_level & res_df$log2FoldChange <= -lfc_cutoff] <- "Down"
  res_df$status[is.na(res_df$status)] <- "NS"
  
  n_up    <- sum(res_df$status == "Up",   na.rm = TRUE)
  n_down  <- sum(res_df$status == "Down", na.rm = TRUE)
  n_total <- n_up + n_down
  
  cat("  ", cn, " -> Up:", n_up, " Down:", n_down, " Total:", n_total, "\n")
  
  # MA plot en PDF (estilo DESeq2, con líneas en ±log2FC cutoff)
  pdf(file.path(deseq_dir, paste0("MAplot_", cn, ".pdf")), width = 7, height = 5)
  DESeq2::plotMA(res_obj, ylim = c(-5, 5), main = "")
  abline(h = c(-lfc_cutoff, lfc_cutoff), col = "gray50", lty = 2)
  dev.off()
}


############################################################
## 8. Guardado de resultados
############################################################

# Decisión: guardamos tanto los resultados completos como las listas
# filtradas de DEGs, además de los counts normalizados, para usar en
# visualizaciones o análisis posteriores (en otros scripts).

deseq_dir <- paths$deseq_dir
fs::dir_create(deseq_dir)

## Counts normalizados
write.csv(
  normalized_counts,
  file = file.path(deseq_dir, "GSE173201_normalized_counts_DESeq2.csv"),
  quote = FALSE
)

## Resultados completos de DESeq2
write.csv(
  as.data.frame(res_A2780_Cis_vs_Ctrl),
  file = file.path(deseq_dir, "DESeq2_A2780_Cisplatino_vs_Control_all_genes.csv"),
  quote = FALSE
)
write.csv(
  as.data.frame(res_A2780cis_Ctrl_vs_A2780_Ctrl),
  file = file.path(deseq_dir, "DESeq2_A2780cis_Control_vs_A2780_Control_all_genes.csv"),
  quote = FALSE
)
write.csv(
  as.data.frame(res_A2780cis_Cis_vs_A2780cis_Ctrl),
  file = file.path(deseq_dir, "DESeq2_A2780cis_Cisplatino_vs_A2780cis_Control_all_genes.csv"),
  quote = FALSE
)

## Listas de DEGs filtradas
write.csv(
  deg_A2780_Cis_vs_Ctrl,
  file = file.path(deseq_dir, "DEGs_A2780_Cisplatino_vs_Control_FDR0.05_log2FC1.csv"),
  quote = FALSE
)
write.csv(
  deg_A2780cis_Ctrl_vs_A2780_Ctrl,
  file = file.path(deseq_dir, "DEGs_A2780cis_Control_vs_A2780_Control_FDR0.05_log2FC1.csv"),
  quote = FALSE
)
write.csv(
  deg_A2780cis_Cis_vs_A2780cis_Ctrl,
  file = file.path(deseq_dir, "DEGs_A2780cis_Cisplatino_vs_A2780cis_Control_FDR0.05_log2FC1.csv"),
  quote = FALSE
)

############################################################
## Fin de la pipeline GSE173201
############################################################

