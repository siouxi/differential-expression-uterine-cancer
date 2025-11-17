# Análisis del dataset GSE57978
# Affymetrix Gene 1.0 ST arrays
# Glioblastoma stem cells tratadas con CBD vs Vehicle

# ==============================================================================
# 1. Cargar librerías necesarias
# ==============================================================================
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Instalar paquetes de Bioconductor si no están instalados
required_packages <- c("oligo", "Biobase", "GEOquery", "arrayQualityMetrics", 
                       "simpleaffy", "affyPLM", "RColorBrewer", "gplots", 
                       "limma", "ggplot2", "gridExtra", "matrixStats",
                       "hugene10sttranscriptcluster.db", "AnnotationDbi")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) BiocManager::install(new_packages)
# Cargar librerías
library(oligo)
library(Biobase)
library(GEOquery)
library(arrayQualityMetrics)
library(simpleaffy)
library(affyPLM)
library(RColorBrewer)
library(gplots)
library(limma)
library(ggplot2)
library(gridExtra)
library(matrixStats)
library(hugene10sttranscriptcluster.db)
library(AnnotationDbi)

# ==============================================================================
# 2. Definir rutas de los archivos
# ==============================================================================
# Obtener la ruta base del proyecto desde la ubicación del script
# Intentar obtener la ruta del script si estamos en RStudio
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

# Si no se pudo obtener la ruta del script, usar alternativas
if (is.null(script_dir) || length(script_dir) == 0 || script_dir == "") {
    # Intentar desde getwd()
    current_dir <- getwd()
    if (basename(current_dir) == "GSE57978") {
        script_dir <- current_dir
    } else {
        # Intentar construir la ruta relativa
        script_dir <- file.path(current_dir, "ANALISIS", "GSE57978")
        if (!dir.exists(script_dir)) {
            # Última opción: asumir que estamos en el directorio raíz del proyecto
            script_dir <- file.path(current_dir, "ANALISIS", "GSE57978")
        }
    }
}

# Obtener la ruta base del proyecto (subir dos niveles desde ANALISIS/GSE57978)
base_dir <- dirname(dirname(script_dir))
data_dir <- file.path(base_dir, "DATA", "GSE57978")
cels_dir <- file.path(data_dir, "CELS")
pheno_file <- file.path(data_dir, "GSE57978_pheno_data.csv")

cat("Rutas configuradas:\n")
cat("  Base directory:", base_dir, "\n")
cat("  Data directory:", data_dir, "\n")
cat("  CELs directory:", cels_dir, "\n")
cat("  Pheno file:", pheno_file, "\n\n")

# Verificar que las rutas existen
if (!dir.exists(cels_dir)) {
    stop(paste("El directorio de CELs no existe:", cels_dir))
}
if (!file.exists(pheno_file)) {
    stop(paste("El archivo de phenodata no existe:", pheno_file))
}

# ==============================================================================
# 3. Cargar phenodata
# ==============================================================================
cat("Cargando phenodata...\n")
pheno_data <- read.csv(pheno_file, stringsAsFactors = FALSE)

# ==============================================================================
# 4. Cargar archivos CEL
# ==============================================================================

# Listar archivos CEL.gz en el directorio
cel_files <- list.files(cels_dir, pattern = "\\.CEL\\.gz$", 
                        full.names = TRUE, ignore.case = TRUE)

print(basename(cel_files))

raw_data <- read.celfiles(cel_files)

# Mostrar información del objeto
cat("\nInformación del objeto raw_data:\n")
print(raw_data)
cat("\nDimensiones:", dim(raw_data), "\n")
cat("Número de features:", nrow(raw_data), "\n")
cat("Número de muestras:", ncol(raw_data), "\n")

# ==============================================================================
# 5. Asociar phenodata con los datos de expresión
# ==============================================================================
cat("\nAsociando phenodata con los datos...\n")

# Extraer nombres de muestras del objeto raw_data
sample_names <- sampleNames(raw_data)
cat("Nombres de muestras en raw_data:\n")
print(sample_names)

gsm_ids <- gsub(".*(GSM\\d+).*", "\\1", basename(cel_files))
cat("\nGSM IDs extraídos:\n")
print(gsm_ids)

# Crear phenodata para el ExpressionSet
# Filtrar pheno_data para incluir solo las muestras presentes
pheno_subset <- pheno_data[match(gsm_ids, pheno_data$geo_accession), ]

# Verificar que todas las muestras tienen phenodata
if (any(is.na(pheno_subset$geo_accession))) {
    warning("Algunas muestras no tienen phenodata correspondiente")
}

# Crear AnnotatedDataFrame para el phenodata
pheno_annotated <- AnnotatedDataFrame(pheno_subset)
rownames(pheno_annotated) <- sample_names

# Asociar phenodata con raw_data
pData(raw_data) <- pheno_subset
rownames(pData(raw_data)) <- sample_names

cat("\nPhenodata asociado correctamente\n")
cat("\nResumen del phenodata asociado:\n")
print(pData(raw_data)[, c("geo_accession", "title", "source_name_ch1")])

# Definir directorio de salida
output_dir <- script_dir
if (!dir.exists(output_dir)) {
    output_dir <- getwd()
}

# ==============================================================================
# 6. Análisis de Calidad de Arrays (QC)
# Objetivo: Detectar arrays malos, outliers o problemas de hibridización
# ==============================================================================
# Crear directorio para guardar gráficos de QC
qc_dir <- file.path(output_dir, "QC")
if (!dir.exists(qc_dir)) {
    dir.create(qc_dir, recursive = TRUE)
}
cat("Directorio de QC creado:", qc_dir, "\n\n")

# ------------------------------------------------------------------------------
# 6.1. Distribución de Intensidades (Histogramas)
# ------------------------------------------------------------------------------
cat("6.1. Generando histogramas de distribución de intensidades...\n")

# Obtener intensidades sin normalizar
intensities <- log2(pm(raw_data) + 1)  # PM probes, log2 transformado

# Crear histograma
png(file.path(qc_dir, "01_histogram_intensities.png"), 
    width = 12, height = 8, units = "in", res = 300)

par(mfrow = c(2, 3))
for (i in 1:ncol(intensities)) {
    hist(intensities[, i], 
         main = paste("Distribución de Intensidades\n", sample_names[i]),
         xlab = "Log2 Intensidad",
         ylab = "Frecuencia",
         breaks = 100,
         col = "lightblue",
         border = "black")
    abline(v = median(intensities[, i], na.rm = TRUE), 
           col = "red", lwd = 2, lty = 2)
    legend("topright", 
           legend = paste("Mediana =", round(median(intensities[, i], na.rm = TRUE), 2)),
           col = "red", lty = 2, lwd = 2)
}
dev.off()
cat("  ✓ Histograma guardado: 01_histogram_intensities.png\n")

# ------------------------------------------------------------------------------
# 6.2. Boxplots de Intensidades
# ------------------------------------------------------------------------------
cat("\n6.2. Generando boxplots de intensidades...\n")

png(file.path(qc_dir, "02_boxplot_intensities.png"), 
    width = 14, height = 8, units = "in", res = 300)

# Boxplot de intensidades por muestra
par(mar = c(8, 4, 4, 2))
boxplot(intensities, 
        main = "Distribución de Intensidades por Muestra (Boxplots)",
        xlab = "",
        ylab = "Log2 Intensidad",
        las = 2,
        col = brewer.pal(n = min(8, ncol(intensities)), name = "Set2"),
        cex.axis = 0.8)
mtext("Muestras", side = 1, line = 6)
abline(h = median(intensities, na.rm = TRUE), 
       col = "red", lwd = 2, lty = 2)
dev.off()
cat("  ✓ Boxplot guardado: 02_boxplot_intensities.png\n")

# ------------------------------------------------------------------------------
# 6.3. MA-plots (para detectar sesgos de intensidad)
# ------------------------------------------------------------------------------
cat("\n6.3. Generando MA-plots...\n")

# Calcular mediana de todas las muestras como referencia
median_intensity <- rowMedians(intensities, na.rm = TRUE)

png(file.path(qc_dir, "03_MAplots.png"), 
    width = 14, height = 10, units = "in", res = 300)

par(mfrow = c(2, 3))
for (i in 1:ncol(intensities)) {
    M <- intensities[, i] - median_intensity
    A <- (intensities[, i] + median_intensity) / 2
    
    # Muestrear para visualización (si hay muchos puntos)
    if (length(M) > 50000) {
        idx <- sample(1:length(M), 50000)
        M <- M[idx]
        A <- A[idx]
    }
    
    smoothScatter(A, M,
                  main = paste("MA-plot\n", sample_names[i]),
                  xlab = "A = (log2(Intensity) + log2(Mediana))/2",
                  ylab = "M = log2(Intensity) - log2(Mediana)",
                  pch = 19,
                  cex = 0.3)
    abline(h = 0, col = "red", lwd = 2, lty = 2)
    abline(h = c(-1, 1), col = "orange", lwd = 1, lty = 3)
    
    # Calcular y mostrar sesgo
    bias <- median(M, na.rm = TRUE)
    legend("topright", 
           legend = paste("Sesgo =", round(bias, 3)),
           col = ifelse(abs(bias) > 0.5, "red", "black"),
           text.col = ifelse(abs(bias) > 0.5, "red", "black"))
}
dev.off()
cat("  ✓ MA-plots guardados: 03_MAplots.png\n")

# ------------------------------------------------------------------------------
# 6.4. Análisis PCA (Principal Component Analysis)
# ------------------------------------------------------------------------------
cat("\n6.4. Realizando análisis PCA...\n")

# Calcular PCA sobre una muestra de probes (para eficiencia)
set.seed(123)
n_probes <- min(10000, nrow(intensities))
selected_probes <- sample(1:nrow(intensities), n_probes)
pca_data <- t(intensities[selected_probes, ])

# Realizar PCA
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
pca_summary <- summary(pca_result)

# Crear gráficos PCA
png(file.path(qc_dir, "04_PCA.png"), 
    width = 14, height = 10, units = "in", res = 300)

par(mfrow = c(2, 2))

# PCA Plot PC1 vs PC2
plot(pca_result$x[, 1], pca_result$x[, 2],
     main = "PCA: PC1 vs PC2",
     xlab = paste("PC1 (", round(pca_summary$importance[2, 1] * 100, 1), "%)"),
     ylab = paste("PC2 (", round(pca_summary$importance[2, 2] * 100, 1), "%)"),
     pch = 19,
     cex = 1.5,
     col = brewer.pal(n = min(8, ncol(intensities)), name = "Set2"))
text(pca_result$x[, 1], pca_result$x[, 2], 
     labels = sample_names, 
     pos = 3, cex = 0.7)
grid()

# PCA Plot PC1 vs PC3
plot(pca_result$x[, 1], pca_result$x[, 3],
     main = "PCA: PC1 vs PC3",
     xlab = paste("PC1 (", round(pca_summary$importance[2, 1] * 100, 1), "%)"),
     ylab = paste("PC3 (", round(pca_summary$importance[2, 3] * 100, 1), "%)"),
     pch = 19,
     cex = 1.5,
     col = brewer.pal(n = min(8, ncol(intensities)), name = "Set2"))
text(pca_result$x[, 1], pca_result$x[, 3], 
     labels = sample_names, 
     pos = 3, cex = 0.7)
grid()

# Scree plot
barplot(pca_summary$importance[2, 1:min(10, ncol(pca_result$x))],
        main = "Varianza Explicada por Componentes Principales",
        xlab = "Componente Principal",
        ylab = "Proporción de Varianza",
        col = "steelblue",
        names.arg = paste("PC", 1:min(10, ncol(pca_result$x)), sep = ""))

# Biplot (primeras dos componentes)
biplot(pca_result, 
       main = "Biplot PCA",
       cex = 0.7,
       xlabs = sample_names)

dev.off()
cat("  ✓ Gráficos PCA guardados: 04_PCA.png\n")

# Guardar resultados de PCA
save(pca_result, pca_summary, file = file.path(qc_dir, "PCA_results.RData"))

# ------------------------------------------------------------------------------
# 6.5. Clustering Jerárquico
# ------------------------------------------------------------------------------
cat("\n6.5. Realizando clustering jerárquico...\n")

# Calcular matriz de distancia (usando correlación de Pearson)
cor_matrix <- cor(intensities, use = "pairwise.complete.obs")
dist_matrix <- as.dist(1 - cor_matrix)

# Realizar clustering
hc <- hclust(dist_matrix, method = "ward.D2")

png(file.path(qc_dir, "05_clustering_hierarchical.png"), 
    width = 12, height = 8, units = "in", res = 300)

par(mar = c(8, 4, 4, 2))
plot(hc,
     main = "Clustering Jerárquico de Muestras\n(Distancia basada en correlación)",
     xlab = "",
     ylab = "Distancia",
     sub = "",
     labels = sample_names,
     cex = 0.8)
rect.hclust(hc, k = 2, border = "red")
dev.off()
cat("  ✓ Dendrograma guardado: 05_clustering_hierarchical.png\n")

# Heatmap de correlación
png(file.path(qc_dir, "06_heatmap_correlation.png"), 
    width = 10, height = 10, units = "in", res = 300)

heatmap.2(cor_matrix,
          main = "Matriz de Correlación entre Muestras",
          trace = "none",
          col = colorRampPalette(c("blue", "white", "red"))(100),
          dendrogram = "both",
          Rowv = as.dendrogram(hc),
          Colv = as.dendrogram(hc),
          key = TRUE,
          keysize = 1.5,
          density.info = "none",
          margins = c(10, 10),
          cexRow = 0.8,
          cexCol = 0.8)
dev.off()
cat("  ✓ Heatmap de correlación guardado: 06_heatmap_correlation.png\n")


# ------------------------------------------------------------------------------
# 6.6. Resumen de Métricas de Calidad
# ------------------------------------------------------------------------------
cat("\n6.7. Generando resumen de métricas de calidad...\n")

# Calcular métricas básicas
qc_metrics <- data.frame(
    Sample = sample_names,
    GSM_ID = gsm_ids,
    Median_Intensity = apply(intensities, 2, median, na.rm = TRUE),
    Mean_Intensity = apply(intensities, 2, mean, na.rm = TRUE),
    SD_Intensity = apply(intensities, 2, sd, na.rm = TRUE),
    IQR_Intensity = apply(intensities, 2, IQR, na.rm = TRUE),
    stringsAsFactors = FALSE
)

# Agregar información del tratamiento
qc_metrics$Treatment <- ifelse(grepl("Veh|VEH", qc_metrics$Sample), "Vehicle", "CBD")

# Guardar métricas
write.csv(qc_metrics, 
          file = file.path(qc_dir, "QC_metrics_summary.csv"), 
          row.names = FALSE)
cat("  ✓ Resumen de métricas guardado: QC_metrics_summary.csv\n")

# Mostrar resumen
cat("\nResumen de Métricas de Calidad:\n")
print(qc_metrics)

# Identificar posibles outliers
cat("\nAnálisis de Outliers:\n")
# Outliers basados en desviación de la mediana
median_deviation <- abs(qc_metrics$Median_Intensity - median(qc_metrics$Median_Intensity))
outlier_threshold <- 2 * sd(qc_metrics$Median_Intensity)
potential_outliers <- qc_metrics[median_deviation > outlier_threshold, ]

if (nrow(potential_outliers) > 0) {
    cat("⚠ Muestras potencialmente problemáticas (desviación > 2 SD):\n")
    print(potential_outliers[, c("Sample", "GSM_ID", "Median_Intensity")])
} else {
    cat("✓ No se detectaron outliers obvios basados en intensidad mediana\n")
}

# ------------------------------------------------------------------------------
# 6.7. Análisis Automatizado con arrayQualityMetrics (opcional, puede ser lento)
# ------------------------------------------------------------------------------
cat("\n6.8. Ejecutando arrayQualityMetrics (análisis automatizado)...\n")
cat("  Nota: Este análisis puede tardar varios minutos...\n")

tryCatch({
    aqm_dir <- file.path(qc_dir, "arrayQualityMetrics")
    if (!dir.exists(aqm_dir)) {
        dir.create(aqm_dir, recursive = TRUE)
    }
    
    # arrayQualityMetrics genera un reporte HTML completo
    arrayQualityMetrics(raw_data, 
                       outdir = aqm_dir,
                       force = TRUE,
                       do.logtransform = TRUE)
    
    cat("  ✓ Reporte HTML de arrayQualityMetrics generado en:", aqm_dir, "\n")
    cat("  Abre el archivo index.html en tu navegador para ver el reporte completo\n")
    
}, error = function(e) {
    cat("  ⚠ No se pudo ejecutar arrayQualityMetrics:", e$message, "\n")
    cat("  Esto puede deberse a problemas de memoria o dependencias\n")
})

cat("\n", rep("=", 70), "\n")
cat("ANÁLISIS DE CALIDAD COMPLETADO\n")
cat(rep("=", 70), "\n")
cat("Todos los gráficos y métricas se han guardado en:", qc_dir, "\n\n")

# ==============================================================================
# 7. Preprocesamiento y Normalización (RMA)
# ==============================================================================
cat("\n", rep("=", 70), "\n")
cat("INICIANDO PREPROCESAMIENTO Y NORMALIZACIÓN\n")
cat(rep("=", 70), "\n\n")

# ------------------------------------------------------------------------------
# 7.1. Normalización RMA
# ------------------------------------------------------------------------------
cat("7.1. Aplicando normalización RMA...\n")
cat("  Nota: RMA incluye background correction, quantile normalization y summarization\n")
cat("  IMPORTANTE: La normalización se aplica a TODOS los probes del array\n")
cat("  (no solo a los anotados). El filtrado por anotación ocurre después.\n")
cat("  Esto puede tardar varios minutos...\n")

# Para Affymetrix Gene 1.0 ST, oligo::rma() usa automáticamente la definición correcta del chip
# RMA normaliza TODOS los probes del array, no solo los anotados
eset <- rma(raw_data)

cat("  ✓ Normalización RMA completada\n")
cat("\nInformación del ExpressionSet normalizado:\n")
print(eset)
cat("\nDimensiones:", dim(eset), "\n")
cat("Número de features:", nrow(eset), "\n")
cat("Número de muestras:", ncol(eset), "\n")

# Extraer matriz de expresión normalizada
exprs_mat <- exprs(eset)
cat("\nRango de valores de expresión (log2):\n")
cat("  Mínimo:", min(exprs_mat, na.rm = TRUE), "\n")
cat("  Máximo:", max(exprs_mat, na.rm = TRUE), "\n")
cat("  Mediana:", median(exprs_mat, na.rm = TRUE), "\n")

# Copiar phenodata al ExpressionSet normalizado
pData(eset) <- pData(raw_data)
cat("\n✓ Phenodata asociado al ExpressionSet normalizado\n")

# ------------------------------------------------------------------------------
# 7.2. Anotación de Probes a Genes
# ------------------------------------------------------------------------------
cat("\n7.2. Mapeando probes a genes usando hugene10sttranscriptcluster.db...\n")

# Obtener IDs de probes
probe_ids <- rownames(exprs_mat)
cat("  Número de probes:", length(probe_ids), "\n")

# Mapear probes a símbolos de genes
# hugene10sttranscriptcluster.db contiene las anotaciones para Affymetrix Human Gene 1.0 ST
probe_to_symbol <- mapIds(hugene10sttranscriptcluster.db,
                          keys = probe_ids,
                          column = "SYMBOL",
                          keytype = "PROBEID",
                          multiVals = "first")

probe_to_entrez <- mapIds(hugene10sttranscriptcluster.db,
                          keys = probe_ids,
                          column = "ENTREZID",
                          keytype = "PROBEID",
                          multiVals = "first")

probe_to_genename <- mapIds(hugene10sttranscriptcluster.db,
                           keys = probe_ids,
                           column = "GENENAME",
                           keytype = "PROBEID",
                           multiVals = "first")

# Crear data frame de anotaciones
annotations_df <- data.frame(
    PROBEID = probe_ids,
    SYMBOL = probe_to_symbol,
    ENTREZID = probe_to_entrez,
    GENENAME = probe_to_genename,
    stringsAsFactors = FALSE
)

# Estadísticas de anotación
n_annotated <- sum(!is.na(annotations_df$SYMBOL))
cat("  Probes anotados con SYMBOL:", n_annotated, "(", 
    round(n_annotated/length(probe_ids)*100, 1), "%)\n")
cat("  Probes sin anotación:", sum(is.na(annotations_df$SYMBOL)), "\n")

# Agregar anotaciones al fData del ExpressionSet
fData(eset) <- annotations_df
cat("  ✓ Anotaciones agregadas al ExpressionSet\n")

# ------------------------------------------------------------------------------
# 7.3. Summarización a nivel de gen
# ------------------------------------------------------------------------------
cat("\n7.3. Creando matriz de expresión a nivel de gen...\n")
cat("  Nota: Solo se usan probes con anotación para la summarización a nivel de gen\n")
cat("  (la normalización ya se aplicó a todos los probes)\n")

# Filtrar probes sin anotación (solo para summarización a nivel de gen)
# La normalización ya se aplicó a todos los probes
annotated_probes <- !is.na(annotations_df$SYMBOL)
exprs_annotated <- exprs_mat[annotated_probes, ]
annotations_annotated <- annotations_df[annotated_probes, ]

cat("  Probes totales normalizados:", nrow(exprs_mat), "\n")
cat("  Probes con anotación (usados para summarización):", nrow(exprs_annotated), "\n")
cat("  Probes sin anotación (excluidos de summarización):", 
    nrow(exprs_mat) - nrow(exprs_annotated), "\n")

# Summarizar a nivel de gen (promediar múltiples probes por gen)
# Usar limma::avereps para promediar probes que mapean al mismo gen
exprs_gene <- avereps(exprs_annotated, 
                      ID = annotations_annotated$SYMBOL)

cat("  Genes únicos después de summarización:", nrow(exprs_gene), "\n")
cat("  Dimensiones de matriz a nivel gen:", dim(exprs_gene), "\n")

# Crear data frame con los datos normalizados y anotados
normalized_counts <- as.data.frame(exprs_gene)
normalized_counts$SYMBOL <- rownames(normalized_counts)

# Reordenar columnas para que SYMBOL esté primero
normalized_counts <- normalized_counts[, c("SYMBOL", 
                                           setdiff(colnames(normalized_counts), "SYMBOL"))]

# Agregar información adicional de anotación (tomar la primera ocurrencia de cada gen)
gene_annotations <- annotations_annotated[!duplicated(annotations_annotated$SYMBOL), 
                                          c("SYMBOL", "ENTREZID", "GENENAME")]
gene_annotations <- gene_annotations[gene_annotations$SYMBOL %in% rownames(exprs_gene), ]

# Merge con normalized_counts
normalized_counts <- merge(gene_annotations, normalized_counts, 
                          by = "SYMBOL", all.y = TRUE)

# Reordenar para que las columnas de expresión estén después de las anotaciones
sample_cols <- setdiff(colnames(normalized_counts), 
                      c("SYMBOL", "ENTREZID", "GENENAME"))
normalized_counts <- normalized_counts[, c("SYMBOL", "ENTREZID", "GENENAME", sample_cols)]

cat("  ✓ Matriz de expresión a nivel de gen creada\n")

# ------------------------------------------------------------------------------
# 7.4. Filtrado de Baja Señal y Baja Varianza
# ------------------------------------------------------------------------------
cat("\n7.4. Aplicando filtrado de baja señal y baja varianza...\n")

# Extraer matriz de expresión (sin columnas de anotación)
exprs_for_filtering <- as.matrix(normalized_counts[, sample_cols])
rownames(exprs_for_filtering) <- normalized_counts$SYMBOL

# Calcular estadísticas para filtrado
cat("  Calculando estadísticas de expresión y varianza...\n")

# 1. Filtrado por baja expresión (background threshold)
# Estimar background como el percentil 5 de todas las intensidades
background_threshold <- quantile(exprs_for_filtering, probs = 0.05, na.rm = TRUE)
cat("  Background threshold (percentil 5):", round(background_threshold, 3), "\n")

# Número mínimo de muestras que deben tener expresión > background
# Por defecto: al menos 50% de las muestras
min_samples_above_bg <- ceiling(ncol(exprs_for_filtering) * 0.5)
cat("  Mínimo de muestras con expresión > background:", min_samples_above_bg, 
    "de", ncol(exprs_for_filtering), "\n")

# Contar cuántas muestras tienen expresión > background para cada gen
samples_above_bg <- rowSums(exprs_for_filtering > background_threshold, na.rm = TRUE)
genes_above_bg <- samples_above_bg >= min_samples_above_bg

cat("  Genes con expresión > background en al menos", min_samples_above_bg, 
    "muestras:", sum(genes_above_bg), "(", 
    round(sum(genes_above_bg)/nrow(exprs_for_filtering)*100, 1), "%)\n")

# 2. Filtrado por baja varianza
# Calcular varianza de cada gen
gene_variances <- apply(exprs_for_filtering, 1, var, na.rm = TRUE)

# Threshold de varianza: quantil 0.2 (eliminar genes con varianza en el 20% más bajo)
variance_threshold <- quantile(gene_variances, probs = 0.2, na.rm = TRUE)
cat("  Threshold de varianza (quantil 0.2):", round(variance_threshold, 3), "\n")

genes_above_variance <- gene_variances > variance_threshold
cat("  Genes con varianza > threshold:", sum(genes_above_variance), "(", 
    round(sum(genes_above_variance)/nrow(exprs_for_filtering)*100, 1), "%)\n")

# 3. Combinar filtros (genes que pasan AMBOS filtros)
genes_passing_filters <- genes_above_bg & genes_above_variance
cat("\n  Genes que pasan ambos filtros:", sum(genes_passing_filters), "(", 
    round(sum(genes_passing_filters)/nrow(exprs_for_filtering)*100, 1), "%)\n")
cat("  Genes eliminados por filtros:", sum(!genes_passing_filters), "(", 
    round(sum(!genes_passing_filters)/nrow(exprs_for_filtering)*100, 1), "%)\n")

# Aplicar filtros
exprs_filtered <- exprs_for_filtering[genes_passing_filters, ]
normalized_counts_filtered <- normalized_counts[genes_passing_filters, ]

# Actualizar exprs_gene también
exprs_gene_filtered <- exprs_gene[genes_passing_filters, ]

cat("\n  Dimensiones antes del filtrado:", dim(exprs_for_filtering), "\n")
cat("  Dimensiones después del filtrado:", dim(exprs_filtered), "\n")

# Crear resumen del filtrado
filtering_summary <- data.frame(
    Criterio = c("Total genes", 
                 "Expresión > background (≥50% muestras)",
                 "Varianza > quantil 0.2",
                 "Pasan ambos filtros",
                 "Eliminados"),
    Numero = c(nrow(exprs_for_filtering),
               sum(genes_above_bg),
               sum(genes_above_variance),
               sum(genes_passing_filters),
               sum(!genes_passing_filters)),
    Porcentaje = c(100,
                   round(sum(genes_above_bg)/nrow(exprs_for_filtering)*100, 1),
                   round(sum(genes_above_variance)/nrow(exprs_for_filtering)*100, 1),
                   round(sum(genes_passing_filters)/nrow(exprs_for_filtering)*100, 1),
                   round(sum(!genes_passing_filters)/nrow(exprs_for_filtering)*100, 1))
)

cat("\n  Resumen del filtrado:\n")
print(filtering_summary)

# Guardar resumen de filtrado
write.csv(filtering_summary,
          file = file.path(output_dir, "GSE57978_filtering_summary.csv"),
          row.names = FALSE)

cat("  ✓ Resumen de filtrado guardado: GSE57978_filtering_summary.csv\n")

# Actualizar objetos con datos filtrados
normalized_counts <- normalized_counts_filtered
exprs_gene <- exprs_gene_filtered

cat("  ✓ Filtrado completado\n")

# ------------------------------------------------------------------------------
# 7.5. Guardar datos normalizados (filtrados)
# ------------------------------------------------------------------------------
cat("\n7.5. Guardando datos normalizados y filtrados...\n")

# Guardar ExpressionSet completo
save(eset, file = file.path(output_dir, "GSE57978_normalized_eset.RData"))
cat("  ✓ ExpressionSet normalizado guardado: GSE57978_normalized_eset.RData\n")

# Guardar matriz de expresión normalizada (probe level)
write.csv(exprs_mat, 
          file = file.path(output_dir, "GSE57978_normalized_probe_level.csv"),
          row.names = TRUE)
cat("  ✓ Matriz normalizada a nivel de probe guardada: GSE57978_normalized_probe_level.csv\n")

# Guardar matriz de expresión normalizada (gene level) como tabla (FILTRADA)
write.csv(normalized_counts,
          file = file.path(output_dir, "GSE57978_normalized_gene_level_filtered.csv"),
          row.names = FALSE)
cat("  ✓ Matriz normalizada y filtrada a nivel de gen guardada: GSE57978_normalized_gene_level_filtered.csv\n")

# También guardar versión sin filtrar para referencia
# Crear data frame sin filtrar antes de aplicar filtros
exprs_unfiltered_df <- as.data.frame(exprs_for_filtering)
exprs_unfiltered_df$SYMBOL <- rownames(exprs_unfiltered_df)
normalized_counts_unfiltered <- merge(gene_annotations, 
                                      exprs_unfiltered_df,
                                      by = "SYMBOL", all.y = TRUE)
normalized_counts_unfiltered <- normalized_counts_unfiltered[, 
    c("SYMBOL", "ENTREZID", "GENENAME", sample_cols)]
write.csv(normalized_counts_unfiltered,
          file = file.path(output_dir, "GSE57978_normalized_gene_level_unfiltered.csv"),
          row.names = FALSE)
cat("  ✓ Matriz normalizada sin filtrar guardada: GSE57978_normalized_gene_level_unfiltered.csv\n")

# Guardar anotaciones
write.csv(annotations_df,
          file = file.path(output_dir, "GSE57978_probe_annotations.csv"),
          row.names = FALSE)
cat("  ✓ Anotaciones de probes guardadas: GSE57978_probe_annotations.csv\n")

cat("\n", rep("=", 70), "\n")
cat("PREPROCESAMIENTO Y NORMALIZACIÓN COMPLETADO\n")
cat(rep("=", 70), "\n\n")

# ==============================================================================
# 8. Análisis de Expresión Diferencial (DEGs)
# ==============================================================================
cat("\n", rep("=", 70), "\n")
cat("INICIANDO ANÁLISIS DE EXPRESIÓN DIFERENCIAL\n")
cat(rep("=", 70), "\n\n")

# Crear directorio para resultados de DEGs
degs_dir <- file.path(output_dir, "DEGs")
if (!dir.exists(degs_dir)) {
    dir.create(degs_dir, recursive = TRUE)
}
cat("Directorio de DEGs creado:", degs_dir, "\n\n")

# ------------------------------------------------------------------------------
# 8.1. Preparar diseño experimental y DIAGNÓSTICO
# ------------------------------------------------------------------------------
cat("8.1. Preparando diseño experimental...\n")
cat("  ⚠ DIAGNÓSTICO: Verificando estructura del experimento...\n\n")

# Crear variable de grupo basada en el tratamiento (CBD vs Vehicle)
treatment_raw <- ifelse(grepl("Veh|VEH", sample_names), "Vehicle", "CBD")
# Crear factor - Vehicle y CBD ya son nombres válidos en R, no necesitan make.names
treatment <- factor(treatment_raw, levels = c("Vehicle", "CBD"))
# No aplicar make.names aquí porque Vehicle y CBD ya son válidos
# levels(treatment) <- make.names(levels(treatment))  # COMENTADO: no necesario

# Extraer línea celular de los nombres de muestra
# Los nombres tienen patrones como: "0609", "CC4121", "4121T3"
cell_line_raw <- rep(NA, length(sample_names))
cell_line_raw[grepl("0609", sample_names, ignore.case = TRUE)] <- "Line_0609"
cell_line_raw[grepl("CC4121|4121", sample_names, ignore.case = TRUE) & 
               !grepl("T3", sample_names, ignore.case = TRUE)] <- "Line_CC4121"
cell_line_raw[grepl("4121T3|T3", sample_names, ignore.case = TRUE)] <- "Line_4121T3"

# Convertir a factor con nombres válidos en R
cell_line <- factor(cell_line_raw)
levels(cell_line) <- make.names(levels(cell_line))

cat("  📊 DIAGNÓSTICO 1: Distribución de muestras\n")
cat("  ===========================================\n")
cat("  Tabla de tratamiento:\n")
print(table(treatment))
cat("\n  Tabla de línea celular:\n")
print(table(cell_line))
cat("\n  Tabla cruzada (Línea x Tratamiento):\n")
print(table(cell_line, treatment))
cat("\n")

# Verificar que hay al menos una réplica por combinación
cat("  📊 DIAGNÓSTICO 2: Réplicas por combinación\n")
cat("  ===========================================\n")
replicates_table <- table(cell_line, treatment)
cat("  Réplicas por combinación:\n")
print(replicates_table)
min_replicates <- min(replicates_table)
cat("  Mínimo de réplicas:", min_replicates, "\n")
if (min_replicates < 2) {
    cat("  ⚠ ADVERTENCIA: Algunas combinaciones tienen menos de 2 réplicas\n")
}
cat("\n")

# Mostrar asignación de muestras
cat("  📊 DIAGNÓSTICO 3: Asignación de muestras\n")
cat("  ===========================================\n")
sample_info <- data.frame(
    Sample = sample_names,
    GSM_ID = gsm_ids,
    Cell_Line = cell_line,
    Treatment = treatment,
    stringsAsFactors = FALSE
)
print(sample_info)
cat("\n")

# ==============================================================================
# ESTRATEGIA DE ANÁLISIS: Diseño Pareado con duplicateCorrelation
# ==============================================================================
# Dado que cada combinación tiene solo 1 réplica, usamos un diseño pareado
# donde cada línea celular tiene su propio control y tratamiento (paired design)
# Esto requiere usar duplicateCorrelation() para modelar la correlación entre
# muestras emparejadas dentro de cada bloque (línea celular)

cat("  📊 ESTRATEGIA DE ANÁLISIS\n")
cat("  ===========================================\n")
cat("  Diseño: Pareado (Paired Design)\n")
cat("  Cada línea celular tiene 1 control (Vehicle) y 1 tratamiento (CBD)\n")
cat("  Usaremos duplicateCorrelation() para modelar correlación entre pares\n")
cat("  Esto es apropiado cuando hay estructura de bloqueo pero pocas réplicas\n\n")

# Crear variable de bloque (cada línea celular es un bloque)
block <- cell_line
cat("  Bloques (líneas celulares):\n")
print(table(block))
cat("\n")

# Crear matriz de diseño simple: solo tratamiento
# En un diseño pareado, no necesitamos incluir el bloque en el diseño
# porque duplicateCorrelation lo maneja
cat("  📊 DIAGNÓSTICO 4: Matriz de diseño\n")
cat("  ===========================================\n")
cat("  Diseño: ~ treatment (diseño simple para paired analysis)\n")
cat("  El bloqueo se manejará con duplicateCorrelation()\n\n")

design <- model.matrix(~ treatment)
# Guardar nombres originales antes de aplicar make.names
design_cols_original <- colnames(design)

cat("  Nombres de columnas ANTES de make.names:\n")
print(design_cols_original)

# Aplicar make.names solo si es necesario
colnames(design) <- make.names(colnames(design))

cat("  Columnas de la matriz de diseño (originales):\n")
print(design_cols_original)
cat("\n  Columnas de la matriz de diseño (después de make.names):\n")
print(colnames(design))
cat("\n  Dimensiones:", dim(design), "\n")
cat("  Rango (debe ser igual al número de columnas):", qr(design)$rank, "\n")
if (qr(design)$rank < ncol(design)) {
    cat("  ⚠ ERROR: Matriz de diseño es singular (columnas colineales)\n")
} else {
    cat("  ✓ Matriz de diseño es de rango completo\n")
}
cat("\n  Primeras filas de la matriz de diseño:\n")
print(head(design))
cat("\n")

# En un diseño ~ treatment, la primera columna es el intercepto (Vehicle)
# y la segunda columna es el efecto de CBD vs Vehicle
# No necesitamos crear contrastes explícitos, el coeficiente de treatmentCBD
# ya representa la diferencia CBD vs Vehicle

cat("  ✓ Diseño experimental creado (diseño pareado)\n")
cat("  El coeficiente de 'treatmentCBD' representa CBD vs Vehicle\n\n")

# ------------------------------------------------------------------------------
# 8.2. Análisis de expresión diferencial con limma (diseño pareado)
# ------------------------------------------------------------------------------
cat("8.2. Realizando análisis de expresión diferencial con limma...\n")
cat("  Usando duplicateCorrelation() para diseño pareado\n\n")

# PASO 1: Estimar correlación entre muestras emparejadas (dentro de cada bloque)
cat("  Paso 1: Estimando correlación entre muestras pareadas...\n")
corfit <- duplicateCorrelation(exprs_gene, design, block = block)
cat("  Correlación estimada entre pares:", round(corfit$consensus.correlation, 3), "\n")
cat("  (Valores cercanos a 1 indican alta correlación entre muestras pareadas)\n\n")

# PASO 2: Ajustar el modelo lineal usando la correlación estimada
cat("  Paso 2: Ajustando modelo lineal con correlación de bloques...\n")
fit <- lmFit(exprs_gene, design, block = block, correlation = corfit$consensus.correlation)
cat("  ✓ Modelo ajustado\n\n")

# PASO 3: Aplicar eBayes para suavizar varianzas
cat("  Paso 3: Aplicando eBayes para suavizar varianzas...\n")
fit2 <- eBayes(fit, trend = TRUE, robust = TRUE)
cat("  ✓ eBayes aplicado\n\n")

# En un diseño ~ treatment, el coeficiente de "treatmentCBD" ya es la comparación
# No necesitamos contrasts.fit, solo extraer el coeficiente directamente
coef_name <- colnames(design)[grepl("CBD|treatment", colnames(design), ignore.case = TRUE)]
if (length(coef_name) > 1) {
    coef_name <- coef_name[grepl("CBD", coef_name, ignore.case = TRUE)][1]
}
if (length(coef_name) == 0) {
    # Si no encontramos, usar la segunda columna (después del intercepto)
    coef_name <- colnames(design)[2]
}
cat("  Coeficiente usado para la comparación:", coef_name, "\n\n")

# Extraer resultados usando el coeficiente directamente
# En diseño ~ treatment, el coeficiente de treatmentCBD es la comparación
results <- topTable(fit2, coef = coef_name, number = Inf, sort.by = "P")

# DIAGNÓSTICO: Verificar resultados sin filtrar
cat("  📊 DIAGNÓSTICO 6: Resultados del modelo (SIN FILTRAR)\n")
cat("  ===========================================\n")
cat("  Top 20 genes por p-value (sin ajustar):\n")
print(head(results[order(results$P.Value), ], 20))
cat("\n  Estadísticas de logFC:\n")
cat("    Mediana:", median(results$logFC, na.rm = TRUE), "\n")
cat("    Media:", mean(results$logFC, na.rm = TRUE), "\n")
cat("    SD:", sd(results$logFC, na.rm = TRUE), "\n")
cat("    Rango:", range(results$logFC, na.rm = TRUE), "\n")
cat("\n  Estadísticas de p-values:\n")
cat("    Genes con P.Value < 0.05:", sum(results$P.Value < 0.05, na.rm = TRUE), "\n")
cat("    Genes con P.Value < 0.01:", sum(results$P.Value < 0.01, na.rm = TRUE), "\n")
cat("    Genes con adj.P.Val < 0.05:", sum(results$adj.P.Val < 0.05, na.rm = TRUE), "\n")
cat("    Genes con adj.P.Val < 0.1:", sum(results$adj.P.Val < 0.1, na.rm = TRUE), "\n")
cat("\n  Estadísticas de logFC absoluto:\n")
cat("    Genes con |logFC| >= 0.5:", sum(abs(results$logFC) >= 0.5, na.rm = TRUE), "\n")
cat("    Genes con |logFC| >= 1.0:", sum(abs(results$logFC) >= 1.0, na.rm = TRUE), "\n")
cat("    Genes con |logFC| >= 1.5:", sum(abs(results$logFC) >= 1.5, na.rm = TRUE), "\n")
cat("\n")

# Agregar anotaciones a los resultados
results$SYMBOL <- rownames(results)
results_annotated <- merge(
    normalized_counts[, c("SYMBOL", "ENTREZID", "GENENAME")],
    results,
    by = "SYMBOL",
    all.y = TRUE
)

# Reordenar columnas
results_annotated <- results_annotated[, c("SYMBOL", "ENTREZID", "GENENAME", 
                                           "logFC", "AveExpr", "t", "P.Value", 
                                           "adj.P.Val", "B")]

cat("  ✓ Análisis completado\n")
cat("  Total de genes analizados:", nrow(results_annotated), "\n\n")

# Definir umbrales para DEGs
fc_threshold <- 0.5     # cambio moderado (~1.41x)
pval_threshold <- 0.1   # FDR moderado  # P-value ajustado (FDR)



cat("  Umbrales para identificar DEGs:\n")
cat("    |logFC| >= ", fc_threshold, " (fold change >= 2x)\n")
cat("    adj.P.Val < ", pval_threshold, " (FDR < 5%)\n\n")

# Identificar DEGs
degs <- results_annotated[
    abs(results_annotated$logFC) >= fc_threshold & 
    results_annotated$adj.P.Val < pval_threshold,
]

# Separar up y down regulated
up_genes <- degs[degs$logFC > 0, ]
down_genes <- degs[degs$logFC < 0, ]

cat("  📊 RESULTADOS FINALES:\n")
cat("  ===========================================\n")
cat("    Total DEGs (|logFC| >=", fc_threshold, "y adj.P.Val <", pval_threshold, "):", 
    nrow(degs), "\n")
cat("    Genes up-regulated (CBD > Vehicle):", nrow(up_genes), "\n")
cat("    Genes down-regulated (CBD < Vehicle):", nrow(down_genes), "\n")
cat("\n")

# Si aún no hay DEGs, mostrar los top genes sin filtrar
if (nrow(degs) == 0) {
    cat("  ⚠ ADVERTENCIA: No se encontraron DEGs con los umbrales seleccionados\n")
    cat("  Mostrando top 10 genes por significancia (sin filtro de FC):\n")
    top_by_pval <- results_annotated[order(results_annotated$adj.P.Val), ][1:min(10, nrow(results_annotated)), ]
    print(top_by_pval[, c("SYMBOL", "logFC", "P.Value", "adj.P.Val")])
    cat("\n  Sugerencia: Considera usar umbrales más permisivos o revisar el diseño experimental\n")
}

# ------------------------------------------------------------------------------
# 8.3. Volcano Plot
# ------------------------------------------------------------------------------
cat("\n8.3. Generando Volcano plot...\n")

png(file.path(degs_dir, "volcano_plot.png"), 
    width = 12, height = 10, units = "in", res = 300)

# Preparar datos para el plot
plot_data <- data.frame(
    logFC = results_annotated$logFC,
    neg_log10_pval = -log10(results_annotated$adj.P.Val),
    SYMBOL = results_annotated$SYMBOL,
    is_DEG = abs(results_annotated$logFC) >= fc_threshold & 
             results_annotated$adj.P.Val < pval_threshold,
    is_up = results_annotated$logFC > 0
)

# Crear volcano plot con ggplot2
volcano_plot <- ggplot(plot_data, aes(x = logFC, y = neg_log10_pval)) +
    geom_point(aes(color = is_DEG, alpha = is_DEG), size = 1.5) +
    scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "red"),
                      labels = c("FALSE" = "No DEG", "TRUE" = "DEG"),
                      name = "") +
    scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 0.8), guide = "none") +
    geom_hline(yintercept = -log10(pval_threshold), 
               linetype = "dashed", color = "blue", linewidth = 1) +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), 
               linetype = "dashed", color = "blue", linewidth = 1) +
    labs(
        title = "Volcano Plot: CBD vs Vehicle",
        subtitle = paste("DEGs:", nrow(degs), "| Up:", nrow(up_genes), 
                        "| Down:", nrow(down_genes)),
        x = "Log2 Fold Change (CBD vs Vehicle)",
        y = "-Log10 Adjusted P-value"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12),
        legend.position = "right"
    ) +
    xlim(c(-max(abs(plot_data$logFC), na.rm = TRUE) * 1.1, 
           max(abs(plot_data$logFC), na.rm = TRUE) * 1.1))

print(volcano_plot)
dev.off()
cat("  ✓ Volcano plot guardado: volcano_plot.png\n")

# ------------------------------------------------------------------------------
# 8.4. MA-plot
# ------------------------------------------------------------------------------
cat("\n8.4. Generando MA-plot...\n")

png(file.path(degs_dir, "MA_plot.png"), 
    width = 12, height = 10, units = "in", res = 300)

ma_plot <- ggplot(results_annotated, aes(x = AveExpr, y = logFC)) +
    geom_point(aes(color = adj.P.Val < pval_threshold & abs(logFC) >= fc_threshold,
                   alpha = adj.P.Val < pval_threshold & abs(logFC) >= fc_threshold),
               size = 1.5) +
    scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "red"),
                      labels = c("FALSE" = "No DEG", "TRUE" = "DEG"),
                      name = "") +
    scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 0.8), guide = "none") +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    geom_hline(yintercept = c(-fc_threshold, fc_threshold), 
               linetype = "dashed", color = "blue", linewidth = 1) +
    labs(
        title = "MA-plot: CBD vs Vehicle",
        subtitle = paste("A = Average Expression, M = Log2 Fold Change"),
        x = "Average Expression (A)",
        y = "Log2 Fold Change (M)"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12),
        legend.position = "right"
    )

print(ma_plot)
dev.off()
cat("  ✓ MA-plot guardado: MA_plot.png\n")

# ------------------------------------------------------------------------------
# 8.5. Heatmap de Top Genes
# ------------------------------------------------------------------------------
cat("\n8.5. Generando heatmap de top genes...\n")

# Seleccionar top genes (por p-value ajustado)
n_top_genes <- min(200, nrow(degs))
if (nrow(degs) > 0) {
    top_genes <- degs[order(degs$adj.P.Val), ][1:min(n_top_genes, nrow(degs)), ]
    top_gene_symbols <- top_genes$SYMBOL
    
    # Extraer expresión de top genes
    top_exprs <- exprs_gene[top_gene_symbols, , drop = FALSE]
    
    # Normalizar por filas (z-score) para mejor visualización
    top_exprs_scaled <- t(scale(t(top_exprs)))
    
    # Crear anotación de colores para grupos
    col_annotation <- data.frame(
        Treatment = treatment,
        row.names = colnames(top_exprs)
    )
    
    # Colores para grupos
    treatment_colors <- list(
        Treatment = c("Vehicle" = "lightblue", "CBD" = "orange")
    )
    
    # Crear heatmap
    png(file.path(degs_dir, "heatmap_top_genes.png"), 
        width = 14, height = 12, units = "in", res = 300)
    
    # Usar pheatmap si está disponible, sino usar heatmap.2
    if (require("pheatmap", quietly = TRUE)) {
        library(pheatmap)
        pheatmap(
            top_exprs_scaled,
            annotation_col = col_annotation,
            annotation_colors = treatment_colors,
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            show_rownames = FALSE,
            show_colnames = TRUE,
            main = paste("Heatmap: Top", nrow(top_exprs), "DEGs (CBD vs Vehicle)"),
            color = colorRampPalette(c("blue", "white", "red"))(100),
            fontsize = 10,
            fontsize_col = 8
        )
    } else {
        # Usar heatmap.2 como alternativa
        heatmap.2(
            top_exprs_scaled,
            main = paste("Heatmap: Top", nrow(top_exprs), "DEGs"),
            trace = "none",
            col = colorRampPalette(c("blue", "white", "red"))(100),
            dendrogram = "both",
            Rowv = TRUE,
            Colv = TRUE,
            scale = "none",
            key = TRUE,
            keysize = 1.5,
            density.info = "none",
            margins = c(10, 8),
            cexRow = 0.6,
            cexCol = 0.8,
            labRow = FALSE
        )
    }
    
    dev.off()
    cat("  ✓ Heatmap guardado: heatmap_top_genes.png (", nrow(top_exprs), "genes)\n")
} else {
    cat("  ⚠ No hay DEGs para generar heatmap\n")
}

# ------------------------------------------------------------------------------
# 8.6. Guardar tablas de DEGs
# ------------------------------------------------------------------------------
cat("\n8.6. Guardando tablas de DEGs...\n")

# Guardar todos los resultados (con anotaciones)
write.table(results_annotated,
            file = file.path(degs_dir, "all_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  ✓ Todos los resultados guardados: all_results.tsv\n")

# Guardar todos los DEGs
write.table(degs,
            file = file.path(degs_dir, "all_DEGs.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  ✓ Todos los DEGs guardados: all_DEGs.tsv\n")

# Guardar genes up-regulated (solo símbolos)
write.table(up_genes$SYMBOL,
            file = file.path(degs_dir, "up_genes.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("  ✓ Genes up-regulated guardados: up_genes.txt (", nrow(up_genes), "genes)\n")

# Guardar genes down-regulated (solo símbolos)
write.table(down_genes$SYMBOL,
            file = file.path(degs_dir, "down_genes.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("  ✓ Genes down-regulated guardados: down_genes.txt (", nrow(down_genes), "genes)\n")

# Guardar tablas completas de up y down con anotaciones
write.table(up_genes,
            file = file.path(degs_dir, "up_genes_complete.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  ✓ Tabla completa up-regulated guardada: up_genes_complete.tsv\n")

write.table(down_genes,
            file = file.path(degs_dir, "down_genes_complete.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("  ✓ Tabla completa down-regulated guardada: down_genes_complete.tsv\n")

# Resumen estadístico
summary_stats <- data.frame(
    Categoria = c("Total genes analizados",
                  "DEGs totales",
                  "Up-regulated (CBD > Vehicle)",
                  "Down-regulated (CBD < Vehicle)",
                  "No significativos"),
    Numero = c(nrow(results_annotated),
               nrow(degs),
               nrow(up_genes),
               nrow(down_genes),
               nrow(results_annotated) - nrow(degs)),
    Porcentaje = c(100,
                   round(nrow(degs)/nrow(results_annotated)*100, 2),
                   round(nrow(up_genes)/nrow(results_annotated)*100, 2),
                   round(nrow(down_genes)/nrow(results_annotated)*100, 2),
                   round((nrow(results_annotated) - nrow(degs))/nrow(results_annotated)*100, 2))
)

write.csv(summary_stats,
          file = file.path(degs_dir, "DEGs_summary.csv"),
          row.names = FALSE)
cat("  ✓ Resumen estadístico guardado: DEGs_summary.csv\n")

cat("\n  Resumen del análisis:\n")
print(summary_stats)

cat("\n", rep("=", 70), "\n")
cat("ANÁLISIS DE EXPRESIÓN DIFERENCIAL COMPLETADO\n")
cat(rep("=", 70), "\n\n")

# ==============================================================================
# 9. Guardar objetos para uso posterior
# ==============================================================================
cat("\nGuardando objetos...\n")
save(raw_data, file = file.path(output_dir, "GSE57978_raw_data.RData"))
cat("Objeto raw_data guardado en:", file.path(output_dir, "GSE57978_raw_data.RData"), "\n")

# Guardar phenodata
save(pheno_data, pheno_subset, file = file.path(output_dir, "GSE57978_phenodata.RData"))
cat("Phenodata guardado en:", file.path(output_dir, "GSE57978_phenodata.RData"), "\n")

# Guardar objetos de QC
save(intensities, qc_metrics, pca_result, pca_summary, cor_matrix, hc,
     file = file.path(output_dir, "GSE57978_QC_objects.RData"))
cat("Objetos de QC guardados en:", file.path(output_dir, "GSE57978_QC_objects.RData"), "\n")

# Guardar objetos de normalización (incluyendo datos filtrados)
save(exprs_mat, exprs_gene, normalized_counts, annotations_df,
     exprs_filtered, filtering_summary, background_threshold, variance_threshold,
     file = file.path(output_dir, "GSE57978_normalization_objects.RData"))
cat("Objetos de normalización guardados en:", 
    file.path(output_dir, "GSE57978_normalization_objects.RData"), "\n")

# Guardar objetos de análisis de expresión diferencial
save(fit, fit2, results_annotated, degs, up_genes, down_genes, 
     design, treatment, cell_line, block, corfit, coef_name, 
     sample_info, summary_stats,
     file = file.path(output_dir, "GSE57978_DEGs_objects.RData"))
cat("Objetos de DEGs guardados en:", 
    file.path(output_dir, "GSE57978_DEGs_objects.RData"), "\n")

# Guardar diagnóstico completo
diagnostic_info <- list(
    sample_info = sample_info,
    design_matrix = design,
    block_structure = table(block),
    correlation_estimate = corfit$consensus.correlation,
    coefficient_used = coef_name,
    final_thresholds = list(FC = fc_threshold, Pval = pval_threshold),
    model_stats = list(
        median_logFC = median(results$logFC, na.rm = TRUE),
        mean_logFC = mean(results$logFC, na.rm = TRUE),
        genes_Pval_005 = sum(results$P.Value < 0.05, na.rm = TRUE),
        genes_adjPval_005 = sum(results$adj.P.Val < 0.05, na.rm = TRUE),
        genes_adjPval_01 = sum(results$adj.P.Val < 0.1, na.rm = TRUE)
    )
)
save(diagnostic_info, file = file.path(degs_dir, "diagnostic_info.RData"))
write.csv(sample_info, file = file.path(degs_dir, "sample_info.csv"), row.names = FALSE)
cat("Información de diagnóstico guardada en:", degs_dir, "\n")

cat("\n", rep("=", 70), "\n")
cat("¡ANÁLISIS COMPLETADO!\n")
cat(rep("=", 70), "\n\n")
cat("RESUMEN DE OBJETOS DISPONIBLES:\n")
cat("  - raw_data: Datos crudos de expresión (ExpressionFeatureSet)\n")
cat("  - eset: Datos normalizados con RMA (ExpressionSet)\n")
cat("  - exprs_mat: Matriz de expresión normalizada a nivel de probe\n")
cat("  - exprs_gene: Matriz de expresión normalizada a nivel de gen\n")
cat("  - normalized_counts: Tabla con expresión normalizada y anotaciones (gene level)\n")
cat("  - annotations_df: Anotaciones de probes a genes\n")
cat("  - pheno_data: Phenodata completo\n")
cat("  - pheno_subset: Phenodata filtrado para las muestras cargadas\n\n")
cat("ARCHIVOS GENERADOS:\n")
cat("  - GSE57978_normalized_gene_level_filtered.csv: Tabla principal con expresión normalizada y filtrada\n")
cat("  - GSE57978_normalized_gene_level_unfiltered.csv: Expresión normalizada sin filtrar (referencia)\n")
cat("  - GSE57978_normalized_probe_level.csv: Expresión normalizada a nivel de probe\n")
cat("  - GSE57978_probe_annotations.csv: Anotaciones de probes\n")
cat("  - GSE57978_filtering_summary.csv: Resumen del filtrado aplicado\n")
cat("  - QC/: Directorio con todos los gráficos de calidad\n")
cat("  - DEGs/: Directorio con resultados de expresión diferencial:\n")
cat("    * volcano_plot.png: Volcano plot de DEGs\n")
cat("    * MA_plot.png: MA-plot de expresión diferencial\n")
cat("    * heatmap_top_genes.png: Heatmap de top genes diferencialmente expresados\n")
cat("    * all_DEGs.tsv: Tabla completa de todos los DEGs\n")
cat("    * up_genes.txt: Lista de genes up-regulated\n")
cat("    * down_genes.txt: Lista de genes down-regulated\n")
cat("    * all_results.tsv: Todos los resultados del análisis\n")
cat("    * DEGs_summary.csv: Resumen estadístico del análisis\n")
cat("  - GSE57978_analysis_summary.html: Reporte HTML completo con resumen de QC y DEGs\n\n")

# ==============================================================================
# 10. Generar Reporte HTML de Resumen
# ==============================================================================
cat("\n", rep("=", 70), "\n")
cat("GENERANDO REPORTE HTML DE RESUMEN\n")
cat(rep("=", 70), "\n\n")

cat("10. Creando reporte HTML...\n")

# Cargar datos necesarios para el reporte
cat("  Cargando datos para el reporte...\n")

# Cargar QC metrics
qc_metrics_file <- file.path(qc_dir, "QC_metrics_summary.csv")
if (file.exists(qc_metrics_file)) {
    qc_metrics_report <- read.csv(qc_metrics_file, stringsAsFactors = FALSE)
} else {
    qc_metrics_report <- qc_metrics
}

# Cargar DEGs summary
degs_summary_file <- file.path(degs_dir, "DEGs_summary.csv")
if (file.exists(degs_summary_file)) {
    degs_summary_report <- read.csv(degs_summary_file, stringsAsFactors = FALSE)
} else {
    degs_summary_report <- summary_stats
}

# Cargar sample info
sample_info_file <- file.path(degs_dir, "sample_info.csv")
if (file.exists(sample_info_file)) {
    sample_info_report <- read.csv(sample_info_file, stringsAsFactors = FALSE)
} else {
    sample_info_report <- sample_info
}

# Cargar DEGs
degs_file <- file.path(degs_dir, "all_DEGs.tsv")
if (file.exists(degs_file)) {
    degs_report <- read.table(degs_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
} else {
    degs_report <- degs
}

# Cargar up y down genes
up_genes_file <- file.path(degs_dir, "up_genes_complete.tsv")
down_genes_file <- file.path(degs_dir, "down_genes_complete.tsv")

if (file.exists(up_genes_file)) {
    up_genes_report <- read.table(up_genes_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
} else {
    up_genes_report <- up_genes
}

if (file.exists(down_genes_file)) {
    down_genes_report <- read.table(down_genes_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
} else {
    down_genes_report <- down_genes
}

# Cargar filtering summary
filtering_summary_file <- file.path(output_dir, "GSE57978_filtering_summary.csv")
if (file.exists(filtering_summary_file)) {
    filtering_summary_report <- read.csv(filtering_summary_file, stringsAsFactors = FALSE)
} else {
    filtering_summary_report <- filtering_summary
}

cat("  ✓ Datos cargados\n")
cat("  Generando contenido HTML...\n")

# Iniciar HTML
html_content <- paste0('<!DOCTYPE html>
<html lang="es">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Resumen de Análisis - GSE57978</title>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: "Segoe UI", Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
        }
        
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 10px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
            padding: 40px;
        }
        
        h1 {
            color: #667eea;
            text-align: center;
            margin-bottom: 10px;
            font-size: 2.5em;
            border-bottom: 3px solid #667eea;
            padding-bottom: 15px;
        }
        
        h2 {
            color: #764ba2;
            margin-top: 40px;
            margin-bottom: 20px;
            font-size: 1.8em;
            border-left: 5px solid #764ba2;
            padding-left: 15px;
        }
        
        h3 {
            color: #555;
            margin-top: 30px;
            margin-bottom: 15px;
            font-size: 1.4em;
        }
        
        h4 {
            color: #666;
            margin-top: 20px;
            margin-bottom: 10px;
            font-size: 1.2em;
        }
        
        .summary-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }
        
        .info-box {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            transition: transform 0.3s;
        }
        
        .info-box:hover {
            transform: translateY(-5px);
        }
        
        .metric-label {
            font-size: 0.9em;
            opacity: 0.9;
            margin-bottom: 10px;
        }
        
        .metric-value {
            font-size: 2em;
            font-weight: bold;
        }
        
        .metric-box {
            background: #f8f9fa;
            border-left: 4px solid #667eea;
            padding: 15px;
            margin: 15px 0;
            border-radius: 4px;
        }
        
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        
        th {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
        }
        
        td {
            padding: 10px;
            border-bottom: 1px solid #ddd;
        }
        
        tr:nth-child(even) {
            background: #f8f9fa;
        }
        
        tr:hover {
            background: #e9ecef;
        }
        
        .link-button {
            display: inline-block;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 10px 20px;
            margin: 5px;
            border-radius: 5px;
            text-decoration: none;
            transition: all 0.3s;
            font-weight: 500;
        }
        
        .link-button:hover {
            background: linear-gradient(135deg, #764ba2 0%, #667eea 100%);
            transform: scale(1.05);
            box-shadow: 0 4px 8px rgba(0,0,0,0.2);
        }
        
        .timestamp {
            background: #e9ecef;
            padding: 20px;
            border-radius: 8px;
            margin-top: 40px;
            text-align: center;
            color: #666;
        }
        
        .section {
            margin: 30px 0;
            padding: 20px;
            background: #f8f9fa;
            border-radius: 8px;
        }
        
        ul {
            margin-left: 20px;
            margin-top: 10px;
        }
        
        li {
            margin: 5px 0;
        }
        
        .warning {
            background: #fff3cd;
            border-left: 4px solid #ffc107;
            padding: 15px;
            margin: 15px 0;
            border-radius: 4px;
        }
        
        .success {
            background: #d4edda;
            border-left: 4px solid #28a745;
            padding: 15px;
            margin: 15px 0;
            border-radius: 4px;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>📊 Resumen de Análisis - GSE57978</h1>
        <div class="metric-box">
            <p><strong>Dataset:</strong> GSE57978</p>
            <p><strong>Plataforma:</strong> Affymetrix Human Gene 1.0 ST Array</p>
            <p><strong>Comparación:</strong> CBD vs Vehicle</p>
            <p><strong>Fecha de análisis:</strong> ', Sys.Date(), '</p>
        </div>')

# Sección de información de muestras
if (!is.null(sample_info_report) && nrow(sample_info_report) > 0) {
    html_content <- paste0(html_content, '
        <h2>🧬 Información de Muestras</h2>
        <div class="section">
            <p><strong>Total de muestras:</strong> ', nrow(sample_info_report), '</p>
            <table>
                <tr>
                    <th>Muestra</th>
                    <th>GSM ID</th>
                    <th>Línea Celular</th>
                    <th>Tratamiento</th>
                </tr>')
    
    for (i in 1:nrow(sample_info_report)) {
        html_content <- paste0(html_content, '<tr>
            <td>', sample_info_report$Sample[i], '</td>
            <td>', sample_info_report$GSM_ID[i], '</td>
            <td>', sample_info_report$Cell_Line[i], '</td>
            <td>', sample_info_report$Treatment[i], '</td>
        </tr>')
    }
    
    html_content <- paste0(html_content, '</table>
        </div>')
}

# Sección de QC
html_content <- paste0(html_content, '
        <h2>🔍 Análisis de Calidad (QC)</h2>')

if (!is.null(qc_metrics_report) && nrow(qc_metrics_report) > 0) {
    html_content <- paste0(html_content, '
        <div class="summary-grid">
            <div class="info-box">
                <div class="metric-label">Mediana de Intensidad (Log2)</div>
                <div class="metric-value">', round(mean(qc_metrics_report$Median_Intensity, na.rm = TRUE), 2), '</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Media de Intensidad (Log2)</div>
                <div class="metric-value">', round(mean(qc_metrics_report$Mean_Intensity, na.rm = TRUE), 2), '</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Desviación Estándar Promedio</div>
                <div class="metric-value">', round(mean(qc_metrics_report$SD_Intensity, na.rm = TRUE), 2), '</div>
            </div>
            <div class="info-box">
                <div class="metric-label">IQR Promedio</div>
                <div class="metric-value">', round(mean(qc_metrics_report$IQR_Intensity, na.rm = TRUE), 2), '</div>
            </div>
        </div>
        
        <h3>Métricas por Muestra</h3>
        <table>
            <tr>
                <th>Muestra</th>
                <th>Tratamiento</th>
                <th>Mediana Intensidad</th>
                <th>Media Intensidad</th>
                <th>SD</th>
                <th>IQR</th>
            </tr>')
    
    for (i in 1:nrow(qc_metrics_report)) {
        html_content <- paste0(html_content, '<tr>
            <td>', substr(qc_metrics_report$Sample[i], 1, 50), ifelse(nchar(qc_metrics_report$Sample[i]) > 50, "...", ""), '</td>
            <td>', qc_metrics_report$Treatment[i], '</td>
            <td>', round(qc_metrics_report$Median_Intensity[i], 3), '</td>
            <td>', round(qc_metrics_report$Mean_Intensity[i], 3), '</td>
            <td>', round(qc_metrics_report$SD_Intensity[i], 3), '</td>
            <td>', round(qc_metrics_report$IQR_Intensity[i], 3), '</td>
        </tr>')
    }
    
    html_content <- paste0(html_content, '</table>')
}

html_content <- paste0(html_content, '
        <h3>Gráficos de QC</h3>
        <p>
            <a href="QC/01_histogram_intensities.png" class="link-button" target="_blank">📊 Histogramas de Intensidades</a>
            <a href="QC/02_boxplot_intensities.png" class="link-button" target="_blank">📦 Boxplots</a>
            <a href="QC/03_MAplots.png" class="link-button" target="_blank">📈 MA-plots</a>
            <a href="QC/04_PCA.png" class="link-button" target="_blank">🔬 Análisis PCA</a>
            <a href="QC/05_clustering_hierarchical.png" class="link-button" target="_blank">🌳 Clustering Jerárquico</a>
            <a href="QC/06_heatmap_correlation.png" class="link-button" target="_blank">🔥 Heatmap Correlación</a>
        </p>
        <p>
            <a href="QC/QC_metrics_summary.csv" class="link-button">📥 Descargar Métricas CSV</a>
            <a href="QC/arrayQualityMetrics/index.html" class="link-button" target="_blank">📋 Reporte arrayQualityMetrics</a>
        </p>')

# Sección de filtrado
if (!is.null(filtering_summary_report) && nrow(filtering_summary_report) > 0) {
    html_content <- paste0(html_content, '
        <h2>🔬 Filtrado de Genes</h2>
        <div class="section">
            <table>
                <tr>
                    <th>Criterio</th>
                    <th>Número</th>
                    <th>Porcentaje</th>
                </tr>')
    
    for (i in 1:nrow(filtering_summary_report)) {
        html_content <- paste0(html_content, '<tr>
            <td>', filtering_summary_report$Criterio[i], '</td>
            <td>', filtering_summary_report$Numero[i], '</td>
            <td>', filtering_summary_report$Porcentaje[i], '%</td>
        </tr>')
    }
    
    html_content <- paste0(html_content, '</table>
        </div>')
}

# Sección de DEGs
html_content <- paste0(html_content, '
        <h2>📈 Análisis de Expresión Diferencial (DEGs)</h2>')

if (!is.null(degs_summary_report) && nrow(degs_summary_report) > 0) {
    total_genes <- degs_summary_report$Numero[degs_summary_report$Categoria == "Total genes analizados"]
    total_degs <- degs_summary_report$Numero[degs_summary_report$Categoria == "DEGs totales"]
    up_degs <- degs_summary_report$Numero[degs_summary_report$Categoria == "Up-regulated (CBD > Vehicle)"]
    down_degs <- degs_summary_report$Numero[degs_summary_report$Categoria == "Down-regulated (CBD < Vehicle)"]
    
    html_content <- paste0(html_content, '
        <div class="summary-grid">
            <div class="info-box">
                <div class="metric-label">Total Genes Analizados</div>
                <div class="metric-value">', total_genes, '</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Total DEGs</div>
                <div class="metric-value">', total_degs, '</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Up-regulated</div>
                <div class="metric-value">', up_degs, '</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Down-regulated</div>
                <div class="metric-value">', down_degs, '</div>
            </div>
        </div>
        
        <div class="section">
            <table>
                <tr>
                    <th>Categoría</th>
                    <th>Número</th>
                    <th>Porcentaje</th>
                </tr>')
    
    for (i in 1:nrow(degs_summary_report)) {
        html_content <- paste0(html_content, '<tr>
            <td>', degs_summary_report$Categoria[i], '</td>
            <td>', degs_summary_report$Numero[i], '</td>
            <td>', degs_summary_report$Porcentaje[i], '%</td>
        </tr>')
    }
    
    html_content <- paste0(html_content, '</table>
        </div>')
}

html_content <- paste0(html_content, '
        <h3>Gráficos de DEGs</h3>
        <p>
            <a href="DEGs/volcano_plot.png" class="link-button" target="_blank">🌋 Volcano Plot</a>
            <a href="DEGs/MA_plot.png" class="link-button" target="_blank">📊 MA-plot</a>')

if (!is.null(degs_report) && nrow(degs_report) > 0) {
    html_content <- paste0(html_content, '
            <a href="DEGs/heatmap_top_genes.png" class="link-button" target="_blank">🔥 Heatmap Top Genes</a>')
}

html_content <- paste0(html_content, '
        </p>')

# Mostrar top DEGs
if (!is.null(degs_report) && nrow(degs_report) > 0) {
    # Ordenar por adj.P.Val
    degs_sorted <- degs_report[order(degs_report$adj.P.Val), ]
    n_top <- min(20, nrow(degs_sorted))
    top_degs <- degs_sorted[1:n_top, ]
    
    html_content <- paste0(html_content, '
        <h3>Top ', n_top, ' DEGs (por FDR)</h3>
        <div class="section">
            <table>
                <tr>
                    <th>Gen</th>
                    <th>ENTREZ ID</th>
                    <th>Log2 FC</th>
                    <th>P-value</th>
                    <th>FDR (adj.P.Val)</th>
                    <th>Regulación</th>
                </tr>')
    
    for (i in 1:nrow(top_degs)) {
        regulation <- ifelse(top_degs$logFC[i] > 0, "Up", "Down")
        regulation_color <- ifelse(regulation == "Up", "#28a745", "#dc3545")
        
        symbol_col <- ifelse("SYMBOL" %in% colnames(top_degs), top_degs$SYMBOL[i], 
                            ifelse("GeneID" %in% colnames(top_degs), top_degs$GeneID[i], "N/A"))
        entrez_col <- ifelse("ENTREZID" %in% colnames(top_degs), top_degs$ENTREZID[i], "N/A")
        
        html_content <- paste0(html_content, '<tr>
            <td><strong>', symbol_col, '</strong></td>
            <td>', entrez_col, '</td>
            <td>', round(top_degs$logFC[i], 3), '</td>
            <td>', format(top_degs$P.Value[i], scientific = TRUE, digits = 3), '</td>
            <td>', format(top_degs$adj.P.Val[i], scientific = TRUE, digits = 3), '</td>
            <td style="color: ', regulation_color, '; font-weight: bold;">', regulation, '</td>
        </tr>')
    }
    
    html_content <- paste0(html_content, '</table>
        </div>')
    
    # Separar up y down
    if (!is.null(up_genes_report) && nrow(up_genes_report) > 0) {
        html_content <- paste0(html_content, '
        <h3>Genes Up-regulated (CBD > Vehicle)</h3>
        <div class="success">
            <p><strong>Total:</strong> ', nrow(up_genes_report), ' genes</p>
        </div>')
    }
    
    if (!is.null(down_genes_report) && nrow(down_genes_report) > 0) {
        html_content <- paste0(html_content, '
        <h3>Genes Down-regulated (CBD < Vehicle)</h3>
        <div class="success">
            <p><strong>Total:</strong> ', nrow(down_genes_report), ' genes</p>
        </div>')
    }
} else {
    html_content <- paste0(html_content, '
        <div class="warning">
            <p><strong>⚠️ Nota:</strong> No se encontraron DEGs con los umbrales seleccionados.</p>
            <p>Revisa los archivos de resultados completos para más detalles.</p>
        </div>')
}

html_content <- paste0(html_content, '
        <h3>Archivos de Resultados</h3>
        <p>
            <a href="DEGs/all_results.tsv" class="link-button">📄 Todos los Resultados (TSV)</a>
            <a href="DEGs/all_DEGs.tsv" class="link-button">📄 Todos los DEGs (TSV)</a>')

if (!is.null(up_genes_report) && nrow(up_genes_report) > 0) {
    html_content <- paste0(html_content, '
            <a href="DEGs/up_genes.txt" class="link-button">📄 Genes Up-regulated</a>
            <a href="DEGs/up_genes_complete.tsv" class="link-button">📄 Up-regulated Completo (TSV)</a>')
}

if (!is.null(down_genes_report) && nrow(down_genes_report) > 0) {
    html_content <- paste0(html_content, '
            <a href="DEGs/down_genes.txt" class="link-button">📄 Genes Down-regulated</a>
            <a href="DEGs/down_genes_complete.tsv" class="link-button">📄 Down-regulated Completo (TSV)</a>')
}

html_content <- paste0(html_content, '
        </p>')

# Sección de archivos generados
html_content <- paste0(html_content, '
        <h2>📁 Archivos Generados</h2>
        <div class="section">
            <h3>Estructura de Directorios</h3>
            <ul>
                <li><strong>QC/</strong> - Análisis de calidad y gráficos
                    <ul>
                        <li>Histogramas, boxplots, MA-plots</li>
                        <li>Análisis PCA y clustering</li>
                        <li>Heatmaps de correlación</li>
                        <li>Reporte arrayQualityMetrics</li>
                    </ul>
                </li>
                <li><strong>DEGs/</strong> - Resultados de expresión diferencial
                    <ul>
                        <li>Volcano plot, MA-plot, Heatmaps</li>
                        <li>Tablas de DEGs (TSV)</li>
                        <li>Listas de genes up/down-regulated</li>
                        <li>Resultados completos del análisis</li>
                    </ul>
                </li>
                <li><strong>Archivos RData/</strong> - Objetos R guardados
                    <ul>
                        <li>Datos crudos y normalizados</li>
                        <li>Objetos de QC y DEGs</li>
                        <li>Phenodata y anotaciones</li>
                    </ul>
                </li>
            </ul>
        </div>')

# Cerrar HTML
html_content <- paste0(html_content, '
        <div class="timestamp">
            <p><strong>✅ Análisis completado exitosamente</strong></p>
            <p>Fecha: ', Sys.Date(), ' | Hora: ', format(Sys.time(), "%H:%M:%S"), '</p>
            <p>Para más detalles, consulta los archivos en los directorios correspondientes.</p>
        </div>
    </div>
</body>
</html>')

# Guardar HTML
html_file <- file.path(output_dir, "GSE57978_analysis_summary.html")
writeLines(html_content, html_file, useBytes = TRUE)
cat("  ✓ Reporte HTML guardado: GSE57978_analysis_summary.html\n")
cat("  📍 Ubicación:", html_file, "\n")
cat("  🌐 Abre el archivo en tu navegador para ver el resumen completo\n\n")

cat(rep("=", 70), "\n")
cat("REPORTE HTML GENERADO EXITOSAMENTE\n")
cat(rep("=", 70), "\n\n")

