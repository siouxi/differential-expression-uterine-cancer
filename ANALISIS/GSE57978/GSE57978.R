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
                       "limma", "ggplot2", "gridExtra", "matrixStats")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) BiocManager::install(new_packages)
BiocManager::install("simpleaffy")
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
# 7. Guardar objetos para uso posterior
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

cat("\n¡Carga de datos completada!\n")
cat("El objeto 'raw_data' contiene los datos crudos de expresión\n")
cat("El objeto 'pheno_data' contiene toda la información del phenodata\n")
cat("El objeto 'pheno_subset' contiene el phenodata filtrado para las muestras cargadas\n")

