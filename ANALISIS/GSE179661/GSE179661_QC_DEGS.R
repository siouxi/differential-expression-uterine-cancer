# Análisis del dataset GSE179661
# RNA-seq (probablemente FPKM o TPM)
# Hepatocarcinoma (HCC) tratado con CBD vs Control
# Estudio: "A Novel Mechanism of Cannabidiol in Suppressing Hepatocellular Carcinoma 
#          by Inducing GSDME-Dependent Pyroptosis"

# ==============================================================================
# 1. Cargar librerías necesarias
# ==============================================================================
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Instalar paquetes de Bioconductor si no están instalados
required_packages <- c("Biobase", "GEOquery", "limma", "edgeR", 
                       "RColorBrewer", "gplots", "ggplot2", "gridExtra",
                       "matrixStats", "AnnotationDbi", "org.Hs.eg.db")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) BiocManager::install(new_packages)

# Cargar librerías
library(Biobase)
library(GEOquery)
library(limma)
library(edgeR)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(gridExtra)
library(matrixStats)
library(AnnotationDbi)
library(org.Hs.eg.db)

# ==============================================================================
# 2. Definir rutas de los archivos
# ==============================================================================
# Obtener la ruta base del proyecto desde la ubicación del script
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
    current_dir <- getwd()
    if (basename(current_dir) == "GSE179661") {
        script_dir <- current_dir
    } else {
        script_dir <- file.path(current_dir, "ANALISIS", "GSE179661")
        if (!dir.exists(script_dir)) {
            script_dir <- file.path(current_dir, "ANALISIS", "GSE179661")
        }
    }
}

# Obtener la ruta base del proyecto (subir dos niveles desde ANALISIS/GSE179661)
base_dir <- dirname(dirname(script_dir))
data_dir <- file.path(base_dir, "DATA", "GSE179661")

# Buscar archivos de expresión (pueden tener diferentes nombres)
possible_count_files <- c(
    file.path(data_dir, "GSE179661_Raw_gene_counts.txt.gz"),
    file.path(data_dir, "GSE179661_Raw_gene_counts.txt"),
    file.path(data_dir, "GSE179661_gene_expression_matrix.txt.gz"),
    file.path(data_dir, "GSE179661_gene_expression_matrix.txt"),
    file.path(data_dir, "GSE179661_FPKM.txt.gz"),
    file.path(data_dir, "GSE179661_FPKM.txt"),
    file.path(data_dir, "GSE179661_TPM.txt.gz"),
    file.path(data_dir, "GSE179661_TPM.txt")
)

counts_file <- NULL
for (f in possible_count_files) {
    if (file.exists(f)) {
        counts_file <- f
        break
    }
}

# Si no se encuentra, buscar cualquier archivo .txt o .txt.gz en el directorio
if (is.null(counts_file)) {
    all_files <- list.files(data_dir, pattern = "\\.txt(\\.gz)?$", full.names = TRUE)
    if (length(all_files) > 0) {
        # Priorizar archivos que contengan "count", "FPKM", "TPM" o "expression"
        priority_files <- all_files[grepl("count|FPKM|TPM|expression", basename(all_files), ignore.case = TRUE)]
        if (length(priority_files) > 0) {
            counts_file <- priority_files[1]
        } else {
            counts_file <- all_files[1]
        }
    }
}

pheno_file <- file.path(data_dir, "GSE179661_pheno_data.csv")

cat("Rutas configuradas:\n")
cat("  Base directory:", base_dir, "\n")
cat("  Data directory:", data_dir, "\n")
if (!is.null(counts_file)) {
    cat("  Expression file:", counts_file, "\n")
} else {
    cat("  ⚠ Expression file: NO ENCONTRADO\n")
}
cat("  Pheno file:", pheno_file, "\n\n")

# Verificar que las rutas existen
if (is.null(counts_file) || !file.exists(counts_file)) {
    stop(paste("El archivo de expresión no existe. Buscado en:", paste(possible_count_files, collapse = ", ")))
}
if (!file.exists(pheno_file)) {
    warning(paste("El archivo de phenodata no existe:", pheno_file))
    cat("  Continuando sin phenodata (se intentará extraer de los nombres de muestras)\n\n")
}

# ==============================================================================
# 3. Cargar phenodata
# ==============================================================================
cat("Cargando phenodata...\n")
pheno_data <- NULL
if (file.exists(pheno_file)) {
    pheno_data <- read.csv(pheno_file, stringsAsFactors = FALSE)
    cat("  Dimensiones del phenodata:", dim(pheno_data), "\n")
    cat("  Número de muestras:", nrow(pheno_data), "\n")
    cat("  Columnas disponibles:", paste(colnames(pheno_data)[1:min(10, ncol(pheno_data))], collapse = ", "), "\n\n")
    
    if ("geo_accession" %in% colnames(pheno_data)) {
        cat("  GSM IDs:\n")
        print(pheno_data$geo_accession)
    }
    if ("title" %in% colnames(pheno_data)) {
        cat("\n  Títulos de muestras:\n")
        print(pheno_data$title)
    }
    cat("\n")
} else {
    cat("  ⚠ Phenodata no encontrado, se usará información de nombres de archivos\n\n")
}

# ==============================================================================
# 4. Cargar matriz de expresión
# ==============================================================================
cat("Cargando matriz de expresión...\n")
cat("  Leyendo archivo:", basename(counts_file), "\n")

# Detectar si está comprimido
if (grepl("\\.gz$", counts_file)) {
    con <- gzfile(counts_file, "r")
} else {
    con <- file(counts_file, "r")
}

# Leer primeras líneas para detectar formato
first_lines <- readLines(con, n = 5)
close(con)

# Leer el archivo
raw_counts <- read.table(counts_file, 
                        header = TRUE, 
                        sep = "\t", 
                        stringsAsFactors = FALSE,
                        check.names = FALSE,
                        quote = "",
                        comment.char = "",
                        row.names = 1)

cat("  Dimensiones de la matriz cruda:", dim(raw_counts), "\n")
cat("  Primeras columnas (muestras):", paste(head(colnames(raw_counts), 5), collapse = ", "), "\n")
cat("  Primeros genes (row.names):", paste(head(rownames(raw_counts), 5), collapse = ", "), "\n\n")

# IDs de genes
gene_ids <- rownames(raw_counts)
exprs_matrix <- as.matrix(raw_counts)

cat("  Matriz de expresión extraída:\n")
cat("    Genes:", nrow(exprs_matrix), "\n")
cat("    Muestras:", ncol(exprs_matrix), "\n")
cat("    Rango de valores:", range(exprs_matrix, na.rm = TRUE), "\n\n")

# Detectar tipo de datos (conteos vs FPKM/TPM)
if (max(exprs_matrix, na.rm = TRUE) > 1000 && any(exprs_matrix == floor(exprs_matrix), na.rm = TRUE)) {
    data_type <- "counts"
    cat("  Tipo de datos detectado: Conteos crudos\n")
} else {
    data_type <- "normalized"
    cat("  Tipo de datos detectado: FPKM/TPM normalizados\n")
}

# ==============================================================================
# 5. Preparar datos para análisis
# ==============================================================================
cat("Preparando datos para análisis...\n")

# Crear objeto DGEList
dge <- DGEList(counts = exprs_matrix, genes = data.frame(GeneID = gene_ids))

cat("  Objeto DGEList creado\n")
cat("  Dimensiones:", dim(dge), "\n")
cat("  Total de conteos:", sum(dge$counts, na.rm = TRUE), "\n")
cat("  Promedio de conteos por gen:", mean(colSums(dge$counts, na.rm = TRUE)), "\n\n")

raw_data <- dge

# Asociar phenodata si está disponible
sample_names <- colnames(exprs_matrix)
if (!is.null(pheno_data) && "geo_accession" %in% colnames(pheno_data)) {
    cat("Intentando asociar phenodata...\n")
    
    # Intentar match entre nombres de columnas y GSM IDs
    sample_mapping <- match(sample_names, pheno_data$geo_accession)
    
    if (any(!is.na(sample_mapping))) {
        valid_indices <- which(!is.na(sample_mapping))
        pheno_subset <- pheno_data[sample_mapping[valid_indices], ]
        rownames(pheno_subset) <- sample_names[valid_indices]
        dge$samples <- cbind(dge$samples, pheno_subset[match(sample_names, rownames(pheno_subset)), ])
        cat("  ✓ Phenodata asociado para", length(valid_indices), "muestras\n\n")
    } else {
        cat("  ⚠ No se pudo hacer match directo, continuando sin phenodata asociado\n\n")
    }
}

# ==============================================================================
# 6. Resumen de datos cargados
# ==============================================================================
cat("\n", rep("=", 70), "\n")
cat("RESUMEN DE DATOS CARGADOS\n")
cat(rep("=", 70), "\n\n")

cat("Dataset: GSE179661\n")
cat("Estudio: A Novel Mechanism of Cannabidiol in Suppressing Hepatocellular Carcinoma\n")
cat("         by Inducing GSDME-Dependent Pyroptosis\n")
cat("Tipo: RNA-seq\n")
cat("Número de genes:", nrow(exprs_matrix), "\n")
cat("Número de muestras:", ncol(exprs_matrix), "\n\n")

cat("Estadísticas de expresión:\n")
cat("  Mínimo:", min(exprs_matrix, na.rm = TRUE), "\n")
cat("  Máximo:", max(exprs_matrix, na.rm = TRUE), "\n")
cat("  Mediana:", median(exprs_matrix, na.rm = TRUE), "\n")
cat("  Media:", mean(exprs_matrix, na.rm = TRUE), "\n")
cat("  Genes con expresión > 0 en todas las muestras:", sum(rowSums(exprs_matrix > 0, na.rm = TRUE) == ncol(exprs_matrix)), "\n")
cat("  Genes con expresión = 0 en todas las muestras:", sum(rowSums(exprs_matrix == 0, na.rm = TRUE) == ncol(exprs_matrix)), "\n\n")

cat("✓ Datos cargados exitosamente\n")
cat("\n", rep("=", 70), "\n\n")

# Definir directorio de salida
output_dir <- script_dir
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# 7. Filtrado de Genes para RNA-seq
# ==============================================================================
cat("\n", rep("=", 70), "\n")
cat("INICIANDO FILTRADO DE GENES\n")
cat(rep("=", 70), "\n\n")

cat("7. Aplicando filtros de calidad a los genes...\n")

exprs_matrix_unfiltered <- exprs_matrix
n_genes_before <- nrow(exprs_matrix)

cat("  Genes antes del filtrado:", n_genes_before, "\n\n")

# Transformar a log2 si es necesario
if (max(exprs_matrix, na.rm = TRUE) > 100) {
    log_exprs_for_filter <- log2(exprs_matrix + 1)
    cat("  Datos transformados a log2 para filtrado\n")
} else {
    log_exprs_for_filter <- exprs_matrix
    cat("  Datos ya en escala log2\n")
}

# 7.1. Filtrado por baja expresión
cat("\n7.1. Filtrado por baja expresión...\n")

# Calcular expresión promedio por gen (en escala log2 para el filtrado)
mean_expression_per_gene <- rowMeans(log_exprs_for_filter, na.rm = TRUE)

# Threshold: percentil 10 de expresión promedio (más estricto que percentil 5)
expression_threshold <- quantile(mean_expression_per_gene, probs = 0.10, na.rm = TRUE)
cat("  Threshold de expresión promedio (percentil 10):", round(expression_threshold, 3), "\n")

# Además, requerir que el gen tenga expresión > threshold en al menos 50% de las muestras
min_samples_with_expression <- ceiling(ncol(exprs_matrix) * 0.5)
cat("  Mínimo de muestras con expresión > threshold:", min_samples_with_expression, 
    "de", ncol(exprs_matrix), "\n")

# Filtro combinado: expresión promedio > threshold Y expresión > threshold en ≥50% muestras
samples_above_threshold <- rowSums(log_exprs_for_filter > expression_threshold, na.rm = TRUE)
genes_above_expression <- (mean_expression_per_gene > expression_threshold) & 
                          (samples_above_threshold >= min_samples_with_expression)

cat("  Genes con expresión promedio > threshold:", sum(mean_expression_per_gene > expression_threshold), "\n")
cat("  Genes con expresión > threshold en al menos", min_samples_with_expression, 
    "muestras:", sum(samples_above_threshold >= min_samples_with_expression), "\n")
cat("  Genes que pasan filtro de expresión (ambos criterios):", sum(genes_above_expression), "(", 
    round(sum(genes_above_expression)/n_genes_before*100, 1), "%)\n")

# 7.2. Filtrado por baja varianza
cat("\n7.2. Filtrado por baja varianza...\n")

# Calcular varianza de cada gen (en escala log2)
gene_variances <- apply(log_exprs_for_filter, 1, var, na.rm = TRUE)

# Threshold: quantil 0.25 (más estricto que 0.2) - eliminar genes con varianza muy baja
variance_threshold <- quantile(gene_variances, probs = 0.25, na.rm = TRUE)
cat("  Threshold de varianza (quantil 0.25):", round(variance_threshold, 6), "\n")
cat("  Rango de varianzas: [", round(min(gene_variances, na.rm = TRUE), 6), 
    ", ", round(max(gene_variances, na.rm = TRUE), 6), "]\n")

genes_above_variance <- gene_variances > variance_threshold
cat("  Genes con varianza > threshold:", sum(genes_above_variance), "(", 
    round(sum(genes_above_variance)/n_genes_before*100, 1), "%)\n")
cat("  Genes eliminados por baja varianza:", sum(!genes_above_variance), "(", 
    round(sum(!genes_above_variance)/n_genes_before*100, 1), "%)\n")

# 7.3. Filtrado de genes no expresados
cat("\n7.3. Filtrado de genes no expresados...\n")
genes_expressed <- rowSums(exprs_matrix > 0, na.rm = TRUE) > 0
n_genes_not_expressed <- sum(!genes_expressed)
cat("  Genes no expresados:", n_genes_not_expressed, 
    "(", round(n_genes_not_expressed/n_genes_before*100, 1), "%)\n")

# 7.4. Combinar filtros
cat("\n7.4. Combinando filtros...\n")

# Mostrar estadísticas de cada filtro
cat("  Resumen de filtros individuales:\n")
cat("    - Filtro de expresión: ", sum(genes_above_expression), " genes pasan\n", sep = "")
cat("    - Filtro de varianza: ", sum(genes_above_variance), " genes pasan\n", sep = "")
cat("    - Filtro de genes expresados: ", sum(genes_expressed), " genes pasan\n", sep = "")

# Combinar todos los filtros (AND lógico)
genes_passing_filters <- genes_above_expression & genes_above_variance & genes_expressed

n_genes_after <- sum(genes_passing_filters)
n_genes_removed <- n_genes_before - n_genes_after

cat("\n  Genes que pasan TODOS los filtros:", n_genes_after, 
    "(", round(n_genes_after/n_genes_before*100, 1), "%)\n")
cat("  Genes eliminados:", n_genes_removed, 
    "(", round(n_genes_removed/n_genes_before*100, 1), "%)\n")

# Mostrar desglose de genes eliminados por cada filtro
genes_removed_by_expression <- sum(!genes_above_expression)
genes_removed_by_variance <- sum(!genes_above_variance)
genes_removed_by_not_expressed <- sum(!genes_expressed)

cat("\n  Desglose de genes eliminados:\n")
cat("    - Eliminados por baja expresión:", genes_removed_by_expression, "\n")
cat("    - Eliminados por baja varianza:", genes_removed_by_variance, "\n")
cat("    - Eliminados por no expresados:", genes_removed_by_not_expressed, "\n")
cat("    - Nota: algunos genes pueden ser eliminados por múltiples criterios\n")

# Aplicar filtros a la matriz
cat("\n  Aplicando filtros a la matriz de expresión...\n")
exprs_matrix_filtered <- exprs_matrix[genes_passing_filters, ]
gene_ids_filtered <- gene_ids[genes_passing_filters]

cat("  Dimensiones antes del filtrado:", dim(exprs_matrix_unfiltered), "\n")
cat("  Dimensiones después del filtrado:", dim(exprs_matrix_filtered), "\n")

# Verificar que el filtrado funcionó
if (nrow(exprs_matrix_filtered) != n_genes_after) {
    warning("¡ADVERTENCIA! El número de filas en la matriz filtrada no coincide con el número de genes que pasan los filtros")
} else {
    cat("  ✓ Filtrado aplicado correctamente\n")
}

# Actualizar variables
exprs_matrix <- exprs_matrix_filtered
gene_ids <- gene_ids_filtered

# Resumen del filtrado
filtering_summary <- data.frame(
    Criterio = c("Total genes iniciales", 
                 "Expresión promedio > threshold (percentil 10)",
                 "Expresión > threshold en ≥50% muestras",
                 "Pasan filtro de expresión (ambos criterios)",
                 "Varianza > threshold (quantil 0.25)",
                 "Genes expresados (al menos 1 muestra)",
                 "Pasan TODOS los filtros",
                 "Genes eliminados"),
    Numero = c(n_genes_before,
               sum(mean_expression_per_gene > expression_threshold),
               sum(samples_above_threshold >= min_samples_with_expression),
               sum(genes_above_expression),
               sum(genes_above_variance),
               sum(genes_expressed),
               n_genes_after,
               n_genes_removed),
    Porcentaje = c(100,
                   round(sum(mean_expression_per_gene > expression_threshold)/n_genes_before*100, 1),
                   round(sum(samples_above_threshold >= min_samples_with_expression)/n_genes_before*100, 1),
                   round(sum(genes_above_expression)/n_genes_before*100, 1),
                   round(sum(genes_above_variance)/n_genes_before*100, 1),
                   round(sum(genes_expressed)/n_genes_before*100, 1),
                   round(n_genes_after/n_genes_before*100, 1),
                   round(n_genes_removed/n_genes_before*100, 1))
)

cat("\n  Resumen del filtrado:\n")
print(filtering_summary)

write.csv(filtering_summary,
          file = file.path(output_dir, "GSE179661_filtering_summary.csv"),
          row.names = FALSE)
cat("\n  ✓ Resumen de filtrado guardado: GSE179661_filtering_summary.csv\n")

# Actualizar objeto DGEList
if (exists("dge")) {
    dge$counts <- exprs_matrix
    dge$genes <- data.frame(GeneID = gene_ids)
}

cat("\n  ✓ Filtrado completado\n")
cat("\n", rep("=", 70), "\n\n")

# ==============================================================================
# 8. Análisis de Calidad (QC) para RNA-seq
# ==============================================================================
cat("\n", rep("=", 70), "\n")
cat("INICIANDO ANÁLISIS DE CALIDAD (QC)\n")
cat(rep("=", 70), "\n\n")

qc_dir <- file.path(output_dir, "QC")
if (!dir.exists(qc_dir)) {
    dir.create(qc_dir, recursive = TRUE)
}
cat("Directorio de QC creado:", qc_dir, "\n\n")

# Transformar a log2 para visualización
if (max(exprs_matrix, na.rm = TRUE) > 100) {
    log_exprs <- log2(exprs_matrix + 1)
} else {
    log_exprs <- exprs_matrix
}

# 8.1. Histogramas de expresión
cat("8.1. Generando histogramas de expresión...\n")
png(file.path(qc_dir, "01_histogram_expression.png"), 
    width = 12, height = 8, units = "in", res = 300)

n_samples <- ncol(log_exprs)
n_cols <- min(3, n_samples)
n_rows <- ceiling(n_samples / n_cols)
par(mfrow = c(n_rows, n_cols))

for (i in 1:n_samples) {
    hist(log_exprs[, i], 
         main = paste("Distribución de Expresión\n", sample_names[i]),
         xlab = "Log2 Expresión",
         ylab = "Frecuencia",
         breaks = 100,
         col = "lightblue",
         border = "black")
    abline(v = median(log_exprs[, i], na.rm = TRUE), 
           col = "red", lwd = 2, lty = 2)
}
dev.off()
cat("  ✓ Histogramas guardados: 01_histogram_expression.png\n")

# 8.2. Boxplots
cat("\n8.2. Generando boxplots...\n")
png(file.path(qc_dir, "02_boxplot_expression.png"), 
    width = 12, height = 8, units = "in", res = 300)

boxplot(log_exprs, 
        main = "Distribución de Expresión por Muestra",
        xlab = "Muestras",
        ylab = "Log2 Expresión",
        las = 2,
        col = brewer.pal(min(12, n_samples), "Set3"))
dev.off()
cat("  ✓ Boxplots guardados: 02_boxplot_expression.png\n")

# 8.3. PCA
cat("\n8.3. Realizando análisis PCA...\n")
pca_result <- prcomp(t(log_exprs), scale. = TRUE, center = TRUE)
pca_summary <- summary(pca_result)

png(file.path(qc_dir, "03_PCA.png"), 
    width = 12, height = 10, units = "in", res = 300)

par(mfrow = c(2, 2))

# Scree plot
plot(pca_result, main = "Scree Plot", type = "l")

# PC1 vs PC2
plot(pca_result$x[, 1], pca_result$x[, 2],
     main = "PCA: PC1 vs PC2",
     xlab = paste("PC1 (", round(pca_summary$importance[2, 1] * 100, 1), "%)"),
     ylab = paste("PC2 (", round(pca_summary$importance[2, 2] * 100, 1), "%)"),
     pch = 19, col = "blue", cex = 1.5)
text(pca_result$x[, 1], pca_result$x[, 2], 
     labels = sample_names, pos = 3, cex = 0.8)

# PC1 vs PC3
plot(pca_result$x[, 1], pca_result$x[, 3],
     main = "PCA: PC1 vs PC3",
     xlab = paste("PC1 (", round(pca_summary$importance[2, 1] * 100, 1), "%)"),
     ylab = paste("PC3 (", round(pca_summary$importance[2, 3] * 100, 1), "%)"),
     pch = 19, col = "red", cex = 1.5)
text(pca_result$x[, 1], pca_result$x[, 3], 
     labels = sample_names, pos = 3, cex = 0.8)

# Varianza explicada
barplot(pca_summary$importance[2, 1:min(10, ncol(pca_result$x))],
        main = "Varianza Explicada por Componentes",
        xlab = "Componente Principal",
        ylab = "Proporción de Varianza",
        col = "steelblue")

dev.off()
cat("  ✓ PCA guardado: 03_PCA.png\n")

# 8.4. Heatmap de correlación
cat("\n8.4. Generando heatmap de correlación...\n")
cor_matrix <- cor(log_exprs, use = "pairwise.complete.obs")

png(file.path(qc_dir, "04_heatmap_correlation.png"), 
    width = 10, height = 10, units = "in", res = 300)

heatmap.2(cor_matrix,
          main = "Matriz de Correlación entre Muestras",
          trace = "none",
          col = colorRampPalette(c("blue", "white", "red"))(100),
          margins = c(10, 10),
          cexRow = 0.8,
          cexCol = 0.8)

dev.off()
cat("  ✓ Heatmap de correlación guardado: 04_heatmap_correlation.png\n")

# 8.5. Métricas de calidad
cat("\n8.5. Calculando métricas de calidad...\n")
qc_metrics <- data.frame(
    Sample = sample_names,
    Total_Counts_Log10 = log10(colSums(exprs_matrix, na.rm = TRUE) + 1),
    Genes_Detected = colSums(exprs_matrix > 0, na.rm = TRUE),
    Genes_Detected_Percent = 100 * colSums(exprs_matrix > 0, na.rm = TRUE) / nrow(exprs_matrix),
    Median_Log2_Expression = apply(log_exprs, 2, median, na.rm = TRUE),
    Mean_Log2_Expression = apply(log_exprs, 2, mean, na.rm = TRUE),
    stringsAsFactors = FALSE
)

png(file.path(qc_dir, "06_QC_metrics.png"), 
    width = 14, height = 10, units = "in", res = 300)

par(mfrow = c(2, 3))

plot(qc_metrics$Total_Counts_Log10, 
     main = "Total Counts (Log10)",
     ylab = "Log10(Total Counts)",
     pch = 19, col = "blue")
text(qc_metrics$Total_Counts_Log10, labels = qc_metrics$Sample, pos = 3, cex = 0.7)

plot(qc_metrics$Genes_Detected_Percent,
     main = "Genes Detectados (%)",
     ylab = "Porcentaje",
     pch = 19, col = "green")
text(qc_metrics$Genes_Detected_Percent, labels = qc_metrics$Sample, pos = 3, cex = 0.7)

plot(qc_metrics$Median_Log2_Expression,
     main = "Mediana de Expresión (Log2)",
     ylab = "Log2 Expresión",
     pch = 19, col = "red")
text(qc_metrics$Median_Log2_Expression, labels = qc_metrics$Sample, pos = 3, cex = 0.7)

barplot(qc_metrics$Total_Counts_Log10, names.arg = qc_metrics$Sample,
        main = "Total Counts por Muestra", las = 2, cex.names = 0.7)

barplot(qc_metrics$Genes_Detected_Percent, names.arg = qc_metrics$Sample,
        main = "Genes Detectados (%)", las = 2, cex.names = 0.7)

barplot(qc_metrics$Median_Log2_Expression, names.arg = qc_metrics$Sample,
        main = "Mediana de Expresión", las = 2, cex.names = 0.7)

dev.off()
cat("  ✓ Métricas de calidad guardadas: 06_QC_metrics.png\n")

# Detectar outliers
outlier_flags <- rep(FALSE, nrow(qc_metrics))
counts_zscore <- abs(scale(qc_metrics$Total_Counts_Log10))
outlier_flags[counts_zscore > 2] <- TRUE
genes_zscore <- abs(scale(qc_metrics$Genes_Detected_Percent))
outlier_flags[genes_zscore > 2] <- TRUE
median_zscore <- abs(scale(qc_metrics$Median_Log2_Expression))
outlier_flags[median_zscore > 2] <- TRUE

qc_metrics$Potential_Outlier <- outlier_flags

if (sum(outlier_flags) > 0) {
    cat("  ⚠ Muestras potencialmente problemáticas:\n")
    print(qc_metrics[outlier_flags, c("Sample", "Total_Counts_Log10", "Genes_Detected_Percent", 
                       "Median_Log2_Expression")])
} else {
    cat("  ✓ No se detectaron outliers obvios\n")
}

write.csv(qc_metrics, 
          file = file.path(qc_dir, "QC_metrics_summary.csv"), 
          row.names = FALSE)

save(log_exprs, cor_matrix, qc_metrics, pca_result, pca_summary,
     file = file.path(qc_dir, "QC_objects.RData"))
cat("  ✓ Objetos de QC guardados: QC_objects.RData\n")

cat("\n", rep("=", 70), "\n")
cat("ANÁLISIS DE CALIDAD COMPLETADO\n")
cat(rep("=", 70), "\n\n")

# ==============================================================================
# 9. Guardar objetos cargados
# ==============================================================================
cat("Guardando objetos cargados...\n")

save(raw_data, file = file.path(output_dir, "GSE179661_raw_data.RData"))
cat("  ✓ Objeto raw_data guardado\n")

save(exprs_matrix, file = file.path(output_dir, "GSE179661_exprs_matrix.RData"))
write.csv(exprs_matrix, 
          file = file.path(output_dir, "GSE179661_exprs_matrix.csv"),
          row.names = TRUE)
cat("  ✓ Matriz de expresión guardada\n")

if (exists("exprs_matrix_unfiltered")) {
    save(exprs_matrix_unfiltered, file = file.path(output_dir, "GSE179661_exprs_matrix_unfiltered.RData"))
    write.csv(exprs_matrix_unfiltered, 
              file = file.path(output_dir, "GSE179661_exprs_matrix_unfiltered.csv"),
              row.names = TRUE)
    cat("  ✓ Matriz sin filtrar guardada\n")
}

if (exists("filtering_summary")) {
    save(filtering_summary, expression_threshold, variance_threshold, 
         min_samples_with_expression, n_genes_before, n_genes_after,
         file = file.path(output_dir, "GSE179661_filtering_objects.RData"))
    cat("  ✓ Objetos de filtrado guardados\n")
}

if (!is.null(pheno_data)) {
    save(pheno_data, file = file.path(output_dir, "GSE179661_phenodata.RData"))
    cat("  ✓ Phenodata guardado\n")
}

cat("\n✓ Carga de datos completada\n\n")

# ==============================================================================
# 10. Análisis de Expresión Diferencial (DEGs)
# ==============================================================================
cat("\n", rep("=", 70), "\n")
cat("INICIANDO ANÁLISIS DE EXPRESIÓN DIFERENCIAL\n")
cat(rep("=", 70), "\n\n")

degs_dir <- file.path(output_dir, "DEGs")
if (!dir.exists(degs_dir)) {
    dir.create(degs_dir, recursive = TRUE)
}
cat("Directorio de DEGs creado:", degs_dir, "\n\n")

# 10.1. Preparar diseño experimental
cat("10.1. Preparando diseño experimental...\n")

# Detectar tratamientos desde nombres de muestras
treatment_info <- rep(NA, length(sample_names))

for (i in 1:length(sample_names)) {
    sample_name <- sample_names[i]
    sample_name_upper <- toupper(sample_name)
    
    # Detectar CBD
    if (grepl("CBD|CANNABIDIOL", sample_name_upper)) {
        treatment_info[i] <- "CBD"
    } else if (grepl("CONTROL|CTRL|VEHICLE|MOCK|UNTREATED|NON-TREATMENT", sample_name_upper)) {
        treatment_info[i] <- "Control"
    } else if (grepl("HCC|HEPATOCARCINOMA|CANCER", sample_name_upper)) {
        # Si solo dice HCC, podría ser control o tratamiento, necesitamos más contexto
        if (grepl("CBD|TREAT", sample_name_upper)) {
            treatment_info[i] <- "CBD"
        } else {
            treatment_info[i] <- "Control"
        }
    }
}

# Intentar desde phenodata si hay NAs
if (any(is.na(treatment_info)) && !is.null(pheno_data)) {
    cat("  Algunas muestras no tienen tratamiento detectado, intentando desde phenodata...\n")
    if ("title" %in% colnames(pheno_data)) {
        for (i in which(is.na(treatment_info))) {
            # Buscar en phenodata
            idx <- which(grepl(sample_names[i], pheno_data$title, ignore.case = TRUE) |
                        grepl(sample_names[i], pheno_data$geo_accession, ignore.case = TRUE))
            if (length(idx) > 0) {
                title_upper <- toupper(pheno_data$title[idx[1]])
                if (grepl("CBD|CANNABIDIOL", title_upper)) {
                    treatment_info[i] <- "CBD"
                } else if (grepl("CONTROL|CTRL|VEHICLE|MOCK", title_upper)) {
                    treatment_info[i] <- "Control"
                }
            }
        }
    }
}

# Si aún hay NAs, asignar basado en posición o dividir en grupos
if (any(is.na(treatment_info))) {
    cat("  ⚠ ADVERTENCIA: Algunas muestras no tienen tratamiento asignado\n")
    cat("  Asignando basado en posición (primera mitad = Control, segunda mitad = CBD)\n")
    n_samples <- length(sample_names)
    treatment_info[is.na(treatment_info)] <- ifelse(
        which(is.na(treatment_info)) <= n_samples/2, "Control", "CBD"
    )
}

# Crear factor
treatment <- factor(treatment_info, levels = c("Control", "CBD"))
treatment <- droplevels(treatment)

cat("  Distribución de tratamientos:\n")
print(table(treatment))
cat("\n")

# Crear data frame de información de muestras
sample_info <- data.frame(
    Sample = sample_names,
    Treatment = treatment,
    stringsAsFactors = FALSE
)

if (!is.null(pheno_data) && "geo_accession" %in% colnames(pheno_data)) {
    sample_info$GSM_ID <- NA
    for (i in 1:nrow(sample_info)) {
        idx <- which(pheno_data$geo_accession == sample_names[i] | 
                    grepl(sample_names[i], pheno_data$title, ignore.case = TRUE))
        if (length(idx) > 0) {
            sample_info$GSM_ID[i] <- pheno_data$geo_accession[idx[1]]
        }
    }
}

cat("  Información de muestras:\n")
print(sample_info)
cat("\n")

# Verificar réplicas
replicates_table <- table(treatment)
min_replicates <- min(replicates_table)
cat("  Réplicas por tratamiento:\n")
print(replicates_table)
if (min_replicates < 2) {
    cat("  ⚠ ADVERTENCIA: Algunos tratamientos tienen menos de 2 réplicas\n")
}
cat("\n")

# Crear matriz de diseño
design <- model.matrix(~ treatment)
colnames(design) <- make.names(colnames(design))

cat("  Matriz de diseño:\n")
cat("    Dimensiones:", dim(design), "\n")
cat("    Columnas:", paste(colnames(design), collapse = ", "), "\n")
cat("    Rango:", qr(design)$rank, "\n")
if (qr(design)$rank < ncol(design)) {
    cat("  ⚠ ERROR: Matriz de diseño es singular\n")
} else {
    cat("  ✓ Matriz de diseño es de rango completo\n")
}
cat("\n")

# 10.2. Análisis con limma
cat("10.2. Realizando análisis de expresión diferencial con limma...\n")

# Transformar a log2 si es necesario
if (max(exprs_matrix, na.rm = TRUE) > 100) {
    cat("  Transformando a log2...\n")
    log_exprs_degs <- log2(exprs_matrix + 1)
} else {
    cat("  Datos ya en escala log2\n")
    log_exprs_degs <- exprs_matrix
}

cat("  Ajustando modelo lineal...\n")
fit <- lmFit(log_exprs_degs, design)
cat("  ✓ Modelo ajustado\n")

cat("  Aplicando eBayes...\n")
fit2 <- eBayes(fit, trend = TRUE, robust = TRUE)
cat("  ✓ eBayes aplicado\n\n")

# Extraer resultados (CBD vs Control)
coef_name <- colnames(design)[grepl("CBD|treatment", colnames(design), ignore.case = TRUE)]
if (length(coef_name) == 0) {
    coef_name <- colnames(design)[-1]  # Todas excepto intercepto
}

cat("  Comparación: CBD vs Control (coeficiente:", coef_name, ")\n\n")

results <- topTable(fit2, coef = coef_name, number = Inf, sort.by = "P")

# Agregar símbolos de genes
results$GeneID <- rownames(results)

# Intentar anotar con símbolos
if (require("org.Hs.eg.db", quietly = TRUE)) {
    tryCatch({
        # Intentar mapear diferentes tipos de IDs
        gene_symbols <- mapIds(org.Hs.eg.db,
                              keys = results$GeneID,
                              column = "SYMBOL",
                              keytype = "ENSEMBL",
                              multiVals = "first")
        if (all(is.na(gene_symbols))) {
            gene_symbols <- mapIds(org.Hs.eg.db,
                                  keys = results$GeneID,
                                  column = "SYMBOL",
                                  keytype = "ENTREZID",
                                  multiVals = "first")
        }
        if (all(is.na(gene_symbols))) {
            gene_symbols <- results$GeneID
        }
        results$SYMBOL <- gene_symbols
    }, error = function(e) {
        results$SYMBOL <- results$GeneID
    })
} else {
    results$SYMBOL <- results$GeneID
}

# Reordenar columnas
desired_cols <- c("GeneID", "SYMBOL", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
available_cols <- intersect(desired_cols, colnames(results))
other_cols <- setdiff(colnames(results), available_cols)
results <- results[, c(available_cols, other_cols)]

# Definir umbrales
fc_threshold <- 0.5  # log2 fold change
pval_threshold <- 0.05  # FDR

# Identificar DEGs
degs <- results[
    abs(results$logFC) >= fc_threshold & 
    results$adj.P.Val < pval_threshold,
]

up_genes <- degs[degs$logFC > 0, ]
down_genes <- degs[degs$logFC < 0, ]

cat("  Total DEGs (|logFC| >=", fc_threshold, "y adj.P.Val <", pval_threshold, "):", nrow(degs), "\n")
cat("  Up-regulated (CBD > Control):", nrow(up_genes), "\n")
cat("  Down-regulated (CBD < Control):", nrow(down_genes), "\n\n")

# 10.3. Análisis específico de genes de piroptosis
cat("10.3. Analizando genes relacionados con piroptosis...\n")

# Genes clave de piroptosis
pyroptosis_genes <- c("GSDME", "CASP3", "CASP1", "CASP4", "CASP5", "CASP11",
                      "GSDMD", "NLRP3", "IL1B", "IL18", "ATF4", "IGFBP1",
                      "AKT1", "AKT2", "AKT3", "CHOP", "DDIT3", "EIF2AK3",
                      "PERK", "ATF6", "XBP1", "IRE1", "ERN1")

# Buscar estos genes en los resultados
pyroptosis_results <- results[results$SYMBOL %in% pyroptosis_genes | 
                              results$GeneID %in% pyroptosis_genes, ]

if (nrow(pyroptosis_results) > 0) {
    cat("  Genes de piroptosis encontrados:", nrow(pyroptosis_results), "\n")
    cat("  Genes encontrados:", paste(pyroptosis_results$SYMBOL, collapse = ", "), "\n\n")
    
    # Guardar resultados de piroptosis
    write.csv(pyroptosis_results,
              file = file.path(degs_dir, "pyroptosis_genes_results.csv"),
              row.names = FALSE)
    cat("  ✓ Resultados de genes de piroptosis guardados\n\n")
} else {
    cat("  ⚠ No se encontraron genes de piroptosis con los símbolos esperados\n")
    cat("  Esto puede deberse a diferentes nomenclaturas de IDs\n\n")
}

# 10.4. Generar gráficos
cat("10.4. Generando gráficos de expresión diferencial...\n")

# Volcano Plot
png(file.path(degs_dir, "volcano_plot.png"), 
    width = 12, height = 10, units = "in", res = 300)

gene_labels <- ifelse(is.na(results$SYMBOL) | results$SYMBOL == "", 
                     results$GeneID, results$SYMBOL)

plot_data <- data.frame(
    logFC = results$logFC,
    neg_log10_pval = -log10(results$adj.P.Val + 1e-300),
    SYMBOL = gene_labels,
    is_DEG = abs(results$logFC) >= fc_threshold & results$adj.P.Val < pval_threshold,
    stringsAsFactors = FALSE
)

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
        title = "Volcano Plot: CBD vs Control (Hepatocarcinoma)",
        subtitle = paste("DEGs:", nrow(degs), "| Up:", nrow(up_genes), "| Down:", nrow(down_genes)),
        x = "Log2 Fold Change (CBD vs Control)",
        y = "-Log10 Adjusted P-value"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.title = element_text(size = 12),
        legend.position = "right"
    )

print(volcano_plot)
dev.off()
cat("  ✓ Volcano plot guardado\n")

# MA-plot
png(file.path(degs_dir, "MA_plot.png"), 
    width = 12, height = 10, units = "in", res = 300)

ma_plot <- ggplot(results, aes(x = AveExpr, y = logFC)) +
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
        title = "MA-plot: CBD vs Control (Hepatocarcinoma)",
        subtitle = "A = Average Expression, M = Log2 Fold Change",
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
cat("  ✓ MA-plot guardado\n")

# Heatmap de top DEGs
if (nrow(degs) > 0) {
    cat("\n10.5. Generando heatmap de top DEGs...\n")
    
    top_n <- min(50, nrow(degs))
    top_degs <- degs[order(degs$adj.P.Val), ][1:top_n, ]
    
    # Obtener expresión de estos genes
    top_gene_ids <- top_degs$GeneID
    top_exprs <- log_exprs_degs[rownames(log_exprs_degs) %in% top_gene_ids, , drop = FALSE]
    
    # Anotar con símbolos si es posible
    rownames(top_exprs) <- top_degs$SYMBOL[match(rownames(top_exprs), top_degs$GeneID)]
    
    # Crear colores para tratamientos
    treatment_colors <- c("Control" = "blue", "CBD" = "red")
    
    png(file.path(degs_dir, "heatmap_top_genes.png"), 
        width = 12, height = 14, units = "in", res = 300)
    
    # Escalar por filas (genes)
    top_exprs_scaled <- t(scale(t(top_exprs)))
    
    # Crear colores para tratamientos
    col_colors <- treatment_colors[as.character(treatment)]
    
    heatmap.2(top_exprs_scaled,
              main = paste("Top", top_n, "DEGs: CBD vs Control"),
              scale = "none",
              col = colorRampPalette(c("blue", "white", "red"))(100),
              trace = "none",
              margins = c(10, 10),
              cexRow = 0.6,
              cexCol = 0.8,
              ColSideColors = col_colors,
              key = TRUE,
              keysize = 1.2,
              density.info = "none")
    
    legend("topright", 
           legend = names(treatment_colors),
           fill = treatment_colors,
           title = "Treatment",
           cex = 0.8)
    
    dev.off()
    cat("  ✓ Heatmap de top DEGs guardado\n")
}

# 10.6. Guardar tablas de resultados
cat("\n10.6. Guardando tablas de resultados...\n")

write.table(results,
            file = file.path(degs_dir, "all_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (nrow(degs) > 0) {
    write.table(degs,
                file = file.path(degs_dir, "all_DEGs.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    up_gene_names <- ifelse(is.na(up_genes$SYMBOL) | up_genes$SYMBOL == "", 
                           up_genes$GeneID, up_genes$SYMBOL)
    write.table(up_gene_names,
                file = file.path(degs_dir, "up_genes.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    down_gene_names <- ifelse(is.na(down_genes$SYMBOL) | down_genes$SYMBOL == "", 
                             down_genes$GeneID, down_genes$SYMBOL)
    write.table(down_gene_names,
                file = file.path(degs_dir, "down_genes.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    write.table(up_genes,
                file = file.path(degs_dir, "up_genes_complete.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    write.table(down_genes,
                file = file.path(degs_dir, "down_genes_complete.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
}

summary_stats <- data.frame(
    Comparison = "CBD_vs_Control",
    Total_DEGs = nrow(degs),
    Up_regulated = nrow(up_genes),
    Down_regulated = nrow(down_genes),
    stringsAsFactors = FALSE
)

write.csv(summary_stats,
          file = file.path(degs_dir, "DEGs_summary.csv"),
          row.names = FALSE)

cat("\n  Resumen de análisis:\n")
print(summary_stats)

cat("\n", rep("=", 70), "\n")
cat("ANÁLISIS DE EXPRESIÓN DIFERENCIAL COMPLETADO\n")
cat(rep("=", 70), "\n\n")

# Guardar objetos
save(fit, fit2, results, degs, up_genes, down_genes, design, treatment, 
     sample_info, summary_stats, pyroptosis_results,
     file = file.path(degs_dir, "DEGs_objects.RData"))
cat("  ✓ Objetos de DEGs guardados\n\n")

cat("✓ Análisis completo finalizado\n")
cat("  Todos los resultados se encuentran en:", output_dir, "\n\n")

# ==============================================================================
# 11. Generar Reporte HTML de Resumen
# ==============================================================================
cat("\n", rep("=", 70), "\n")
cat("GENERANDO REPORTE HTML DE RESUMEN\n")
cat(rep("=", 70), "\n\n")

cat("11. Creando reporte HTML...\n")

# Crear contenido HTML
html_content <- paste0('<!DOCTYPE html>
<html lang="es">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Resumen del Análisis - GSE179661</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        h1 {
            color: #2c3e50;
            border-bottom: 3px solid #e74c3c;
            padding-bottom: 10px;
        }
        h2 {
            color: #34495e;
            margin-top: 30px;
            border-left: 4px solid #e74c3c;
            padding-left: 10px;
        }
        h3 {
            color: #555;
            margin-top: 20px;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }
        th, td {
            border: 1px solid #ddd;
            padding: 12px;
            text-align: left;
        }
        th {
            background-color: #e74c3c;
            color: white;
            font-weight: bold;
        }
        tr:nth-child(even) {
            background-color: #f2f2f2;
        }
        tr:hover {
            background-color: #ffe8e8;
        }
        .metric-box {
            background-color: #ecf0f1;
            border-left: 4px solid #e74c3c;
            padding: 15px;
            margin: 15px 0;
            border-radius: 5px;
        }
        .metric-label {
            font-weight: bold;
            color: #2c3e50;
        }
        .metric-value {
            font-size: 1.2em;
            color: #27ae60;
            margin-top: 5px;
        }
        .comparison-section {
            background-color: #f8f9fa;
            padding: 20px;
            margin: 20px 0;
            border-radius: 5px;
            border: 1px solid #dee2e6;
        }
        .link-button {
            display: inline-block;
            padding: 10px 20px;
            background-color: #e74c3c;
            color: white;
            text-decoration: none;
            border-radius: 5px;
            margin: 5px;
        }
        .link-button:hover {
            background-color: #c0392b;
        }
        .summary-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }
        .info-box {
            background-color: #fff;
            border: 1px solid #ddd;
            padding: 15px;
            border-radius: 5px;
        }
        .timestamp {
            color: #7f8c8d;
            font-style: italic;
            margin-top: 30px;
        }
        .highlight-box {
            background-color: #fff3cd;
            border-left: 4px solid #ffc107;
            padding: 15px;
            margin: 20px 0;
            border-radius: 5px;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>📊 Resumen del Análisis - GSE179661</h1>
        
        <div class="highlight-box">
            <strong>Estudio:</strong> A Novel Mechanism of Cannabidiol in Suppressing Hepatocellular Carcinoma 
            by Inducing GSDME-Dependent Pyroptosis<br>
            <strong>Contexto:</strong> Análisis de expresión génica en hepatocarcinoma (HCC) tratado con CBD (cannabidiol)
        </div>
        
        <div class="timestamp">
            Generado el: ', Sys.time(), '
        </div>
        
        <h2>📋 Información del Dataset</h2>
        <div class="summary-grid">
            <div class="info-box">
                <div class="metric-label">Dataset</div>
                <div class="metric-value">GSE179661</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Tipo de Datos</div>
                <div class="metric-value">RNA-seq</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Tipo de Datos Detectado</div>
                <div class="metric-value">', data_type, '</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Comparación</div>
                <div class="metric-value">CBD vs Control</div>
            </div>
        </div>
        
        <h2>🧪 Información de Muestras</h2>
        <div class="metric-box">
            <div class="metric-label">Total de Muestras:</div>
            <div class="metric-value">', ncol(exprs_matrix), '</div>
        </div>
        <div class="metric-box">
            <div class="metric-label">Total de Genes (después de filtrado):</div>
            <div class="metric-value">', nrow(exprs_matrix), '</div>
        </div>')

# Agregar sección de filtrado
if (exists("filtering_summary")) {
    html_content <- paste0(html_content, '
        <h2>🔧 Filtrado de Genes</h2>
        <div class="metric-box">
            <h3>Resumen del Filtrado</h3>
            <table>
                <tr>
                    <th>Criterio</th>
                    <th>Número de Genes</th>
                    <th>Porcentaje</th>
                </tr>')
    
    for (i in 1:nrow(filtering_summary)) {
        html_content <- paste0(html_content, '<tr>
            <td>', filtering_summary$Criterio[i], '</td>
            <td>', filtering_summary$Numero[i], '</td>
            <td>', filtering_summary$Porcentaje[i], '%</td>
        </tr>')
    }
    
    html_content <- paste0(html_content, '</table>
            <p>
                <a href="GSE179661_filtering_summary.csv" class="link-button">Descargar Resumen de Filtrado (CSV)</a>
            </p>
        </div>')
    
    if (exists("expression_threshold") && exists("variance_threshold")) {
        html_content <- paste0(html_content, '
        <div class="summary-grid">
            <div class="info-box">
                <div class="metric-label">Threshold de Expresión</div>
                <div class="metric-value">', round(expression_threshold, 3), '</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Threshold de Varianza</div>
                <div class="metric-value">', format(variance_threshold, scientific = TRUE, digits = 3), '</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Mínimo de Muestras con Expresión</div>
                <div class="metric-value">', if(exists("min_samples_with_expression")) min_samples_with_expression else ceiling(ncol(exprs_matrix) * 0.5), '</div>
            </div>
        </div>')
    }
}

# Agregar tabla de muestras
if (exists("sample_info")) {
    html_content <- paste0(html_content, '
        <h3>Distribución de Muestras por Tratamiento</h3>
        <table>
            <tr>
                <th>Muestra</th>
                <th>Tratamiento</th>')
    if ("GSM_ID" %in% colnames(sample_info)) {
        html_content <- paste0(html_content, '<th>GSM ID</th>')
    }
    html_content <- paste0(html_content, '</tr>')
    
    for (i in 1:nrow(sample_info)) {
        html_content <- paste0(html_content, '<tr><td>', sample_info$Sample[i], '</td><td>', 
                              sample_info$Treatment[i], '</td>')
        if ("GSM_ID" %in% colnames(sample_info) && !is.na(sample_info$GSM_ID[i])) {
            html_content <- paste0(html_content, '<td>', sample_info$GSM_ID[i], '</td>')
        }
        html_content <- paste0(html_content, '</tr>')
    }
    html_content <- paste0(html_content, '</table>')
    
    # Tabla de distribución
    treatment_table <- table(sample_info$Treatment)
    html_content <- paste0(html_content, '
        <h3>Resumen por Tratamiento</h3>
        <table>
            <tr><th>Tratamiento</th><th>Número de Muestras</th></tr>')
    for (i in 1:length(treatment_table)) {
        html_content <- paste0(html_content, '<tr><td>', names(treatment_table)[i], '</td><td>', 
                              treatment_table[i], '</td></tr>')
    }
    html_content <- paste0(html_content, '</table>')
}

# Agregar sección de QC
html_content <- paste0(html_content, '
        <h2>🔍 Análisis de Calidad (QC)</h2>')

if (exists("qc_metrics")) {
    html_content <- paste0(html_content, '
        <div class="summary-grid">
            <div class="info-box">
                <div class="metric-label">Promedio de Conteos Totales (Log10)</div>
                <div class="metric-value">', round(mean(qc_metrics$Total_Counts_Log10, na.rm = TRUE), 2), '</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Promedio de Genes Detectados (%)</div>
                <div class="metric-value">', round(mean(qc_metrics$Genes_Detected_Percent, na.rm = TRUE), 2), '%</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Mediana de Expresión (Log2)</div>
                <div class="metric-value">', round(median(qc_metrics$Median_Log2_Expression, na.rm = TRUE), 2), '</div>
            </div>
        </div>
        
        <h3>Gráficos de QC</h3>
        <p>
            <a href="QC/01_histogram_expression.png" class="link-button" target="_blank">Histogramas de Expresión</a>
            <a href="QC/02_boxplot_expression.png" class="link-button" target="_blank">Boxplots</a>
            <a href="QC/03_PCA.png" class="link-button" target="_blank">Análisis PCA</a>
            <a href="QC/04_heatmap_correlation.png" class="link-button" target="_blank">Heatmap Correlación</a>
            <a href="QC/06_QC_metrics.png" class="link-button" target="_blank">Métricas de Calidad</a>
        </p>
        <p>
            <a href="QC/QC_metrics_summary.csv" class="link-button">Descargar Métricas CSV</a>
        </p>')
}

# Agregar sección de DEGs
html_content <- paste0(html_content, '
        <h2>📈 Análisis de Expresión Diferencial (DEGs)</h2>')

if (exists("summary_stats") && nrow(summary_stats) > 0) {
    html_content <- paste0(html_content, '
        <h3>Resumen de Comparación</h3>
        <table>
            <tr>
                <th>Comparación</th>
                <th>Total DEGs</th>
                <th>Up-regulated</th>
                <th>Down-regulated</th>
            </tr>')
    
    html_content <- paste0(html_content, '<tr>
            <td><strong>CBD vs Control</strong></td>
            <td>', summary_stats$Total_DEGs[1], '</td>
            <td>', summary_stats$Up_regulated[1], '</td>
            <td>', summary_stats$Down_regulated[1], '</td>
        </tr>')
    
    html_content <- paste0(html_content, '</table>')
    
    # Sección de comparación CBD vs Control
    html_content <- paste0(html_content, '
        <div class="comparison-section">
            <h3>🔬 CBD vs Control (Hepatocarcinoma)</h3>
            <div class="summary-grid">
                <div class="info-box">
                    <div class="metric-label">Total DEGs</div>
                    <div class="metric-value">', nrow(degs), '</div>
                </div>
                <div class="info-box">
                    <div class="metric-label">Up-regulated</div>
                    <div class="metric-value">', nrow(up_genes), '</div>
                </div>
                <div class="info-box">
                    <div class="metric-label">Down-regulated</div>
                    <div class="metric-value">', nrow(down_genes), '</div>
                </div>
            </div>
            
            <h4>Gráficos</h4>
            <p>
                <a href="DEGs/volcano_plot.png" class="link-button" target="_blank">Volcano Plot</a>
                <a href="DEGs/MA_plot.png" class="link-button" target="_blank">MA-plot</a>')
    
    if (nrow(degs) > 0) {
        html_content <- paste0(html_content, '
                <a href="DEGs/heatmap_top_genes.png" class="link-button" target="_blank">Heatmap Top Genes</a>')
    }
    
    html_content <- paste0(html_content, '
            </p>
            
            <h4>Archivos de Resultados</h4>
            <p>
                <a href="DEGs/all_results.tsv" class="link-button">Todos los Resultados (TSV)</a>
                <a href="DEGs/all_DEGs.tsv" class="link-button">Todos los DEGs (TSV)</a>')
    
    if (nrow(up_genes) > 0) {
        html_content <- paste0(html_content, '
                <a href="DEGs/up_genes.txt" class="link-button">Genes Up-regulated</a>
                <a href="DEGs/down_genes.txt" class="link-button">Genes Down-regulated</a>')
    }
    
    html_content <- paste0(html_content, '
            </p>')
    
    # Mostrar top 10 DEGs
    if (nrow(degs) > 0) {
        top_degs_table <- degs[order(degs$adj.P.Val), ][1:min(10, nrow(degs)), ]
        html_content <- paste0(html_content, '
            <h4>Top 10 DEGs (por FDR)</h4>
            <table>
                <tr>
                    <th>Gen</th>
                    <th>Log2 FC</th>
                    <th>P-value</th>
                    <th>FDR</th>
                </tr>')
        
        for (j in 1:nrow(top_degs_table)) {
            html_content <- paste0(html_content, '<tr>
                <td>', ifelse("SYMBOL" %in% colnames(top_degs_table), top_degs_table$SYMBOL[j], top_degs_table$GeneID[j]), '</td>
                <td>', round(top_degs_table$logFC[j], 3), '</td>
                <td>', format(top_degs_table$P.Value[j], scientific = TRUE, digits = 3), '</td>
                <td>', format(top_degs_table$adj.P.Val[j], scientific = TRUE, digits = 3), '</td>
            </tr>')
        }
        html_content <- paste0(html_content, '</table>')
    }
    
    html_content <- paste0(html_content, '</div>')
}

# Agregar sección de genes de piroptosis
if (exists("pyroptosis_results") && nrow(pyroptosis_results) > 0) {
    html_content <- paste0(html_content, '
        <h2>🔥 Análisis de Genes de Piroptosis</h2>
        <div class="highlight-box">
            <strong>Genes clave analizados:</strong> GSDME, CASP3, CASP1, ATF4, IGFBP1, AKT, CHOP, y otros relacionados con piroptosis y estrés celular
        </div>
        <h3>Genes de Piroptosis Encontrados</h3>
        <table>
            <tr>
                <th>Gen</th>
                <th>Log2 FC</th>
                <th>P-value</th>
                <th>FDR</th>
                <th>Regulación</th>
            </tr>')
    
    for (i in 1:nrow(pyroptosis_results)) {
        regulation <- ifelse(pyroptosis_results$logFC[i] > 0, "Up", "Down")
        regulation_color <- ifelse(regulation == "Up", "#e74c3c", "#3498db")
        html_content <- paste0(html_content, '<tr>
            <td><strong>', ifelse("SYMBOL" %in% colnames(pyroptosis_results), pyroptosis_results$SYMBOL[i], pyroptosis_results$GeneID[i]), '</strong></td>
            <td>', round(pyroptosis_results$logFC[i], 3), '</td>
            <td>', format(pyroptosis_results$P.Value[i], scientific = TRUE, digits = 3), '</td>
            <td>', format(pyroptosis_results$adj.P.Val[i], scientific = TRUE, digits = 3), '</td>
            <td style="color: ', regulation_color, '; font-weight: bold;">', regulation, '</td>
        </tr>')
    }
    
    html_content <- paste0(html_content, '</table>
        <p>
            <a href="DEGs/pyroptosis_genes_results.csv" class="link-button">Descargar Resultados de Piroptosis (CSV)</a>
        </p>')
} else {
    html_content <- paste0(html_content, '
        <h2>🔥 Análisis de Genes de Piroptosis</h2>
        <div class="highlight-box">
            <strong>Nota:</strong> Los genes de piroptosis buscados (GSDME, CASP3, ATF4, etc.) no se encontraron con los símbolos esperados.
            Esto puede deberse a diferentes nomenclaturas de IDs de genes. Revisa el archivo de resultados completo para identificarlos.
        </div>')
}

# Cerrar HTML
html_content <- paste0(html_content, '
        <h2>📁 Archivos Generados</h2>
        <div class="metric-box">
            <h3>Estructura de Directorios</h3>
            <ul>
                <li><strong>QC/</strong> - Análisis de calidad y gráficos
                    <ul>
                        <li>Histogramas de expresión</li>
                        <li>Boxplots</li>
                        <li>Análisis PCA</li>
                        <li>Heatmap de correlación</li>
                        <li>Métricas de calidad</li>
                    </ul>
                </li>
                <li><strong>DEGs/</strong> - Resultados de expresión diferencial
                    <ul>
                        <li>Volcano plot</li>
                        <li>MA-plot</li>
                        <li>Heatmap de top genes</li>
                        <li>Tablas de resultados (TSV)</li>
                        <li>Listas de genes up/down regulados</li>
                        <li>Resultados de genes de piroptosis</li>
                    </ul>
                </li>
            </ul>
        </div>
        
        <div class="highlight-box">
            <h3>📝 Notas del Análisis</h3>
            <ul>
                <li><strong>Umbrales utilizados:</strong> |logFC| ≥ 0.5, FDR < 0.05</li>
                <li><strong>Método:</strong> limma con modelo lineal y eBayes</li>
                <li><strong>Filtrado:</strong> Genes con baja expresión, baja varianza o no expresados fueron eliminados</li>
                <li><strong>Piroptosis:</strong> Se buscaron genes clave relacionados con GSDME, caspasa-3, ATF4, y rutas de estrés celular</li>
            </ul>
        </div>
        
        <div class="timestamp">
            <p>Análisis completado exitosamente</p>
            <p>Para más detalles, consulta los archivos en los directorios correspondientes.</p>
        </div>
    </div>
</body>
</html>')

# Guardar HTML
html_file <- file.path(output_dir, "GSE179661_analysis_summary.html")
writeLines(html_content, html_file)
cat("  ✓ Reporte HTML guardado: GSE179661_analysis_summary.html\n")
cat("  Abre el archivo en tu navegador para ver el resumen completo\n\n")

# ==============================================================================
# OBJETOS DISPONIBLES:
# ==============================================================================
# - raw_data: Objeto DGEList con los datos crudos
# - exprs_matrix: Matriz de expresión (genes x muestras)
# - pheno_data: Phenodata completo (si disponible)
# - dge: Objeto DGEList con phenodata asociado (si fue posible)
# - results: Resultados completos de expresión diferencial
# - degs: DEGs identificados (CBD vs Control)
# - up_genes: Genes up-regulated
# - down_genes: Genes down-regulated
# - pyroptosis_results: Resultados de genes de piroptosis
# ==============================================================================

