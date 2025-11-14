# ==============================================================================
# Análisis del dataset GSE223827
# BRB-seq (3'-end, UMI) - Testículos fetales humanos expuestos a CBD/THC
# Estudio: "Cannabis exposure disrupts the endocannabinoid system in fetal testis"
# Pipeline especializada para BRB-seq con análisis pareado por donante
# ==============================================================================

cat("\n", rep("=", 80), "\n")
cat("ANÁLISIS GSE223827 - BRB-seq (Testículos Fetales Humanos)\n")
cat(rep("=", 80), "\n\n")

# ==============================================================================
# 1. Cargar librerías necesarias
# ==============================================================================
cat("1. Cargando librerías...\n")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Instalar paquetes de Bioconductor si no están instalados
required_bioc <- c("DESeq2", "limma", "edgeR", "Biobase", "GEOquery",
                   "AnnotationDbi", "org.Hs.eg.db", "clusterProfiler",
                   "enrichplot", "DOSE")
required_cran <- c("RColorBrewer", "gplots", "ggplot2", "gridExtra",
                   "matrixStats", "pheatmap", "uwot", "FactoMineR",
                   "dplyr", "tidyr", "viridis")

new_bioc <- required_bioc[!(required_bioc %in% installed.packages()[,"Package"])]
new_cran <- required_cran[!(required_cran %in% installed.packages()[,"Package"])]

if(length(new_bioc)) {
    cat("  Instalando paquetes Bioconductor:", paste(new_bioc, collapse=", "), "\n")
    BiocManager::install(new_bioc)
}
if(length(new_cran)) {
    cat("  Instalando paquetes CRAN:", paste(new_cran, collapse=", "), "\n")
    install.packages(new_cran)
}

install.packages("uwot")
library(uwot)



# Cargar librerías
library(DESeq2)
library(limma)
library(edgeR)
library(Biobase)
library(GEOquery)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(gridExtra)
library(matrixStats)
library(pheatmap)
library(uwot)
library(FactoMineR)
library(dplyr)
library(tidyr)
library(viridis)
library(AnnotationDbi)
library(org.Hs.eg.db)

cat("  ✓ Librerías cargadas\n\n")

# ==============================================================================
# 2. Definir rutas de los archivos
# ==============================================================================
cat("2. Configurando rutas...\n")

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
    if (basename(current_dir) == "GSE223827") {
        script_dir <- current_dir
    } else {
        script_dir <- file.path(current_dir, "ANALISIS", "GSE223827")
        if (!dir.exists(script_dir)) {
            script_dir <- file.path(current_dir, "ANALISIS", "GSE223827")
        }
    }
}

# Obtener la ruta base del proyecto
base_dir <- dirname(dirname(script_dir))
data_dir <- file.path(base_dir, "DATA", "GSE223827")
output_dir <- script_dir
qc_dir <- file.path(output_dir, "QC")
degs_dir <- file.path(output_dir, "DEGs")

# Crear directorios si no existen
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(degs_dir, showWarnings = FALSE, recursive = TRUE)

# Archivos de datos
counts_file <- file.path(data_dir, "GSE223827_BRBSeq.unq.refseq.umi.dat.txt.gz")
pheno_file <- file.path(data_dir, "GSE223827_pheno_data.csv")

cat("  Base directory:", base_dir, "\n")
cat("  Data directory:", data_dir, "\n")
cat("  Output directory:", output_dir, "\n")
cat("  Counts file:", basename(counts_file), "\n")
cat("  Pheno file:", basename(pheno_file), "\n\n")

# Verificar que los archivos existen
if (!file.exists(counts_file)) {
    stop(paste("El archivo de conteos no existe:", counts_file))
}
if (!file.exists(pheno_file)) {
    stop(paste("El archivo de phenodata no existe:", pheno_file))
}

# ==============================================================================
# 3. Cargar phenodata y procesar metadatos
# ==============================================================================
cat("3. Cargando y procesando metadatos...\n")

pheno_data <- read.csv(pheno_file, stringsAsFactors = FALSE, check.names = FALSE)
cat("  Dimensiones del phenodata:", dim(pheno_data), "\n")

# Extraer información relevante
# Las columnas importantes son: "exposure:ch1", "age (pcw):ch1", "developmental stage:ch1", "outlier:ch1"
sample_info <- data.frame(
    Sample = pheno_data$title,
    GSM_ID = pheno_data$geo_accession,
    Exposure = pheno_data[["exposure:ch1"]],
    Age_PCW = pheno_data[["age (pcw):ch1"]],
    DevStage = pheno_data[["developmental stage:ch1"]],
    Outlier = pheno_data[["outlier:ch1"]],
    stringsAsFactors = FALSE
)

# Filtrar outliers si es necesario
if (any(sample_info$Outlier == "yes", na.rm = TRUE)) {
    cat("  ⚠ Detectados", sum(sample_info$Outlier == "yes", na.rm = TRUE), "outliers\n")
    sample_info <- sample_info[sample_info$Outlier != "yes" | is.na(sample_info$Outlier), ]
    cat("  Muestras después de filtrar outliers:", nrow(sample_info), "\n")
}

# Crear variables de condición y edad
sample_info$Condition <- sample_info$Exposure
sample_info$Condition[sample_info$Condition == "DMSO"] <- "Control"

# Extraer edad en semanas (simplificar grupos)
sample_info$Age_Group <- ifelse(grepl("10-12GW|8-10PCW", sample_info$DevStage), "10-12GW", 
                                ifelse(grepl("12-14GW|10-12PCW", sample_info$DevStage), "12-14GW", 
                                       "Other"))

# Extraer donante del nombre de muestra (formato: testis_XX-XXGW_Cond_DMSO_EH####)
sample_info$Donor <- gsub(".*_EH([0-9]+)$", "EH\\1", sample_info$Sample)
sample_info$Donor[!grepl("^EH", sample_info$Donor)] <- paste0("Donor_", seq_along(sample_info$Donor[!grepl("^EH", sample_info$Donor)]))

cat("  Total de muestras:", nrow(sample_info), "\n")
cat("  Condiciones:", paste(unique(sample_info$Condition), collapse = ", "), "\n")
cat("  Grupos de edad:", paste(unique(sample_info$Age_Group), collapse = ", "), "\n")
cat("  Donantes únicos:", length(unique(sample_info$Donor)), "\n\n")

# Guardar sample_info
write.csv(sample_info, file = file.path(degs_dir, "sample_info.csv"), row.names = FALSE)
cat("  ✓ Metadatos guardados\n\n")

# ==============================================================================
# 4. Cargar matriz de conteos UMI
# ==============================================================================
cat("4. Cargando matriz de conteos UMI (BRB-seq)...\n")

# El archivo .dat tiene formato: gene_id \t sample1 \t sample2 \t ...
# Leer archivo comprimido
con <- gzfile(counts_file, "r")
first_lines <- readLines(con, n = 5)
close(con)

cat("  Primeras líneas del archivo:\n")
for (i in 1:min(3, length(first_lines))) {
    cat("    ", first_lines[i], "\n")
}

# Leer matriz de conteos
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
cat("  Primeros genes:", paste(head(rownames(raw_counts), 5), collapse = ", "), "\n")

# Convertir a matriz numérica
counts_matrix <- as.matrix(raw_counts)
counts_matrix[is.na(counts_matrix)] <- 0
counts_matrix <- round(counts_matrix)  # Asegurar enteros

cat("  Rango de valores:", range(counts_matrix, na.rm = TRUE), "\n")
cat("  Genes con expresión > 0:", sum(rowSums(counts_matrix) > 0), "\n")
cat("  Muestras con conteos > 0:", sum(colSums(counts_matrix) > 0), "\n\n")

# Filtrar muestras que no están en sample_info
common_samples <- intersect(colnames(counts_matrix), sample_info$Sample)
if (length(common_samples) < nrow(sample_info)) {
    cat("  ⚠ Algunas muestras del phenodata no están en la matriz de conteos\n")
    cat("    Muestras en phenodata:", nrow(sample_info), "\n")
    cat("    Muestras en matriz:", ncol(counts_matrix), "\n")
    cat("    Muestras comunes:", length(common_samples), "\n")
}

# Filtrar matriz y sample_info para mantener solo muestras comunes
counts_matrix <- counts_matrix[, common_samples, drop = FALSE]
sample_info <- sample_info[sample_info$Sample %in% common_samples, ]
sample_info <- sample_info[match(common_samples, sample_info$Sample), ]

cat("  ✓ Matriz filtrada:", dim(counts_matrix), "\n\n")

# ==============================================================================
# 5. Filtrado de genes (método del paper: mediana > 0)
# ==============================================================================
cat("5. Filtrado de genes (mediana > 0, según paper)...\n")

genes_before <- nrow(counts_matrix)
gene_medians <- rowMedians(counts_matrix)
genes_to_keep <- gene_medians > 0

counts_matrix <- counts_matrix[genes_to_keep, , drop = FALSE]
genes_after <- nrow(counts_matrix)

cat("  Genes antes del filtrado:", genes_before, "\n")
cat("  Genes después del filtrado (mediana > 0):", genes_after, "\n")
cat("  Genes eliminados:", genes_before - genes_after, 
    paste0("(", round(100 * (genes_before - genes_after) / genes_before, 2), "%)\n"))

# Guardar resumen de filtrado
filtering_summary <- data.frame(
    Criterio = c("Total genes", "Mediana > 0", "Genes finales"),
    Numero = c(genes_before, sum(genes_to_keep), genes_after),
    Porcentaje = c(100, round(100 * sum(genes_to_keep) / genes_before, 2), 
                   round(100 * genes_after / genes_before, 2))
)
write.csv(filtering_summary, file = file.path(output_dir, "GSE223827_filtering_summary.csv"), 
          row.names = FALSE)

cat("  ✓ Filtrado completado\n\n")

# ==============================================================================
# 6. Control de Calidad (QC)
# ==============================================================================
cat("6. Realizando control de calidad...\n")

# Calcular métricas de QC
qc_metrics <- data.frame(
    Sample = colnames(counts_matrix),
    Total_UMIs = colSums(counts_matrix),
    Genes_Detected = colSums(counts_matrix > 0),
    Median_Expression = colMedians(counts_matrix),
    Mean_Expression = colMeans(counts_matrix)
)

# Agregar información de muestra
qc_metrics <- merge(qc_metrics, sample_info, by = "Sample", all.x = TRUE)

# Guardar métricas
write.csv(qc_metrics, file = file.path(qc_dir, "QC_metrics_summary.csv"), row.names = FALSE)

cat("  ✓ Métricas de QC calculadas\n")

# Gráficos de QC
cat("  Generando gráficos de QC...\n")

# 6.1 Histograma de expresión
png(file.path(qc_dir, "01_histogram_expression.png"), width = 2000, height = 1200, res = 300)
par(mfrow = c(1, 2))
hist(log10(counts_matrix[counts_matrix > 0] + 1), 
     breaks = 100, 
     main = "Distribución de expresión (log10)",
     xlab = "log10(UMI counts + 1)",
     col = "steelblue",
     border = "white")
hist(rowSums(counts_matrix), 
     breaks = 100,
     main = "Distribución de suma de conteos por gen",
     xlab = "Total UMI counts",
     col = "coral",
     border = "white")
dev.off()
cat("    ✓ Histograma de expresión\n")

# 6.2 Boxplot de expresión por muestra
png(file.path(qc_dir, "02_boxplot_expression.png"), width = 3000, height = 1200, res = 300)
par(mar = c(12, 4, 2, 2))
boxplot(log10(counts_matrix + 1), 
        las = 2,
        main = "Distribución de expresión por muestra (log10)",
        ylab = "log10(UMI counts + 1)",
        col = ifelse(qc_metrics$Condition == "Control", "lightblue", "lightcoral"))
dev.off()
cat("    ✓ Boxplot de expresión\n")

# 6.3 Library size y genes detectados
png(file.path(qc_dir, "03_library_size.png"), width = 2000, height = 1200, res = 300)
par(mfrow = c(1, 2))
plot(qc_metrics$Total_UMIs, 
     qc_metrics$Genes_Detected,
     pch = 19,
     col = ifelse(qc_metrics$Condition == "Control", "blue", "red"),
     main = "Library Size vs Genes Detectados",
     xlab = "Total UMI counts",
     ylab = "Genes detectados")
legend("bottomright", 
       legend = c("Control", "Treatment"),
       col = c("blue", "red"),
       pch = 19)

barplot(qc_metrics$Total_UMIs,
        names.arg = qc_metrics$Sample,
        las = 2,
        cex.names = 0.6,
        main = "Library Size por Muestra",
        ylab = "Total UMI counts",
        col = ifelse(qc_metrics$Condition == "Control", "lightblue", "lightcoral"))
dev.off()
cat("    ✓ Library size\n")

cat("  ✓ Gráficos de QC generados\n\n")

# Guardar objetos de QC
save(counts_matrix, sample_info, qc_metrics, filtering_summary,
     file = file.path(qc_dir, "QC_objects.RData"))
cat("  ✓ Objetos de QC guardados\n\n")

# ==============================================================================
# 7. Normalización con rlog (DESeq2) - según paper
# ==============================================================================
cat("7. Normalización con rlog (DESeq2)...\n")

# Crear objeto DESeqDataSet
# Diseño: ~ Donor + Condition (pareado por donante)
dds <- DESeqDataSetFromMatrix(
    countData = counts_matrix,
    colData = sample_info,
    design = ~ Donor + Condition
)

# Pre-filtrado (genes con al menos 1 conteo en al menos 1 muestra)
keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep, ]

cat("  Genes después de pre-filtrado:", nrow(dds), "\n")

# Normalización con rlog (según paper)
rld <- rlog(dds, blind = FALSE)
normalized_matrix <- assay(rld)

cat("  Dimensiones de matriz normalizada:", dim(normalized_matrix), "\n")
cat("  Rango de valores normalizados:", range(normalized_matrix, na.rm = TRUE), "\n")
cat("  ✓ Normalización completada\n\n")

# Guardar matriz normalizada
write.csv(normalized_matrix, 
          file = file.path(output_dir, "GSE223827_normalized_matrix.csv"),
          row.names = TRUE)
save(rld, dds, normalized_matrix,
     file = file.path(output_dir, "GSE223827_normalization_objects.RData"))

# ==============================================================================
# 8. Análisis exploratorio (PCA, UMAP)
# ==============================================================================
cat("8. Análisis exploratorio (PCA, UMAP)...\n")

# 8.1 PCA
cat("  Calculando PCA...\n")
pca_result <- prcomp(t(normalized_matrix), scale. = TRUE, center = TRUE)
pca_variance <- summary(pca_result)$importance[2, 1:5] * 100

# Guardar resultados PCA
save(pca_result, pca_variance, file = file.path(qc_dir, "PCA_results.RData"))

# Gráfico PCA - por condición
png(file.path(qc_dir, "04_PCA_condition.png"), width = 2000, height = 1600, res = 300)
par(mfrow = c(2, 2))
plot(pca_result$x[, 1], pca_result$x[, 2],
     pch = 19,
     col = as.numeric(as.factor(sample_info$Condition)),
     main = "PCA: PC1 vs PC2 (por Condición)",
     xlab = paste0("PC1 (", round(pca_variance[1], 2), "%)"),
     ylab = paste0("PC2 (", round(pca_variance[2], 2), "%)"))
legend("topright", 
       legend = levels(as.factor(sample_info$Condition)),
       col = 1:length(levels(as.factor(sample_info$Condition))),
       pch = 19)

plot(pca_result$x[, 1], pca_result$x[, 3],
     pch = 19,
     col = as.numeric(as.factor(sample_info$Condition)),
     main = "PCA: PC1 vs PC3 (por Condición)",
     xlab = paste0("PC1 (", round(pca_variance[1], 2), "%)"),
     ylab = paste0("PC3 (", round(pca_variance[3], 2), "%)"))

plot(pca_result$x[, 1], pca_result$x[, 2],
     pch = 19,
     col = as.numeric(as.factor(sample_info$Age_Group)),
     main = "PCA: PC1 vs PC2 (por Edad)",
     xlab = paste0("PC1 (", round(pca_variance[1], 2), "%)"),
     ylab = paste0("PC2 (", round(pca_variance[2], 2), "%)"))
legend("topright", 
       legend = levels(as.factor(sample_info$Age_Group)),
       col = 1:length(levels(as.factor(sample_info$Age_Group))),
       pch = 19)

plot(pca_result$x[, 1], pca_result$x[, 2],
     pch = 19,
     col = as.numeric(as.factor(sample_info$Donor)),
     main = "PCA: PC1 vs PC2 (por Donante)",
     xlab = paste0("PC1 (", round(pca_variance[1], 2), "%)"),
     ylab = paste0("PC2 (", round(pca_variance[2], 2), "%)"))
dev.off()
cat("    ✓ PCA generado\n")

# 8.2 UMAP (si hay suficientes muestras) - usando uwot (más robusto, usado por Seurat)
if (ncol(normalized_matrix) >= 4) {
    cat("  Calculando UMAP (uwot)...\n")
    tryCatch({
        # uwot::umap() devuelve directamente una matriz (n_samples x 2)
        n_neighbors <- min(15, ncol(normalized_matrix) - 1)
        umap_result <- umap(t(normalized_matrix), 
                           n_neighbors = n_neighbors,
                           n_components = 2,
                           metric = "cosine",
                           verbose = FALSE)
        
        # uwot devuelve una matriz directamente, no un objeto con $layout
        png(file.path(qc_dir, "05_UMAP.png"), width = 2000, height = 1600, res = 300)
        par(mfrow = c(2, 2))
        plot(umap_result[, 1], umap_result[, 2],
             pch = 19,
             col = as.numeric(as.factor(sample_info$Condition)),
             main = "UMAP (por Condición)",
             xlab = "UMAP1", ylab = "UMAP2")
        legend("topright", 
               legend = levels(as.factor(sample_info$Condition)),
               col = 1:length(levels(as.factor(sample_info$Condition))),
               pch = 19)
        
        plot(umap_result[, 1], umap_result[, 2],
             pch = 19,
             col = as.numeric(as.factor(sample_info$Age_Group)),
             main = "UMAP (por Edad)",
             xlab = "UMAP1", ylab = "UMAP2")
        legend("topright", 
               legend = levels(as.factor(sample_info$Age_Group)),
               col = 1:length(levels(as.factor(sample_info$Age_Group))),
               pch = 19)
        
        plot(umap_result[, 1], umap_result[, 2],
             pch = 19,
             col = as.numeric(as.factor(sample_info$Donor)),
             main = "UMAP (por Donante)",
             xlab = "UMAP1", ylab = "UMAP2")
        legend("topright", 
               legend = levels(as.factor(sample_info$Donor)),
               col = 1:length(levels(as.factor(sample_info$Donor))),
               pch = 19,
               cex = 0.7)
        
        # Scree plot de varianza explicada (no aplica a UMAP, pero mostramos densidad)
        plot(density(umap_result[, 1]), 
             main = "Distribución UMAP1",
             xlab = "UMAP1",
             col = "steelblue",
             lwd = 2)
        dev.off()
        cat("    ✓ UMAP generado (uwot)\n")
    }, error = function(e) {
        cat("    ⚠ Error en UMAP:", e$message, "\n")
    })
}

# 8.3 Heatmap de correlación
cat("  Generando heatmap de correlación...\n")
cor_matrix <- cor(normalized_matrix)
png(file.path(qc_dir, "06_heatmap_correlation.png"), width = 2000, height = 2000, res = 300)
pheatmap(cor_matrix,
         main = "Matriz de Correlación entre Muestras",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation")
dev.off()
cat("    ✓ Heatmap de correlación\n")

cat("  ✓ Análisis exploratorio completado\n\n")

# ==============================================================================
# 9. Identificación de DEGs
# ==============================================================================
cat("9. Identificación de Genes Diferencialmente Expresados (DEGs)...\n")

# Definir comparaciones según el paper:
# - Por edad (10-12GW vs 12-14GW)
# - Por condición (CBD, THC, CBD+THC vs Control)
# - Pareado por donante

# Obtener condiciones únicas (excluyendo Control)
conditions <- unique(sample_info$Condition)
conditions <- conditions[conditions != "Control"]
age_groups <- unique(sample_info$Age_Group)
age_groups <- age_groups[age_groups != "Other"]

cat("  Condiciones a analizar:", paste(conditions, collapse = ", "), "\n")
cat("  Grupos de edad:", paste(age_groups, collapse = ", "), "\n\n")

# Almacenar todos los resultados
all_degs <- list()
all_results <- list()

# 9.1 Método del paper: Paired t-test con FC ≥ 1.3
cat("  9.1. Método del paper: Paired t-test (FC ≥ 1.3, BH adj p ≤ 0.05)...\n")

for (age_group in age_groups) {
    for (condition in conditions) {
        comp_name <- paste0(condition, "_vs_Control_", age_group)
        cat("    Analizando:", comp_name, "\n")
        
        # Filtrar muestras por edad y condición
        idx_age <- sample_info$Age_Group == age_group
        idx_cond <- sample_info$Condition %in% c("Control", condition)
        idx_keep <- idx_age & idx_cond
        
        if (sum(idx_keep) < 4) {
            cat("      ⚠ Muy pocas muestras, saltando comparación\n")
            next
        }
        
        sub_matrix <- normalized_matrix[, idx_keep, drop = FALSE]
        sub_info <- sample_info[idx_keep, ]
        
        # Verificar que hay muestras pareadas por donante
        donor_counts <- table(sub_info$Donor, sub_info$Condition)
        paired_donors <- rownames(donor_counts)[rowSums(donor_counts > 0) == 2]
        
        if (length(paired_donors) < 2) {
            cat("      ⚠ Pocos donantes pareados, usando t-test no pareado\n")
            paired <- FALSE
        } else {
            paired <- TRUE
            # Filtrar solo donantes pareados
            sub_info <- sub_info[sub_info$Donor %in% paired_donors, ]
            sub_matrix <- sub_matrix[, sub_info$Sample, drop = FALSE]
        }
        
        # Separar grupos
        control_idx <- sub_info$Condition == "Control"
        treat_idx <- sub_info$Condition == condition
        
        if (sum(control_idx) < 2 || sum(treat_idx) < 2) {
            cat("      ⚠ Muy pocas muestras por grupo, saltando\n")
            next
        }
        
        # Calcular fold change y estadísticas
        control_mean <- rowMeans(sub_matrix[, control_idx, drop = FALSE])
        treat_mean <- rowMeans(sub_matrix[, treat_idx, drop = FALSE])
        
        # FC en escala log2 (rlog ya está en log2)
        log2FC <- treat_mean - control_mean
        FC <- 2^log2FC
        
        # Paired t-test
        p_values <- numeric(nrow(sub_matrix))
        for (i in 1:nrow(sub_matrix)) {
            if (paired && length(paired_donors) >= 2) {
                # Organizar datos pareados
                control_vals <- sub_matrix[i, control_idx]
                treat_vals <- sub_matrix[i, treat_idx]
                
                # Intentar parear por donante
                if (length(control_vals) == length(treat_vals)) {
                    test_result <- tryCatch(
                        t.test(control_vals, treat_vals, paired = TRUE),
                        error = function(e) NULL
                    )
                    if (!is.null(test_result)) {
                        p_values[i] <- test_result$p.value
                    } else {
                        p_values[i] <- 1
                    }
                } else {
                    test_result <- tryCatch(
                        t.test(control_vals, treat_vals, paired = FALSE),
                        error = function(e) NULL
                    )
                    if (!is.null(test_result)) {
                        p_values[i] <- test_result$p.value
                    } else {
                        p_values[i] <- 1
                    }
                }
            } else {
                test_result <- tryCatch(
                    t.test(sub_matrix[i, control_idx], sub_matrix[i, treat_idx], paired = FALSE),
                    error = function(e) NULL
                )
                if (!is.null(test_result)) {
                    p_values[i] <- test_result$p.value
                } else {
                    p_values[i] <- 1
                }
            }
        }
        
        # Ajustar p-values (Benjamini-Hochberg)
        adj_p_values <- p.adjust(p_values, method = "BH")
        
        # Crear tabla de resultados
        results_table <- data.frame(
            GeneID = rownames(sub_matrix),
            log2FC = log2FC,
            FC = FC,
            Control_Mean = control_mean,
            Treat_Mean = treat_mean,
            PValue = p_values,
            AdjPValue = adj_p_values,
            stringsAsFactors = FALSE
        )
        
        # Filtrar según criterios del paper:
        # - Mediana > 0 (ya filtrado)
        # - |FC| ≥ 1.3
        # - AdjP ≤ 0.05
        results_table$Significant <- abs(results_table$FC) >= 1.3 & results_table$AdjPValue <= 0.05
        results_table$Direction <- ifelse(results_table$log2FC > 0, "Up", "Down")
        results_table$Direction[!results_table$Significant] <- "NS"
        
        # Separar DEGs
        degs_all <- results_table[results_table$Significant, ]
        degs_up <- degs_all[degs_all$log2FC > 0, ]
        degs_down <- degs_all[degs_all$log2FC < 0, ]
        
        cat("      Total DEGs (FC≥1.3, adjP≤0.05):", nrow(degs_all), "\n")
        cat("      Up-regulated:", nrow(degs_up), "\n")
        cat("      Down-regulated:", nrow(degs_down), "\n")
        
        # Guardar resultados
        all_results[[comp_name]] <- results_table
        all_degs[[comp_name]] <- list(
            all = degs_all,
            up = degs_up,
            down = degs_down
        )
        
        # Guardar archivos
        comp_dir <- file.path(degs_dir, comp_name)
        dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)
        
        write.table(results_table,
                   file = file.path(comp_dir, "all_results.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
        
        write.table(degs_all,
                   file = file.path(comp_dir, "all_DEGs.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
        
        write.table(degs_up$GeneID,
                   file = file.path(comp_dir, "up_genes.txt"),
                   sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        
        write.table(degs_down$GeneID,
                   file = file.path(comp_dir, "down_genes.txt"),
                   sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
}

cat("  ✓ Método del paper completado\n\n")

# 9.2 Método DESeq2 (más robusto)
cat("  9.2. Método DESeq2 (diseño pareado: ~ Donor + Condition)...\n")

# Ejecutar DESeq2
tryCatch({
    dds_de <- DESeq(dds)
    
    for (age_group in age_groups) {
        for (condition in conditions) {
            comp_name <- paste0(condition, "_vs_Control_", age_group, "_DESeq2")
            cat("    Analizando:", comp_name, "\n")
            
            # Filtrar muestras
            idx_age <- sample_info$Age_Group == age_group
            idx_cond <- sample_info$Condition %in% c("Control", condition)
            idx_keep <- idx_age & idx_cond
            
            if (sum(idx_keep) < 4) {
                cat("      ⚠ Muy pocas muestras, saltando\n")
                next
            }
            
            # Crear sub-DESeqDataSet
            dds_sub <- dds[, idx_keep]
            dds_sub$Condition <- relevel(dds_sub$Condition, ref = "Control")
            
            # Re-ejecutar DESeq en subconjunto
            dds_sub <- DESeq(dds_sub)
            
            # Obtener resultados
            res <- results(dds_sub, contrast = c("Condition", condition, "Control"))
            res_df <- as.data.frame(res)
            res_df$GeneID <- rownames(res_df)
            
            # Filtrar según criterios del paper
            res_df$Significant <- abs(res_df$log2FoldChange) >= log2(1.3) & 
                                  res_df$padj <= 0.05 & 
                                  !is.na(res_df$padj)
            res_df$Direction <- ifelse(res_df$log2FoldChange > 0, "Up", "Down")
            res_df$Direction[!res_df$Significant] <- "NS"
            
            # Separar DEGs
            degs_all <- res_df[res_df$Significant, ]
            degs_up <- degs_all[degs_all$log2FoldChange > 0, ]
            degs_down <- degs_all[degs_all$log2FoldChange < 0, ]
            
            cat("      Total DEGs (FC≥1.3, padj≤0.05):", nrow(degs_all), "\n")
            cat("      Up-regulated:", nrow(degs_up), "\n")
            cat("      Down-regulated:", nrow(degs_down), "\n")
            
            # Guardar
            all_results[[comp_name]] <- res_df
            all_degs[[comp_name]] <- list(
                all = degs_all,
                up = degs_up,
                down = degs_down
            )
            
            comp_dir <- file.path(degs_dir, comp_name)
            dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)
            
            write.table(res_df,
                       file = file.path(comp_dir, "all_results.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
            
            write.table(degs_all,
                       file = file.path(comp_dir, "all_DEGs.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
        }
    }
    
    cat("  ✓ DESeq2 completado\n\n")
}, error = function(e) {
    cat("  ⚠ Error en DESeq2:", e$message, "\n")
    cat("  Continuando con resultados del método del paper\n\n")
})

# Resumen de DEGs
summary_stats <- data.frame(
    Comparison = names(all_degs),
    Total_DEGs = sapply(all_degs, function(x) nrow(x$all)),
    Up_regulated = sapply(all_degs, function(x) nrow(x$up)),
    Down_regulated = sapply(all_degs, function(x) nrow(x$down)),
    stringsAsFactors = FALSE
)

write.csv(summary_stats,
          file = file.path(degs_dir, "DEGs_summary.csv"),
          row.names = FALSE)

cat("  Resumen de análisis:\n")
print(summary_stats)

# Guardar objetos
save(all_results, all_degs, summary_stats,
     file = file.path(degs_dir, "DEGs_objects.RData"))

cat("  ✓ DEGs identificados y guardados\n\n")

# ==============================================================================
# 10. Visualizaciones de DEGs
# ==============================================================================
cat("10. Generando visualizaciones de DEGs...\n")

for (comp_name in names(all_degs)) {
    if (nrow(all_degs[[comp_name]]$all) == 0) next
    
    cat("  Generando gráficos para:", comp_name, "\n")
    
    results_table <- all_results[[comp_name]]
    
    # Volcano plot
    png(file.path(degs_dir, paste0(comp_name, "_volcano_plot.png")), 
        width = 2000, height = 2000, res = 300)
    plot(results_table$log2FC, 
         -log10(results_table$AdjPValue),
         pch = 19,
         col = ifelse(results_table$Significant & results_table$log2FC > 0, "red",
                     ifelse(results_table$Significant & results_table$log2FC < 0, "blue", "gray")),
         main = paste("Volcano Plot -", comp_name),
         xlab = "log2 Fold Change",
         ylab = "-log10(Adjusted P-value)")
    abline(h = -log10(0.05), lty = 2, col = "black")
    abline(v = c(-log2(1.3), log2(1.3)), lty = 2, col = "black")
    legend("topright",
           legend = c("Up", "Down", "NS"),
           col = c("red", "blue", "gray"),
           pch = 19)
    dev.off()
    
    # MA plot
    mean_expr <- (results_table$Control_Mean + results_table$Treat_Mean) / 2
    if (all(is.na(mean_expr))) {
        mean_expr <- rowMeans(normalized_matrix[results_table$GeneID, ], na.rm = TRUE)
    }
    
    png(file.path(degs_dir, paste0(comp_name, "_MA_plot.png")), 
        width = 2000, height = 2000, res = 300)
    plot(mean_expr,
         results_table$log2FC,
         pch = 19,
         col = ifelse(results_table$Significant & results_table$log2FC > 0, "red",
                     ifelse(results_table$Significant & results_table$log2FC < 0, "blue", "gray")),
         main = paste("MA Plot -", comp_name),
         xlab = "Mean Expression",
         ylab = "log2 Fold Change")
    abline(h = c(-log2(1.3), log2(1.3)), lty = 2, col = "black")
    abline(h = 0, lty = 1, col = "black")
    legend("topright",
           legend = c("Up", "Down", "NS"),
           col = c("red", "blue", "gray"),
           pch = 19)
    dev.off()
    
    # Heatmap de top genes
    if (nrow(all_degs[[comp_name]]$all) > 0) {
        top_genes <- head(all_degs[[comp_name]]$all[order(abs(all_degs[[comp_name]]$all$log2FC), 
                                                          decreasing = TRUE), ], 50)
        if (nrow(top_genes) > 0) {
            # Obtener muestras relevantes
            idx_comp <- sample_info$Condition %in% c("Control", 
                                                     gsub("_vs_Control.*", "", comp_name))
            if (sum(idx_comp) > 0) {
                heatmap_data <- normalized_matrix[top_genes$GeneID, idx_comp, drop = FALSE]
                
                png(file.path(degs_dir, paste0(comp_name, "_heatmap_top_genes.png")), 
                    width = 2500, height = 2000, res = 300)
                pheatmap(heatmap_data,
                        main = paste("Top 50 DEGs -", comp_name),
                        color = colorRampPalette(c("blue", "white", "red"))(100),
                        scale = "row",
                        clustering_distance_rows = "correlation",
                        clustering_distance_cols = "correlation",
                        show_rownames = FALSE,
                        fontsize_col = 8)
                dev.off()
            }
        }
    }
}

cat("  ✓ Visualizaciones generadas\n\n")

# ==============================================================================
# 11. Generar Reporte HTML
# ==============================================================================
cat("11. Generando reporte HTML...\n")

# Cargar datos para el reporte
html_content <- paste0('<!DOCTYPE html>
<html lang="es">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Resumen del Análisis - GSE223827</title>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }
        .container {
            max-width: 1400px;
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
        img {
            max-width: 100%;
            height: auto;
            margin: 10px 0;
            border: 1px solid #ddd;
            border-radius: 5px;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>📊 Resumen del Análisis - GSE223827</h1>
        
        <div class="highlight-box">
            <strong>Estudio:</strong> Cannabis exposure disrupts the endocannabinoid system in fetal testis<br>
            <strong>Contexto:</strong> Análisis BRB-seq (3\'-end, UMI) de testículos fetales humanos expuestos a CBD/THC<br>
            <strong>Tipo de datos:</strong> BRB-seq (3\'-end sequencing con UMI deduplication)
        </div>
        
        <div class="timestamp">
            Generado el: ', Sys.time(), '
        </div>
        
        <h2>📋 Información del Dataset</h2>
        <div class="summary-grid">
            <div class="info-box">
                <div class="metric-label">Dataset</div>
                <div class="metric-value">GSE223827</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Tipo de Datos</div>
                <div class="metric-value">BRB-seq (UMI counts)</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Total de Muestras</div>
                <div class="metric-value">', nrow(sample_info), '</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Total de Genes</div>
                <div class="metric-value">', nrow(counts_matrix), '</div>
            </div>
        </div>
        
        <h2>🧪 Información de Muestras</h2>
        <div class="metric-box">
            <div class="metric-label">Condiciones:</div>
            <div class="metric-value">', paste(unique(sample_info$Condition), collapse = ", "), '</div>
        </div>
        <div class="metric-box">
            <div class="metric-label">Grupos de Edad:</div>
            <div class="metric-value">', paste(unique(sample_info$Age_Group), collapse = ", "), '</div>
        </div>
        <div class="metric-box">
            <div class="metric-label">Donantes Únicos:</div>
            <div class="metric-value">', length(unique(sample_info$Donor)), '</div>
        </div>')

# Tabla de muestras
html_content <- paste0(html_content, '
        <h3>Tabla de Muestras</h3>
        <table>
            <tr>
                <th>Muestra</th>
                <th>Condición</th>
                <th>Edad</th>
                <th>Donante</th>
                <th>Total UMIs</th>
                <th>Genes Detectados</th>
            </tr>')

for (i in 1:min(20, nrow(qc_metrics))) {
    html_content <- paste0(html_content, '<tr>
        <td>', qc_metrics$Sample[i], '</td>
        <td>', qc_metrics$Condition[i], '</td>
        <td>', qc_metrics$Age_Group[i], '</td>
        <td>', qc_metrics$Donor[i], '</td>
        <td>', round(qc_metrics$Total_UMIs[i], 0), '</td>
        <td>', qc_metrics$Genes_Detected[i], '</td>
    </tr>')
}

if (nrow(qc_metrics) > 20) {
    html_content <- paste0(html_content, '<tr><td colspan="6"><em>... y ', (nrow(qc_metrics) - 20), ' muestras más</em></td></tr>')
}

html_content <- paste0(html_content, '</table>')

# Filtrado
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
        </div>')

# QC
html_content <- paste0(html_content, '
        <h2>📊 Control de Calidad</h2>
        <div class="comparison-section">
            <h3>Gráficos de QC</h3>
            <p>
                <a href="QC/01_histogram_expression.png" class="link-button" target="_blank">📈 Histograma de Expresión</a>
                <a href="QC/02_boxplot_expression.png" class="link-button" target="_blank">📦 Boxplot de Expresión</a>
                <a href="QC/03_library_size.png" class="link-button" target="_blank">📚 Library Size</a>
                <a href="QC/04_PCA_condition.png" class="link-button" target="_blank">🔍 PCA</a>
                <a href="QC/06_heatmap_correlation.png" class="link-button" target="_blank">🔥 Heatmap Correlación</a>
            </p>
        </div>')

# DEGs
html_content <- paste0(html_content, '
        <h2>🧬 Genes Diferencialmente Expresados (DEGs)</h2>
        <div class="metric-box">
            <h3>Resumen de Comparaciones</h3>
            <table>
                <tr>
                    <th>Comparación</th>
                    <th>Total DEGs</th>
                    <th>Up-regulated</th>
                    <th>Down-regulated</th>
                </tr>')

for (i in 1:nrow(summary_stats)) {
    html_content <- paste0(html_content, '<tr>
        <td>', summary_stats$Comparison[i], '</td>
        <td>', summary_stats$Total_DEGs[i], '</td>
        <td>', summary_stats$Up_regulated[i], '</td>
        <td>', summary_stats$Down_regulated[i], '</td>
    </tr>')
}

html_content <- paste0(html_content, '</table>
        </div>')

# Métodos
html_content <- paste0(html_content, '
        <div class="highlight-box">
            <h3>Métodos de Análisis</h3>
            <ul>
                <li><strong>Método del Paper:</strong> Paired t-test con FC ≥ 1.3 y BH adj p ≤ 0.05</li>
                <li><strong>Método DESeq2:</strong> Diseño pareado (~ Donor + Condition) con FC ≥ 1.3 y padj ≤ 0.05</li>
                <li><strong>Normalización:</strong> rlog (DESeq2) según paper</li>
                <li><strong>Filtrado:</strong> Genes con mediana de expresión > 0</li>
            </ul>
        </div>')

# Archivos generados
html_content <- paste0(html_content, '
        <h2>📁 Archivos Generados</h2>
        <div class="section">
            <h3>Estructura de Directorios</h3>
            <ul>
                <li><strong>QC/</strong> - Análisis de calidad y gráficos
                    <ul>
                        <li>Histogramas, boxplots, library size</li>
                        <li>Análisis PCA</li>
                        <li>Heatmaps de correlación</li>
                        <li>Métricas de QC</li>
                    </ul>
                </li>
                <li><strong>DEGs/</strong> - Resultados de expresión diferencial
                    <ul>
                        <li>Volcano plots, MA-plots, Heatmaps</li>
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
html_file <- file.path(output_dir, "GSE223827_analysis_summary.html")
writeLines(html_content, html_file, useBytes = TRUE)
cat("  ✓ Reporte HTML guardado: GSE223827_analysis_summary.html\n")
cat("  📍 Ubicación:", html_file, "\n")
cat("  🌐 Abre el archivo en tu navegador para ver el resumen completo\n\n")

# ==============================================================================
# 12. Resumen final
# ==============================================================================
cat("\n", rep("=", 80), "\n")
cat("ANÁLISIS COMPLETADO EXITOSAMENTE\n")
cat(rep("=", 80), "\n\n")

cat("Resumen de archivos generados:\n")
cat("  - GSE223827_normalized_matrix.csv: Matriz normalizada (rlog)\n")
cat("  - GSE223827_filtering_summary.csv: Resumen del filtrado\n")
cat("  - sample_info.csv: Información de muestras\n")
cat("  - QC/: Directorio con todos los gráficos de calidad\n")
cat("  - DEGs/: Directorio con resultados de expresión diferencial\n")
cat("  - GSE223827_analysis_summary.html: Reporte HTML completo\n\n")

cat("Todos los resultados se encuentran en:", output_dir, "\n\n")

cat(rep("=", 80), "\n")
cat("FIN DEL ANÁLISIS\n")
cat(rep("=", 80), "\n\n")

