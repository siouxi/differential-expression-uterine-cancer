# Análisis del dataset GSE285498
# RNA-seq bulk (Illumina HiSeq4000)
# TopHat + Cufflinks → FPKM
# A549 cells tratadas con Mock, CBD, Etoposide, y combinación

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
    if (basename(current_dir) == "GSE285498") {
        script_dir <- current_dir
    } else {
        # Intentar construir la ruta relativa
        script_dir <- file.path(current_dir, "ANALISIS", "GSE285498")
        if (!dir.exists(script_dir)) {
            # Última opción: asumir que estamos en el directorio raíz del proyecto
            script_dir <- file.path(current_dir, "ANALISIS", "GSE285498")
        }
    }
}

# Obtener la ruta base del proyecto (subir dos niveles desde ANALISIS/GSE285498)
base_dir <- dirname(dirname(script_dir))
data_dir <- file.path(base_dir, "DATA", "GSE285498")
counts_file <- file.path(data_dir, "GSE285498_Raw_gene_counts.txt.gz")
pheno_file <- file.path(data_dir, "GSE285498_pheno_data.csv")

cat("Rutas configuradas:\n")
cat("  Base directory:", base_dir, "\n")
cat("  Data directory:", data_dir, "\n")
cat("  Counts file:", counts_file, "\n")
cat("  Pheno file:", pheno_file, "\n\n")

# Verificar que las rutas existen
if (!file.exists(counts_file)) {
    stop(paste("El archivo de conteos no existe:", counts_file))
}
if (!file.exists(pheno_file)) {
    stop(paste("El archivo de phenodata no existe:", pheno_file))
}

# ==============================================================================
# 3. Cargar phenodata
# ==============================================================================
cat("Cargando phenodata...\n")
pheno_data <- read.csv(pheno_file, stringsAsFactors = FALSE)

cat("  Dimensiones del phenodata:", dim(pheno_data), "\n")
cat("  Número de muestras:", nrow(pheno_data), "\n")
cat("  Columnas disponibles:", paste(colnames(pheno_data)[1:min(10, ncol(pheno_data))], collapse = ", "), "\n\n")

# Mostrar información básica de las muestras
if ("geo_accession" %in% colnames(pheno_data)) {
    cat("  GSM IDs:\n")
    print(pheno_data$geo_accession)
}
if ("treatment:ch1" %in% colnames(pheno_data)) {
    cat("\n  Tratamientos:\n")
    print(table(pheno_data[["treatment:ch1"]]))
}
cat("\n")

# ==============================================================================
# 4. Cargar matriz de expresión (conteos o FPKM)
# ==============================================================================
cat("Cargando matriz de expresión...\n")
cat("  Leyendo archivo:", basename(counts_file), "\n")

# Leer el archivo comprimido
# El archivo es tab-delimited con IDs de genes en la primera columna (que se convertirá en row.names)
raw_counts <- read.table(counts_file, 
                        header = TRUE, 
                        sep = "\t", 
                        stringsAsFactors = FALSE,
                        check.names = FALSE,
                        quote = "",
                        comment.char = "",
                        row.names = 1)  # La primera columna se usa como nombres de fila

cat("  Dimensiones de la matriz cruda:", dim(raw_counts), "\n")
cat("  Primeras columnas (muestras):", paste(head(colnames(raw_counts), 5), collapse = ", "), "\n")
cat("  Primeros genes (row.names):", paste(head(rownames(raw_counts), 5), collapse = ", "), "\n\n")

# Los IDs de genes están en los nombres de fila (row.names)
gene_ids <- rownames(raw_counts)

# La matriz de expresión son todas las columnas (todas son muestras)
exprs_matrix <- as.matrix(raw_counts)

# Los IDs de genes ya están como nombres de fila
cat("  IDs de genes extraídos de row.names:", length(gene_ids), "genes\n")

cat("  Matriz de expresión extraída:\n")
cat("    Genes:", nrow(exprs_matrix), "\n")
cat("    Muestras:", ncol(exprs_matrix), "\n")
cat("    Rango de valores:", range(exprs_matrix, na.rm = TRUE), "\n")
cat("    Valores únicos en primera fila:", length(unique(exprs_matrix[1, ])), "\n\n")

# Mostrar nombres de columnas (muestras)
cat("  Nombres de muestras en la matriz:\n")
print(colnames(exprs_matrix))
cat("\n")

# Verificar si los nombres de columnas coinciden con GSM IDs o Library names del phenodata
if ("geo_accession" %in% colnames(pheno_data)) {
    gsm_ids <- pheno_data$geo_accession
    matching_samples_gsm <- intersect(colnames(exprs_matrix), gsm_ids)
    cat("  Muestras que coinciden con GSM IDs:", length(matching_samples_gsm), "de", length(gsm_ids), "\n")
    
    # Intentar extraer Library names del phenodata
    library_names <- NULL
    if ("description" %in% colnames(pheno_data)) {
        # Extraer "Library name: XXX" de la columna description
        library_names <- gsub(".*Library name: ([^,]+).*", "\\1", pheno_data$description, ignore.case = TRUE)
        # Limpiar espacios y convertir a minúsculas para comparación
        library_names_clean <- tolower(trimws(library_names))
        matrix_names_clean <- tolower(trimws(colnames(exprs_matrix)))
        matching_samples_lib <- sum(library_names_clean %in% matrix_names_clean)
        cat("  Muestras que coinciden con Library names:", matching_samples_lib, "de", length(library_names), "\n")
    }
    
    # También intentar extraer del título (formato: "A549 cells, CBD_1")
    if ("title" %in% colnames(pheno_data)) {
        title_names <- gsub(".*,\\s*([^,]+)$", "\\1", pheno_data$title, ignore.case = TRUE)
        title_names_clean <- tolower(trimws(title_names))
        matrix_names_clean <- tolower(trimws(colnames(exprs_matrix)))
        matching_samples_title <- sum(title_names_clean %in% matrix_names_clean)
        cat("  Muestras que coinciden con nombres del título:", matching_samples_title, "de", length(title_names), "\n")
    }
    
    if (length(matching_samples_gsm) == 0 && 
        (is.null(library_names) || matching_samples_lib == 0) &&
        (!exists("matching_samples_title") || matching_samples_title == 0)) {
        cat("  ⚠ ADVERTENCIA: No se encontraron coincidencias directas\n")
        cat("    Intentando match flexible...\n")
    }
    cat("\n")
}

# ==============================================================================
# 5. Preparar datos para análisis
# ==============================================================================
cat("Preparando datos para análisis...\n")

# Crear objeto DGEList de edgeR (útil incluso si trabajamos con FPKM)
# Esto nos permite manejar los datos de manera consistente
dge <- DGEList(counts = exprs_matrix, genes = data.frame(GeneID = gene_ids))

cat("  Objeto DGEList creado\n")
cat("  Dimensiones:", dim(dge), "\n")
cat("  Total de conteos:", sum(dge$counts, na.rm = TRUE), "\n")
cat("  Promedio de conteos por gen:", mean(colSums(dge$counts, na.rm = TRUE)), "\n\n")

# Guardar objeto raw_data para referencia
raw_data <- dge

# Asociar phenodata con los datos
# Intentar hacer match entre nombres de columnas y diferentes identificadores
sample_mapping <- NULL
matching_method <- NULL

if ("geo_accession" %in% colnames(pheno_data)) {
    # Método 1: Intentar match directo con GSM IDs
    sample_mapping <- match(colnames(exprs_matrix), pheno_data$geo_accession)
    if (any(!is.na(sample_mapping)) && sum(!is.na(sample_mapping)) == length(colnames(exprs_matrix))) {
        matching_method <- "GSM_ID"
        cat("  ✓ Match encontrado usando GSM IDs\n")
    } else {
        # Método 2: Intentar con Library names (de description)
        if ("description" %in% colnames(pheno_data)) {
            library_names <- gsub(".*Library name: ([^,]+).*", "\\1", pheno_data$description, ignore.case = TRUE)
            library_names <- trimws(library_names)
            # Hacer match case-insensitive
            matrix_names_lower <- tolower(colnames(exprs_matrix))
            library_names_lower <- tolower(library_names)
            sample_mapping <- match(matrix_names_lower, library_names_lower)
            if (any(!is.na(sample_mapping)) && sum(!is.na(sample_mapping)) == length(colnames(exprs_matrix))) {
                matching_method <- "Library_name"
                cat("  ✓ Match encontrado usando Library names (description)\n")
            }
        }
        
        # Método 3: Intentar con nombres del título
        if (is.null(matching_method) && "title" %in% colnames(pheno_data)) {
            title_names <- gsub(".*,\\s*([^,]+)$", "\\1", pheno_data$title, ignore.case = TRUE)
            title_names <- trimws(title_names)
            # Hacer match case-insensitive
            matrix_names_lower <- tolower(colnames(exprs_matrix))
            title_names_lower <- tolower(title_names)
            sample_mapping <- match(matrix_names_lower, title_names_lower)
            if (any(!is.na(sample_mapping)) && sum(!is.na(sample_mapping)) == length(colnames(exprs_matrix))) {
                matching_method <- "Title"
                cat("  ✓ Match encontrado usando nombres del título\n")
            }
        }
        
        # Método 4: Match parcial (buscar patrones como "CBD_1" en cualquier parte)
        if (is.null(matching_method) || sum(!is.na(sample_mapping)) < length(colnames(exprs_matrix))) {
            cat("  Intentando match parcial...\n")
            sample_mapping <- rep(NA, length(colnames(exprs_matrix)))
            matrix_names_lower <- tolower(colnames(exprs_matrix))
            
            for (i in 1:length(matrix_names_lower)) {
                # Buscar en Library names
                if ("description" %in% colnames(pheno_data)) {
                    library_names <- gsub(".*Library name: ([^,]+).*", "\\1", pheno_data$description, ignore.case = TRUE)
                    library_names <- tolower(trimws(library_names))
                    idx <- which(library_names == matrix_names_lower[i])
                    if (length(idx) > 0) {
                        sample_mapping[i] <- idx[1]
                        next
                    }
                }
                # Buscar en títulos
                if ("title" %in% colnames(pheno_data)) {
                    title_names <- gsub(".*,\\s*([^,]+)$", "\\1", pheno_data$title, ignore.case = TRUE)
                    title_names <- tolower(trimws(title_names))
                    idx <- which(title_names == matrix_names_lower[i])
                    if (length(idx) > 0) {
                        sample_mapping[i] <- idx[1]
                        next
                    }
                }
                # Buscar patrón parcial (ej: "CBD_1" en "A549 cells, CBD_1")
                if ("title" %in% colnames(pheno_data)) {
                    for (j in 1:nrow(pheno_data)) {
                        if (grepl(matrix_names_lower[i], tolower(pheno_data$title[j]), fixed = TRUE)) {
                            sample_mapping[i] <- j
                            break
                        }
                    }
                }
            }
            
            if (sum(!is.na(sample_mapping)) > 0) {
                matching_method <- "Partial"
                cat("  ✓ Match parcial encontrado para", sum(!is.na(sample_mapping)), "muestras\n")
            }
        }
    }
    
    # Si encontramos algún match, asociar phenodata
    if (!is.null(sample_mapping) && any(!is.na(sample_mapping))) {
        # Filtrar phenodata para incluir solo las muestras presentes
        valid_indices <- which(!is.na(sample_mapping))
        pheno_subset <- pheno_data[sample_mapping[valid_indices], ]
        rownames(pheno_subset) <- colnames(exprs_matrix)[valid_indices]
        
        # Agregar phenodata al objeto DGEList
        dge$samples <- cbind(dge$samples, pheno_subset[match(colnames(dge$counts), rownames(pheno_subset)), ])
        
        cat("  ✓ Phenodata asociado con los datos (método:", matching_method, ")\n")
        cat("  Muestras con phenodata:", nrow(pheno_subset), "de", ncol(exprs_matrix), "\n")
        
        # Mostrar mapeo
        cat("\n  Mapeo de muestras:\n")
        for (i in valid_indices) {
            cat("    ", colnames(exprs_matrix)[i], " -> ", pheno_data$geo_accession[sample_mapping[i]], 
                " (", pheno_data$title[sample_mapping[i]], ")\n", sep = "")
        }
        cat("\n")
    } else {
        cat("  ⚠ No se pudo hacer match entre nombres de columnas y phenodata\n")
        cat("  Nombres en matriz:", paste(colnames(exprs_matrix), collapse = ", "), "\n")
        if ("description" %in% colnames(pheno_data)) {
            library_names <- gsub(".*Library name: ([^,]+).*", "\\1", pheno_data$description, ignore.case = TRUE)
            cat("  Library names en phenodata:", paste(trimws(library_names), collapse = ", "), "\n")
        }
        cat("  Continuando sin phenodata asociado\n\n")
    }
} else {
    cat("  ⚠ No se encontró columna 'geo_accession' en phenodata\n")
    cat("  Continuando sin phenodata asociado\n\n")
}

# ==============================================================================
# 6. Resumen de datos cargados
# ==============================================================================
cat("\n", rep("=", 70), "\n")
cat("RESUMEN DE DATOS CARGADOS\n")
cat(rep("=", 70), "\n\n")

cat("Dataset: GSE285498\n")
cat("Tipo: RNA-seq bulk (Illumina HiSeq4000)\n")
cat("Cuantificación: TopHat + Cufflinks → FPKM\n")
cat("Número de genes:", nrow(exprs_matrix), "\n")
cat("Número de muestras:", ncol(exprs_matrix), "\n\n")

if ("treatment:ch1" %in% colnames(pheno_data)) {
    cat("Distribución de tratamientos:\n")
    print(table(pheno_data[["treatment:ch1"]]))
    cat("\n")
}

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
# Objetivo: Eliminar genes con baja expresión, baja varianza o no informativos
# ==============================================================================
cat("\n", rep("=", 70), "\n")
cat("INICIANDO FILTRADO DE GENES\n")
cat(rep("=", 70), "\n\n")

cat("7. Aplicando filtros de calidad a los genes...\n")

# Guardar matriz original antes de filtrar
exprs_matrix_unfiltered <- exprs_matrix
n_genes_before <- nrow(exprs_matrix)

cat("  Genes antes del filtrado:", n_genes_before, "\n\n")

# Transformar a log2 si es necesario para los cálculos
if (max(exprs_matrix, na.rm = TRUE) > 100) {
    log_exprs_for_filter <- log2(exprs_matrix + 1)
    cat("  Datos transformados a log2 para filtrado\n")
} else {
    log_exprs_for_filter <- exprs_matrix
    cat("  Datos ya en escala log2\n")
}

# ------------------------------------------------------------------------------
# 7.1. Filtrado por baja expresión
# ------------------------------------------------------------------------------
cat("\n7.1. Filtrado por baja expresión...\n")

# Estimar threshold de expresión como percentil 5
expression_threshold <- quantile(exprs_matrix, probs = 0.05, na.rm = TRUE)
cat("  Threshold de expresión (percentil 5):", round(expression_threshold, 3), "\n")

# Número mínimo de muestras que deben tener expresión > threshold
# Por defecto: al menos 50% de las muestras
min_samples_with_expression <- ceiling(ncol(exprs_matrix) * 0.5)
cat("  Mínimo de muestras con expresión > threshold:", min_samples_with_expression, 
    "de", ncol(exprs_matrix), "\n")

# Contar cuántas muestras tienen expresión > threshold para cada gen
samples_above_threshold <- rowSums(exprs_matrix > expression_threshold, na.rm = TRUE)
genes_above_expression <- samples_above_threshold >= min_samples_with_expression

cat("  Genes con expresión > threshold en al menos", min_samples_with_expression, 
    "muestras:", sum(genes_above_expression), "(", 
    round(sum(genes_above_expression)/n_genes_before*100, 1), "%)\n")

# ------------------------------------------------------------------------------
# 7.2. Filtrado por baja varianza
# ------------------------------------------------------------------------------
cat("\n7.2. Filtrado por baja varianza...\n")

# Calcular varianza de cada gen (en escala log2)
gene_variances <- apply(log_exprs_for_filter, 1, var, na.rm = TRUE)

# Threshold de varianza: quantil 0.2 (eliminar genes con varianza en el 20% más bajo)
variance_threshold <- quantile(gene_variances, probs = 0.2, na.rm = TRUE)
cat("  Threshold de varianza (quantil 0.2):", round(variance_threshold, 6), "\n")

genes_above_variance <- gene_variances > variance_threshold
cat("  Genes con varianza > threshold:", sum(genes_above_variance), "(", 
    round(sum(genes_above_variance)/n_genes_before*100, 1), "%)\n")

# ------------------------------------------------------------------------------
# 7.3. Filtrado de genes no expresados en ninguna muestra
# ------------------------------------------------------------------------------
cat("\n7.3. Filtrado de genes no expresados...\n")

# Genes con expresión = 0 en todas las muestras
genes_expressed <- rowSums(exprs_matrix > 0, na.rm = TRUE) > 0
n_genes_not_expressed <- sum(!genes_expressed)
cat("  Genes no expresados (expresión = 0 en todas las muestras):", n_genes_not_expressed, 
    "(", round(n_genes_not_expressed/n_genes_before*100, 1), "%)\n")

# ------------------------------------------------------------------------------
# 7.4. Combinar todos los filtros
# ------------------------------------------------------------------------------
cat("\n7.4. Combinando filtros...\n")

# Genes que pasan TODOS los filtros
genes_passing_filters <- genes_above_expression & genes_above_variance & genes_expressed

n_genes_after <- sum(genes_passing_filters)
n_genes_removed <- n_genes_before - n_genes_after

cat("  Genes que pasan todos los filtros:", n_genes_after, 
    "(", round(n_genes_after/n_genes_before*100, 1), "%)\n")
cat("  Genes eliminados:", n_genes_removed, 
    "(", round(n_genes_removed/n_genes_before*100, 1), "%)\n")

# Aplicar filtros
exprs_matrix <- exprs_matrix[genes_passing_filters, ]
gene_ids <- gene_ids[genes_passing_filters]

cat("\n  Dimensiones antes del filtrado:", dim(exprs_matrix_unfiltered), "\n")
cat("  Dimensiones después del filtrado:", dim(exprs_matrix), "\n")

# Crear resumen del filtrado
filtering_summary <- data.frame(
    Criterio = c("Total genes iniciales", 
                 "Expresión > threshold (≥50% muestras)",
                 "Varianza > quantil 0.2",
                 "Genes expresados (al menos 1 muestra)",
                 "Pasan todos los filtros",
                 "Genes eliminados"),
    Numero = c(n_genes_before,
               sum(genes_above_expression),
               sum(genes_above_variance),
               sum(genes_expressed),
               n_genes_after,
               n_genes_removed),
    Porcentaje = c(100,
                   round(sum(genes_above_expression)/n_genes_before*100, 1),
                   round(sum(genes_above_variance)/n_genes_before*100, 1),
                   round(sum(genes_expressed)/n_genes_before*100, 1),
                   round(n_genes_after/n_genes_before*100, 1),
                   round(n_genes_removed/n_genes_before*100, 1))
)

cat("\n  Resumen del filtrado:\n")
print(filtering_summary)

# Guardar resumen de filtrado
write.csv(filtering_summary,
          file = file.path(output_dir, "GSE285498_filtering_summary.csv"),
          row.names = FALSE)
cat("\n  ✓ Resumen de filtrado guardado: GSE285498_filtering_summary.csv\n")

# Actualizar objeto DGEList con datos filtrados
if (exists("dge")) {
    dge$counts <- exprs_matrix
    dge$genes <- data.frame(GeneID = gene_ids)
}

cat("\n  ✓ Filtrado completado\n")
cat("\n", rep("=", 70), "\n\n")

# ==============================================================================
# 8. Análisis de Calidad (QC) para RNA-seq
# Objetivo: Detectar muestras problemáticas, outliers o problemas de secuenciación
# ==============================================================================
cat("\n", rep("=", 70), "\n")
cat("INICIANDO ANÁLISIS DE CALIDAD (QC)\n")
cat(rep("=", 70), "\n\n")

# Crear directorio para guardar gráficos de QC
qc_dir <- file.path(output_dir, "QC")
if (!dir.exists(qc_dir)) {
    dir.create(qc_dir, recursive = TRUE)
}
cat("Directorio de QC creado:", qc_dir, "\n\n")

# Obtener nombres de muestras
sample_names <- colnames(exprs_matrix)

# Transformar a log2 para visualización (agregar pseudocount pequeño para evitar log(0))
log_exprs <- log2(exprs_matrix + 1)

# ------------------------------------------------------------------------------
# 7.1. Distribución de Expresión (Histogramas)
# ------------------------------------------------------------------------------
cat("7.1. Generando histogramas de distribución de expresión...\n")

png(file.path(qc_dir, "01_histogram_expression.png"), 
    width = 14, height = 10, units = "in", res = 300)

n_samples <- ncol(log_exprs)
n_cols <- ceiling(sqrt(n_samples))
n_rows <- ceiling(n_samples / n_cols)

par(mfrow = c(n_rows, n_cols))
for (i in 1:n_samples) {
    hist(log_exprs[, i], 
         main = paste("Distribución de Expresión\n", sample_names[i]),
         xlab = "Log2(Expresión + 1)",
         ylab = "Frecuencia",
         breaks = 100,
         col = "lightblue",
         border = "black")
    abline(v = median(log_exprs[, i], na.rm = TRUE), 
           col = "red", lwd = 2, lty = 2)
    legend("topright", 
           legend = paste("Mediana =", round(median(log_exprs[, i], na.rm = TRUE), 2)),
           col = "red", lty = 2, lwd = 2,
           cex = 0.7)
}
dev.off()
cat("  ✓ Histograma guardado: 01_histogram_expression.png\n")

# ------------------------------------------------------------------------------
# 7.2. Boxplots de Expresión
# ------------------------------------------------------------------------------
cat("\n7.2. Generando boxplots de expresión...\n")

png(file.path(qc_dir, "02_boxplot_expression.png"), 
    width = 14, height = 8, units = "in", res = 300)

par(mar = c(8, 4, 4, 2))
boxplot(log_exprs, 
        main = "Distribución de Expresión por Muestra (Boxplots)",
        xlab = "",
        ylab = "Log2(Expresión + 1)",
        las = 2,
        col = brewer.pal(n = min(12, n_samples), name = "Set3"),
        cex.axis = 0.8)
mtext("Muestras", side = 1, line = 6)
abline(h = median(log_exprs, na.rm = TRUE), 
       col = "red", lwd = 2, lty = 2)
dev.off()
cat("  ✓ Boxplot guardado: 02_boxplot_expression.png\n")

# ------------------------------------------------------------------------------
# 7.3. Análisis PCA (Principal Component Analysis)
# ------------------------------------------------------------------------------
cat("\n7.3. Realizando análisis PCA...\n")

# Filtrar genes con baja varianza para PCA (más eficiente)
gene_variances <- apply(log_exprs, 1, var, na.rm = TRUE)
variance_threshold <- quantile(gene_variances, probs = 0.5, na.rm = TRUE)  # Top 50% más variables
high_var_genes <- gene_variances >= variance_threshold

cat("  Genes usados para PCA:", sum(high_var_genes), "de", nrow(log_exprs), "\n")

# Calcular PCA
pca_data <- t(log_exprs[high_var_genes, ])
pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
pca_summary <- summary(pca_result)

# Crear gráficos PCA
png(file.path(qc_dir, "03_PCA.png"), 
    width = 18, height = 6, units = "in", res = 300)

par(mfrow = c(1, 3))

# Obtener colores por tratamiento si está disponible
if (exists("pheno_subset") && "treatment:ch1" %in% colnames(pheno_subset)) {
    treatments <- pheno_subset[["treatment:ch1"]]
    treatment_colors <- brewer.pal(n = min(8, length(unique(treatments))), name = "Set2")
    names(treatment_colors) <- unique(treatments)
    sample_colors <- treatment_colors[treatments[match(sample_names, rownames(pheno_subset))]]
    sample_colors[is.na(sample_colors)] <- "gray"
} else {
    sample_colors <- brewer.pal(n = min(12, n_samples), name = "Set3")
}

# PCA Plot PC1 vs PC2
plot(pca_result$x[, 1], pca_result$x[, 2],
     main = "PCA: PC1 vs PC2",
     xlab = paste("PC1 (", round(pca_summary$importance[2, 1] * 100, 1), "%)"),
     ylab = paste("PC2 (", round(pca_summary$importance[2, 2] * 100, 1), "%)"),
     pch = 19,
     cex = 1.5,
     col = sample_colors)
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
     col = sample_colors)
text(pca_result$x[, 1], pca_result$x[, 3], 
     labels = sample_names, 
     pos = 3, cex = 0.7)
grid()

# Scree plot
n_pcs <- min(10, ncol(pca_result$x))
barplot(pca_summary$importance[2, 1:n_pcs],
        main = "Varianza Explicada por Componentes Principales",
        xlab = "Componente Principal",
        ylab = "Proporción de Varianza",
        col = "steelblue",
        names.arg = paste("PC", 1:n_pcs, sep = ""))

dev.off()
cat("  ✓ Gráficos PCA guardados: 03_PCA.png\n")

# Guardar resultados de PCA
save(pca_result, pca_summary, file = file.path(qc_dir, "PCA_results.RData"))

# ------------------------------------------------------------------------------
# 7.4. Heatmap de Correlación
# ------------------------------------------------------------------------------
cat("\n7.4. Generando heatmap de correlación...\n")

# Calcular matriz de correlación (usando correlación de Pearson)
cor_matrix <- cor(log_exprs, use = "pairwise.complete.obs")

# Heatmap de correlación (sin clustering jerárquico)
png(file.path(qc_dir, "04_heatmap_correlation.png"), 
    width = 12, height = 10, units = "in", res = 300)

heatmap.2(cor_matrix,
          main = "Matriz de Correlación entre Muestras",
          trace = "none",
          col = colorRampPalette(c("blue", "white", "red"))(100),
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          key = TRUE,
          keysize = 1.5,
          density.info = "none",
          margins = c(10, 10),
          cexRow = 0.8,
          cexCol = 0.8)
dev.off()
cat("  ✓ Heatmap de correlación guardado: 04_heatmap_correlation.png\n")

# ------------------------------------------------------------------------------
# 7.5. Métricas de Calidad para RNA-seq
# ------------------------------------------------------------------------------
cat("\n7.5. Generando resumen de métricas de calidad...\n")

# Calcular métricas específicas de RNA-seq
qc_metrics <- data.frame(
    Sample = sample_names,
    Total_Counts = colSums(exprs_matrix, na.rm = TRUE),
    Total_Counts_Log10 = log10(colSums(exprs_matrix, na.rm = TRUE) + 1),
    Genes_Detected = colSums(exprs_matrix > 0, na.rm = TRUE),
    Genes_Detected_Percent = 100 * colSums(exprs_matrix > 0, na.rm = TRUE) / nrow(exprs_matrix),
    Median_Expression = apply(exprs_matrix, 2, median, na.rm = TRUE),
    Mean_Expression = apply(exprs_matrix, 2, mean, na.rm = TRUE),
    Median_Log2_Expression = apply(log_exprs, 2, median, na.rm = TRUE),
    Mean_Log2_Expression = apply(log_exprs, 2, mean, na.rm = TRUE),
    SD_Log2_Expression = apply(log_exprs, 2, sd, na.rm = TRUE),
    IQR_Log2_Expression = apply(log_exprs, 2, IQR, na.rm = TRUE),
    stringsAsFactors = FALSE
)

# Agregar información del tratamiento si está disponible
if (exists("pheno_subset") && "treatment:ch1" %in% colnames(pheno_subset)) {
    qc_metrics$Treatment <- pheno_subset[["treatment:ch1"]][match(sample_names, rownames(pheno_subset))]
} else if ("geo_accession" %in% colnames(pheno_data)) {
    # Intentar extraer tratamiento del título
    qc_metrics$Treatment <- NA
    for (i in 1:nrow(qc_metrics)) {
        idx <- which(pheno_data$geo_accession == sample_names[i])
        if (length(idx) > 0 && "treatment:ch1" %in% colnames(pheno_data)) {
            qc_metrics$Treatment[i] <- pheno_data[["treatment:ch1"]][idx[1]]
        }
    }
}

# Guardar métricas
write.csv(qc_metrics, 
          file = file.path(qc_dir, "QC_metrics_summary.csv"), 
          row.names = FALSE)
cat("  ✓ Resumen de métricas guardado: QC_metrics_summary.csv\n")

# Mostrar resumen
cat("\nResumen de Métricas de Calidad:\n")
print(qc_metrics[, c("Sample", "Total_Counts", "Genes_Detected", "Genes_Detected_Percent", 
                     "Median_Log2_Expression")])

# ------------------------------------------------------------------------------
# 7.6. Gráficos de Métricas de Calidad
# ------------------------------------------------------------------------------
cat("\n7.6. Generando gráficos de métricas de calidad...\n")

png(file.path(qc_dir, "06_QC_metrics.png"), 
    width = 16, height = 10, units = "in", res = 300)

par(mfrow = c(2, 3))

# Total de conteos
barplot(qc_metrics$Total_Counts_Log10,
        main = "Total de Conteos por Muestra (Log10)",
        xlab = "Muestras",
        ylab = "Log10(Total Conteos + 1)",
        names.arg = sample_names,
        las = 2,
        cex.names = 0.7,
        col = brewer.pal(n = min(12, n_samples), name = "Set3"))

# Genes detectados
barplot(qc_metrics$Genes_Detected_Percent,
        main = "Porcentaje de Genes Detectados",
        xlab = "Muestras",
        ylab = "Porcentaje de Genes",
        names.arg = sample_names,
        las = 2,
        cex.names = 0.7,
        col = brewer.pal(n = min(12, n_samples), name = "Set3"))

# Mediana de expresión
barplot(qc_metrics$Median_Log2_Expression,
        main = "Mediana de Expresión (Log2)",
        xlab = "Muestras",
        ylab = "Log2(Mediana + 1)",
        names.arg = sample_names,
        las = 2,
        cex.names = 0.7,
        col = brewer.pal(n = min(12, n_samples), name = "Set3"))

# Scatter: Total conteos vs Genes detectados
plot(qc_metrics$Total_Counts_Log10, qc_metrics$Genes_Detected_Percent,
     main = "Total Conteos vs Genes Detectados",
     xlab = "Log10(Total Conteos + 1)",
     ylab = "Porcentaje de Genes Detectados",
     pch = 19,
     col = brewer.pal(n = min(12, n_samples), name = "Set3"))
text(qc_metrics$Total_Counts_Log10, qc_metrics$Genes_Detected_Percent,
     labels = sample_names, pos = 3, cex = 0.6)
grid()

# Scatter: Mediana vs Media
plot(qc_metrics$Median_Log2_Expression, qc_metrics$Mean_Log2_Expression,
     main = "Mediana vs Media de Expresión",
     xlab = "Mediana Log2(Expresión + 1)",
     ylab = "Media Log2(Expresión + 1)",
     pch = 19,
     col = brewer.pal(n = min(12, n_samples), name = "Set3"))
abline(0, 1, col = "red", lty = 2)
text(qc_metrics$Median_Log2_Expression, qc_metrics$Mean_Log2_Expression,
     labels = sample_names, pos = 3, cex = 0.6)
grid()

# Boxplot de métricas normalizadas
metrics_scaled <- scale(qc_metrics[, c("Total_Counts_Log10", "Genes_Detected_Percent", 
                                       "Median_Log2_Expression")])
boxplot(metrics_scaled,
        main = "Métricas de Calidad Normalizadas",
        names = c("Total Conteos", "Genes Detectados", "Mediana Expresión"),
        ylab = "Valor Normalizado (Z-score)",
        col = brewer.pal(3, "Set2"))

dev.off()
cat("  ✓ Gráficos de métricas guardados: 06_QC_metrics.png\n")

# ------------------------------------------------------------------------------
# 7.7. Detección de Outliers
# ------------------------------------------------------------------------------
cat("\n7.7. Detectando outliers...\n")

# Outliers basados en desviación de métricas clave
outlier_flags <- rep(FALSE, nrow(qc_metrics))

# Outliers en total de conteos
counts_zscore <- abs(scale(qc_metrics$Total_Counts_Log10))
outlier_flags[counts_zscore > 2] <- TRUE

# Outliers en genes detectados
genes_zscore <- abs(scale(qc_metrics$Genes_Detected_Percent))
outlier_flags[genes_zscore > 2] <- TRUE

# Outliers en mediana de expresión
median_zscore <- abs(scale(qc_metrics$Median_Log2_Expression))
outlier_flags[median_zscore > 2] <- TRUE

qc_metrics$Potential_Outlier <- outlier_flags

if (sum(outlier_flags) > 0) {
    cat("  ⚠ Muestras potencialmente problemáticas (desviación > 2 SD en alguna métrica):\n")
    outliers <- qc_metrics[outlier_flags, ]
    print(outliers[, c("Sample", "Total_Counts_Log10", "Genes_Detected_Percent", 
                       "Median_Log2_Expression")])
} else {
    cat("  ✓ No se detectaron outliers obvios basados en métricas de calidad\n")
}

# Guardar métricas actualizadas
write.csv(qc_metrics, 
          file = file.path(qc_dir, "QC_metrics_summary.csv"), 
          row.names = FALSE)

# Guardar objetos de QC
save(log_exprs, cor_matrix, qc_metrics, pca_result, pca_summary,
     file = file.path(qc_dir, "QC_objects.RData"))
cat("  ✓ Objetos de QC guardados: QC_objects.RData\n")

cat("\n", rep("=", 70), "\n")
cat("ANÁLISIS DE CALIDAD COMPLETADO\n")
cat(rep("=", 70), "\n")
cat("Todos los gráficos y métricas se han guardado en:", qc_dir, "\n\n")

# ==============================================================================
# 8. Guardar objetos cargados
# ==============================================================================
cat("Guardando objetos cargados...\n")

# Guardar objeto DGEList
save(raw_data, file = file.path(output_dir, "GSE285498_raw_data.RData"))
cat("  ✓ Objeto raw_data guardado: GSE285498_raw_data.RData\n")

# Guardar matriz de expresión (filtrada)
save(exprs_matrix, file = file.path(output_dir, "GSE285498_exprs_matrix.RData"))
write.csv(exprs_matrix, 
          file = file.path(output_dir, "GSE285498_exprs_matrix.csv"),
          row.names = TRUE)
cat("  ✓ Matriz de expresión guardada: GSE285498_exprs_matrix.csv\n")

# Guardar también matriz sin filtrar para referencia
if (exists("exprs_matrix_unfiltered")) {
    save(exprs_matrix_unfiltered, file = file.path(output_dir, "GSE285498_exprs_matrix_unfiltered.RData"))
    write.csv(exprs_matrix_unfiltered, 
              file = file.path(output_dir, "GSE285498_exprs_matrix_unfiltered.csv"),
              row.names = TRUE)
    cat("  ✓ Matriz de expresión sin filtrar guardada: GSE285498_exprs_matrix_unfiltered.csv\n")
}

# Guardar objetos de filtrado
if (exists("filtering_summary")) {
    save(filtering_summary, expression_threshold, variance_threshold, 
         min_samples_with_expression, n_genes_before, n_genes_after,
         file = file.path(output_dir, "GSE285498_filtering_objects.RData"))
    cat("  ✓ Objetos de filtrado guardados: GSE285498_filtering_objects.RData\n")
}

# Guardar phenodata
save(pheno_data, file = file.path(output_dir, "GSE285498_phenodata.RData"))
cat("  ✓ Phenodata guardado: GSE285498_phenodata.RData\n")

cat("\n✓ Carga de datos completada\n")
cat("  Todos los archivos guardados en:", output_dir, "\n\n")

# ==============================================================================
# 9. Análisis de Expresión Diferencial (DEGs)
# Nota: El análisis se realiza sobre los datos filtrados
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
# 9.1. Preparar diseño experimental
# ------------------------------------------------------------------------------
cat("9.1. Preparando diseño experimental...\n")

# Extraer información de tratamiento directamente de los nombres de las muestras
# Los nombres son: CBD_1, CBD_2, CBD_3, Eto_1, Eto_2, Eto_3, Merge_1, Merge_2, Merge_3, Mock_1, Mock_2, Mock_3
treatment_info <- rep(NA, length(sample_names))

for (i in 1:length(sample_names)) {
    sample_name <- sample_names[i]
    
    # Extraer la parte antes del guión bajo (o el nombre completo si no hay guión)
    # Patrón: TRATAMIENTO_NUMERO
    if (grepl("_", sample_name)) {
        treatment_part <- strsplit(sample_name, "_")[[1]][1]
    } else {
        treatment_part <- sample_name
    }
    
    # Detectar tratamiento (case-insensitive)
    treatment_part_upper <- toupper(treatment_part)
    
    if (treatment_part_upper == "CBD") {
        treatment_info[i] <- "CBD"
    } else if (treatment_part_upper == "ETO") {
        treatment_info[i] <- "Etoposide"
    } else if (treatment_part_upper == "MERGE") {
        treatment_info[i] <- "Combination"
    } else if (treatment_part_upper == "MOCK") {
        treatment_info[i] <- "Mock"
    } else {
        # Fallback: buscar patrones más flexibles
        if (grepl("^CBD", sample_name, ignore.case = TRUE)) {
            treatment_info[i] <- "CBD"
        } else if (grepl("^Eto", sample_name, ignore.case = TRUE)) {
            treatment_info[i] <- "Etoposide"
        } else if (grepl("^Merge", sample_name, ignore.case = TRUE)) {
            treatment_info[i] <- "Combination"
        } else if (grepl("^Mock|mock|non-treatment", sample_name, ignore.case = TRUE)) {
            treatment_info[i] <- "Mock"
        }
    }
}

# Si aún hay NAs, intentar desde phenodata como respaldo
if (any(is.na(treatment_info))) {
    cat("  Algunas muestras no tienen tratamiento detectado, intentando desde phenodata...\n")
    if (exists("pheno_subset") && "treatment:ch1" %in% colnames(pheno_subset)) {
        for (i in which(is.na(treatment_info))) {
            idx <- match(sample_names[i], rownames(pheno_subset))
            if (!is.na(idx)) {
                treatment_info[i] <- pheno_subset[["treatment:ch1"]][idx]
            }
        }
    } else if ("treatment:ch1" %in% colnames(pheno_data)) {
        for (i in which(is.na(treatment_info))) {
            idx <- which(pheno_data$geo_accession == sample_names[i] | 
                        grepl(sample_names[i], pheno_data$title, ignore.case = TRUE))
            if (length(idx) > 0) {
                treatment_info[i] <- pheno_data[["treatment:ch1"]][idx[1]]
            }
        }
    }
}

# Limpiar y estandarizar nombres de tratamiento
treatment_info <- gsub(".*non-treatment.*", "Mock", treatment_info, ignore.case = TRUE)
treatment_info <- gsub(".*Cannabidiol.*", "CBD", treatment_info, ignore.case = TRUE)
treatment_info <- gsub(".*Etoposide.*", "Etoposide", treatment_info, ignore.case = TRUE)
treatment_info <- gsub(".*Combination.*|.*Merge.*", "Combination", treatment_info, ignore.case = TRUE)

# Verificar que todas las muestras tienen tratamiento asignado
if (any(is.na(treatment_info))) {
    cat("  ⚠ ADVERTENCIA: Algunas muestras no tienen tratamiento asignado:\n")
    print(sample_names[is.na(treatment_info)])
    cat("  Asignando 'Unknown' a estas muestras\n")
    treatment_info[is.na(treatment_info)] <- "Unknown"
}

# Crear factor de tratamiento
treatment <- factor(treatment_info, levels = c("Mock", "CBD", "Etoposide", "Combination"))
treatment <- droplevels(treatment)  # Eliminar niveles no usados

cat("  Distribución de tratamientos:\n")
print(table(treatment))
cat("\n")

# Crear data frame de información de muestras
sample_info <- data.frame(
    Sample = sample_names,
    Treatment = treatment,
    stringsAsFactors = FALSE
)

# Agregar GSM IDs si están disponibles
if (exists("pheno_subset") && "geo_accession" %in% colnames(pheno_subset)) {
    sample_info$GSM_ID <- pheno_subset$geo_accession[match(sample_names, rownames(pheno_subset))]
} else if ("geo_accession" %in% colnames(pheno_data)) {
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

# Verificar que hay réplicas suficientes
replicates_table <- table(treatment)
min_replicates <- min(replicates_table)
cat("  Réplicas por tratamiento:\n")
print(replicates_table)
cat("  Mínimo de réplicas:", min_replicates, "\n")
if (min_replicates < 2) {
    cat("  ⚠ ADVERTENCIA: Algunos tratamientos tienen menos de 2 réplicas\n")
    cat("  El análisis puede no ser confiable\n")
}
cat("\n")

# Crear matriz de diseño
# Usar Mock como referencia (intercepto)
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

# ------------------------------------------------------------------------------
# 9.2. Análisis de expresión diferencial con limma
# ------------------------------------------------------------------------------
cat("9.2. Realizando análisis de expresión diferencial con limma...\n")
cat("  Nota: Usando limma con datos FPKM (puede tener limitaciones)\n")
cat("  Para mejores resultados con RNA-seq, se recomienda usar conteos crudos con edgeR/DESeq2\n\n")

# Transformar a log2 si no está ya transformado
# Asumir que los datos pueden estar en escala lineal (FPKM) o log2
if (max(exprs_matrix, na.rm = TRUE) > 100) {
    cat("  Detectando datos en escala lineal, transformando a log2...\n")
    log_exprs_degs <- log2(exprs_matrix + 1)
} else {
    cat("  Datos parecen estar ya en escala log2, usando directamente...\n")
    log_exprs_degs <- exprs_matrix
}

# Ajustar modelo lineal
cat("  Ajustando modelo lineal...\n")
fit <- lmFit(log_exprs_degs, design)
cat("  ✓ Modelo ajustado\n\n")

# Aplicar eBayes para suavizar varianzas
cat("  Aplicando eBayes para suavizar varianzas...\n")
fit2 <- eBayes(fit, trend = TRUE, robust = TRUE)
cat("  ✓ eBayes aplicado\n\n")

# Determinar qué comparaciones hacer
# Comparar cada tratamiento vs Mock (control)
comparisons <- colnames(design)[grepl("treatment", colnames(design), ignore.case = TRUE)]
if (length(comparisons) == 0) {
    # Si no hay columnas de tratamiento, usar todas excepto intercepto
    comparisons <- colnames(design)[-1]
}

cat("  Comparaciones a realizar:\n")
for (coef in comparisons) {
    treatment_name <- gsub("treatment", "", coef, ignore.case = TRUE)
    treatment_name <- gsub("\\.", "", treatment_name)
    cat("    -", treatment_name, "vs Mock (coeficiente:", coef, ")\n")
}
cat("\n")

# Lista para almacenar todos los resultados
all_results <- list()
all_degs <- list()

# Realizar análisis para cada comparación
for (coef in comparisons) {
    treatment_name <- gsub("treatment", "", coef, ignore.case = TRUE)
    treatment_name <- gsub("\\.", "", treatment_name)
    if (treatment_name == "") treatment_name <- "Treatment"
    
    cat("  Analizando:", treatment_name, "vs Mock...\n")
    
    # Extraer resultados
    results <- topTable(fit2, coef = coef, number = Inf, sort.by = "P")
    
    # Verificar columnas disponibles
    cat("    Columnas disponibles en results:", paste(colnames(results), collapse = ", "), "\n")
    
    # Agregar símbolos de genes (usar rownames)
    results$GeneID <- rownames(results)
    
    # Intentar anotar con símbolos de genes si es posible
    if (require("org.Hs.eg.db", quietly = TRUE)) {
        tryCatch({
            # Intentar mapear IDs a símbolos
            gene_symbols <- mapIds(org.Hs.eg.db,
                                  keys = results$GeneID,
                                  column = "SYMBOL",
                                  keytype = "ENSEMBL",
                                  multiVals = "first")
            # Si no funciona con ENSEMBL, intentar directamente
            if (all(is.na(gene_symbols))) {
                gene_symbols <- mapIds(org.Hs.eg.db,
                                      keys = results$GeneID,
                                      column = "SYMBOL",
                                      keytype = "ENTREZID",
                                      multiVals = "first")
            }
            # Si aún no funciona, usar el ID como símbolo
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
    
    # Reordenar columnas (solo las que existen)
    desired_cols <- c("GeneID", "SYMBOL", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
    available_cols <- intersect(desired_cols, colnames(results))
    missing_cols <- setdiff(desired_cols, colnames(results))
    
    if (length(missing_cols) > 0) {
        cat("    ⚠ ADVERTENCIA: Columnas faltantes:", paste(missing_cols, collapse = ", "), "\n")
    }
    
    # Reordenar solo con las columnas disponibles
    if (length(available_cols) > 0) {
        # Mantener todas las columnas, pero poner las deseadas primero
        other_cols <- setdiff(colnames(results), available_cols)
        results <- results[, c(available_cols, other_cols)]
    }
    
    # Guardar resultados
    all_results[[treatment_name]] <- results
    
    # Definir umbrales para DEGs
    fc_threshold <- 0.5  # log2 fold change (aproximadamente 1.41x)
    pval_threshold <- 0.05  # FDR
    
    # Identificar DEGs
    degs <- results[
        abs(results$logFC) >= fc_threshold & 
        results$adj.P.Val < pval_threshold,
    ]
    
    # Separar up y down regulated
    up_genes <- degs[degs$logFC > 0, ]
    down_genes <- degs[degs$logFC < 0, ]
    
    cat("    Total DEGs (|logFC| >=", fc_threshold, "y adj.P.Val <", pval_threshold, "):", nrow(degs), "\n")
    cat("    Up-regulated (", treatment_name, " > Mock):", nrow(up_genes), "\n")
    cat("    Down-regulated (", treatment_name, " < Mock):", nrow(down_genes), "\n\n")
    
    all_degs[[treatment_name]] <- list(
        all = degs,
        up = up_genes,
        down = down_genes,
        results = results
    )
}

# ------------------------------------------------------------------------------
# 9.3. Generar gráficos para cada comparación
# ------------------------------------------------------------------------------
cat("9.3. Generando gráficos de expresión diferencial...\n")

for (treatment_name in names(all_degs)) {
    cat("  Generando gráficos para:", treatment_name, "vs Mock...\n")
    
    results <- all_degs[[treatment_name]]$results
    degs <- all_degs[[treatment_name]]$all
    up_genes <- all_degs[[treatment_name]]$up
    down_genes <- all_degs[[treatment_name]]$down
    
    # Crear subdirectorio para esta comparación
    comp_dir <- file.path(degs_dir, treatment_name)
    if (!dir.exists(comp_dir)) {
        dir.create(comp_dir, recursive = TRUE)
    }
    
    fc_threshold <- 0.5
    pval_threshold <- 0.05
    
    # Volcano Plot
    png(file.path(comp_dir, "volcano_plot.png"), 
        width = 12, height = 10, units = "in", res = 300)
    
    # Usar SYMBOL si existe y tiene el mismo número de filas, sino usar GeneID
    if ("SYMBOL" %in% colnames(results) && length(results$SYMBOL) == nrow(results)) {
        gene_labels <- ifelse(is.na(results$SYMBOL) | results$SYMBOL == "", 
                              results$GeneID, results$SYMBOL)
    } else {
        gene_labels <- results$GeneID
    }
    
    plot_data <- data.frame(
        logFC = results$logFC,
        neg_log10_pval = -log10(results$adj.P.Val + 1e-300),  # Evitar log(0)
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
            title = paste("Volcano Plot:", treatment_name, "vs Mock"),
            subtitle = paste("DEGs:", nrow(degs), "| Up:", nrow(up_genes), "| Down:", nrow(down_genes)),
            x = paste("Log2 Fold Change (", treatment_name, " vs Mock)"),
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
    cat("    ✓ Volcano plot guardado\n")
    
    # MA-plot
    png(file.path(comp_dir, "MA_plot.png"), 
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
            title = paste("MA-plot:", treatment_name, "vs Mock"),
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
    cat("    ✓ MA-plot guardado\n")
    
    
# ------------------------------------------------------------------------------
# 9.4. Guardar tablas de resultados
# ------------------------------------------------------------------------------
cat("\n9.4. Guardando tablas de resultados...\n")

for (treatment_name in names(all_degs)) {
    comp_dir <- file.path(degs_dir, treatment_name)
    
    # Guardar todos los resultados
    write.table(all_degs[[treatment_name]]$results,
                file = file.path(comp_dir, "all_results.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Guardar todos los DEGs
    if (nrow(all_degs[[treatment_name]]$all) > 0) {
        write.table(all_degs[[treatment_name]]$all,
                    file = file.path(comp_dir, "all_DEGs.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        
        # Guardar genes up-regulated (usar SYMBOL si existe, sino GeneID)
        up_genes_list <- all_degs[[treatment_name]]$up
        if ("SYMBOL" %in% colnames(up_genes_list) && length(up_genes_list$SYMBOL) == nrow(up_genes_list)) {
            up_gene_names <- ifelse(is.na(up_genes_list$SYMBOL) | up_genes_list$SYMBOL == "", 
                                   up_genes_list$GeneID, up_genes_list$SYMBOL)
        } else {
            up_gene_names <- up_genes_list$GeneID
        }
        write.table(up_gene_names,
                    file = file.path(comp_dir, "up_genes.txt"),
                    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        
        # Guardar genes down-regulated (usar SYMBOL si existe, sino GeneID)
        down_genes_list <- all_degs[[treatment_name]]$down
        if ("SYMBOL" %in% colnames(down_genes_list) && length(down_genes_list$SYMBOL) == nrow(down_genes_list)) {
            down_gene_names <- ifelse(is.na(down_genes_list$SYMBOL) | down_genes_list$SYMBOL == "", 
                                     down_genes_list$GeneID, down_genes_list$SYMBOL)
        } else {
            down_gene_names <- down_genes_list$GeneID
        }
        write.table(down_gene_names,
                    file = file.path(comp_dir, "down_genes.txt"),
                    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        
        # Guardar tablas completas
        write.table(all_degs[[treatment_name]]$up,
                    file = file.path(comp_dir, "up_genes_complete.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        
        write.table(all_degs[[treatment_name]]$down,
                    file = file.path(comp_dir, "down_genes_complete.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
    }
    
    cat("  ✓ Resultados guardados para:", treatment_name, "\n")
}

# Resumen general
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

cat("\n  Resumen de análisis:\n")
print(summary_stats)

cat("\n", rep("=", 70), "\n")
cat("ANÁLISIS DE EXPRESIÓN DIFERENCIAL COMPLETADO\n")
cat(rep("=", 70), "\n\n")

# Guardar objetos de DEGs
save(fit, fit2, all_results, all_degs, design, treatment, sample_info, summary_stats,
     file = file.path(degs_dir, "DEGs_objects.RData"))
cat("  ✓ Objetos de DEGs guardados: DEGs_objects.RData\n\n")

# ==============================================================================
# 10. Generar Reporte HTML de Resumen
# ==============================================================================
cat("\n", rep("=", 70), "\n")
cat("GENERANDO REPORTE HTML DE RESUMEN\n")
cat(rep("=", 70), "\n\n")

cat("10. Creando reporte HTML...\n")

# Crear contenido HTML
html_content <- paste0('<!DOCTYPE html>
<html lang="es">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Resumen del Análisis - GSE285498</title>
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
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }
        h2 {
            color: #34495e;
            margin-top: 30px;
            border-left: 4px solid #3498db;
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
            background-color: #3498db;
            color: white;
            font-weight: bold;
        }
        tr:nth-child(even) {
            background-color: #f2f2f2;
        }
        tr:hover {
            background-color: #e8f4f8;
        }
        .metric-box {
            background-color: #ecf0f1;
            border-left: 4px solid #3498db;
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
            background-color: #3498db;
            color: white;
            text-decoration: none;
            border-radius: 5px;
            margin: 5px;
        }
        .link-button:hover {
            background-color: #2980b9;
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
    </style>
</head>
<body>
    <div class="container">
        <h1>📊 Resumen del Análisis - GSE285498</h1>
        
        <div class="timestamp">
            Generado el: ', Sys.time(), '
        </div>
        
        <h2>📋 Información del Dataset</h2>
        <div class="summary-grid">
            <div class="info-box">
                <div class="metric-label">Dataset</div>
                <div class="metric-value">GSE285498</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Tipo de Datos</div>
                <div class="metric-value">RNA-seq bulk</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Plataforma</div>
                <div class="metric-value">Illumina HiSeq4000</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Cuantificación</div>
                <div class="metric-value">TopHat + Cufflinks (FPKM)</div>
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
                <a href="GSE285498_filtering_summary.csv" class="link-button">Descargar Resumen de Filtrado (CSV)</a>
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
            <td><strong>', summary_stats$Comparison[i], ' vs Mock</strong></td>
            <td>', summary_stats$Total_DEGs[i], '</td>
            <td>', summary_stats$Up_regulated[i], '</td>
            <td>', summary_stats$Down_regulated[i], '</td>
        </tr>')
    }
    html_content <- paste0(html_content, '</table>')
    
    # Agregar secciones por comparación
    for (treatment_name in names(all_degs)) {
        degs <- all_degs[[treatment_name]]$all
        up_genes <- all_degs[[treatment_name]]$up
        down_genes <- all_degs[[treatment_name]]$down
        
        html_content <- paste0(html_content, '
        <div class="comparison-section">
            <h3>🔬 ', treatment_name, ' vs Mock</h3>
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
                <a href="DEGs/', treatment_name, '/volcano_plot.png" class="link-button" target="_blank">Volcano Plot</a>
                <a href="DEGs/', treatment_name, '/MA_plot.png" class="link-button" target="_blank">MA-plot</a>')
        
        if (nrow(degs) > 0) {
            html_content <- paste0(html_content, '
                <a href="DEGs/', treatment_name, '/heatmap_top_genes.png" class="link-button" target="_blank">Heatmap Top Genes</a>')
        }
        
        html_content <- paste0(html_content, '
            </p>
            
            <h4>Archivos de Resultados</h4>
            <p>
                <a href="DEGs/', treatment_name, '/all_results.tsv" class="link-button">Todos los Resultados (TSV)</a>
                <a href="DEGs/', treatment_name, '/all_DEGs.tsv" class="link-button">Todos los DEGs (TSV)</a>')
        
        if (nrow(up_genes) > 0) {
            html_content <- paste0(html_content, '
                <a href="DEGs/', treatment_name, '/up_genes.txt" class="link-button">Genes Up-regulated</a>
                <a href="DEGs/', treatment_name, '/down_genes.txt" class="link-button">Genes Down-regulated</a>')
        }
        
        html_content <- paste0(html_content, '
            </p>')
        
        # Mostrar top 10 DEGs
        if (nrow(degs) > 0) {
            top_degs <- degs[order(degs$adj.P.Val), ][1:min(10, nrow(degs)), ]
            html_content <- paste0(html_content, '
            <h4>Top 10 DEGs (por FDR)</h4>
            <table>
                <tr>
                    <th>Gen</th>
                    <th>Log2 FC</th>
                    <th>P-value</th>
                    <th>FDR</th>
                </tr>')
            
            for (j in 1:nrow(top_degs)) {
                html_content <- paste0(html_content, '<tr>
                    <td>', ifelse("SYMBOL" %in% colnames(top_degs), top_degs$SYMBOL[j], top_degs$GeneID[j]), '</td>
                    <td>', round(top_degs$logFC[j], 3), '</td>
                    <td>', format(top_degs$P.Value[j], scientific = TRUE, digits = 3), '</td>
                    <td>', format(top_degs$adj.P.Val[j], scientific = TRUE, digits = 3), '</td>
                </tr>')
            }
            html_content <- paste0(html_content, '</table>')
        }
        
        html_content <- paste0(html_content, '</div>')
    }
}

# Cerrar HTML
html_content <- paste0(html_content, '
        <h2>📁 Archivos Generados</h2>
        <div class="metric-box">
            <h3>Estructura de Directorios</h3>
            <ul>
                <li><strong>QC/</strong> - Análisis de calidad y gráficos</li>
                <li><strong>DEGs/</strong> - Resultados de expresión diferencial
                    <ul>
                        <li><strong>CBD/</strong> - Resultados CBD vs Mock</li>
                        <li><strong>Etoposide/</strong> - Resultados Etoposide vs Mock</li>
                        <li><strong>Combination/</strong> - Resultados Combination vs Mock</li>
                    </ul>
                </li>
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
html_file <- file.path(output_dir, "GSE285498_analysis_summary.html")
writeLines(html_content, html_file)
cat("  ✓ Reporte HTML guardado: GSE285498_analysis_summary.html\n")
cat("  Abre el archivo en tu navegador para ver el resumen completo\n\n")

# ==============================================================================
# OBJETOS DISPONIBLES:
# ==============================================================================
# - raw_data: Objeto DGEList con los datos crudos
# - exprs_matrix: Matriz de expresión (genes x muestras)
# - pheno_data: Phenodata completo
# - dge: Objeto DGEList con phenodata asociado (si fue posible)
# - all_results: Lista con todos los resultados de expresión diferencial
# - all_degs: Lista con DEGs identificados para cada comparación
# ==============================================================================

