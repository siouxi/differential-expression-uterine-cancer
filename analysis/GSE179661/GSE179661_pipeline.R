required_pkgs <- c(
  "tidyverse", "here", "jsonlite", "glue", "fs",
  "BiocManager", "GEOquery", "DESeq2", "edgeR",
  "SummarizedExperiment", "apeglm", "EnhancedVolcano",
  "pheatmap", "ComplexHeatmap", "RColorBrewer", "plotly",
  "clusterProfiler", "enrichplot", "msigdbr", "fgsea",
  "org.Hs.eg.db", "AnnotationDbi", "biomaRt", "GSVA",
  "corrplot", "vsn",
  "viridis", "VennDiagram", "ggrepel"
)

install_and_load <- function(pkgs) {
  # Verificar qué paquetes ya están instalados (chequeo rápido)
  installed <- rownames(installed.packages())
  new_pkgs <- pkgs[!pkgs %in% installed]
  
  if (length(new_pkgs) > 0) {
    cat("⚠️  Instalando paquetes faltantes (esto puede tardar un poco):", paste(new_pkgs, collapse = ", "), "\n")
    BiocManager::install(new_pkgs, ask = FALSE, update = FALSE)
    cat("✅ Instalación de paquetes completada\n\n")
  } else {
    cat("✅ Todos los paquetes requeridos ya están instalados\n\n")
  }
  
  # Cargar paquetes con mensajes de progreso para los grandes
  cat("Cargando paquetes...\n")
  for (pkg in pkgs) {
    if (pkg %in% c("org.Hs.eg.db", "clusterProfiler", "biomaRt", "msigdbr")) {
      cat("  Cargando", pkg, "(esto puede tardar un momento)...\n")
    }
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE, quietly = TRUE)
    )
  }
  cat("✅ Todos los paquetes cargados\n\n")
}

install_and_load(required_pkgs)


geo_id <- "GSE179661"

# Prefijo común para títulos de gráficos: ej. "GSE179661 / <anotación(gset)>"
plot_title_prefix <- tryCatch({
  if (exists("gset")) {
    ann <- tryCatch(as.character(annotation(gset)), error = function(e) "")
    if (!is.null(ann) && ann != "") {
      paste(geo_id, "/", ann)
    } else {
      geo_id
    }
  } else {
    geo_id
  }
}, error = function(e) geo_id)

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


base_dir <- dirname(dirname(script_dir))

paths <- list(
  data_dir = file.path(base_dir, "DATA", geo_id),
  results_dir = file.path(base_dir, "results", geo_id),
  qc_dir = file.path(base_dir, "results", geo_id, "qc"),
  deseq_dir = file.path(base_dir, "results", geo_id, "deseq2"),
  enrich_dir = file.path(base_dir, "results", geo_id, "enrichment")
)

fs::dir_create(unlist(paths))

cat("Rutas del proyecto configuradas:\n")
cat("  Directorio base:", base_dir, "\n")
cat("  Directorio de datos:", paths$data_dir, "\n")
cat("  Directorio de resultados:", paths$results_dir, "\n\n")


# 3. Cargar y procesar datos de fenotipo ----------------------------------------------------------------
cat("Cargando datos de fenotipo...\n")
pheno_file <- file.path(paths$data_dir, paste0(geo_id, "_pheno_data.csv"))

if (!file.exists(pheno_file)) {
  stop(paste("Archivo de fenotipo no encontrado:", pheno_file))
}

pheno_raw <- readr::read_csv(pheno_file, show_col_types = FALSE)

# Extraer y procesar metadatos de muestras
# Del pheno_data, tenemos: title, geo_accession, treatment:ch1, source_name_ch1
sample_sheet <- pheno_raw %>%
  dplyr::select(
    sample_id = title,
    geo_accession,
    treatment = `treatment:ch1`,
    cell_line = source_name_ch1
  ) %>%
  dplyr::mutate(
    # Extraer condición de la columna treatment
    condition = dplyr::case_when(
      grepl("none treated", treatment, ignore.case = TRUE) ~ "Control",
      grepl("CBD", treatment, ignore.case = TRUE) ~ "CBD",
      TRUE ~ "Unknown"
    ),
    # Extraer tipo de línea celular
    cell_line = dplyr::case_when(
      grepl("MHCC97H", cell_line, ignore.case = TRUE) ~ "MHCC97H",
      grepl("HepG2", cell_line, ignore.case = TRUE) ~ "HepG2",
      TRUE ~ cell_line
    ),
    # Extraer número de réplica
    replicate = stringr::str_extract(sample_id, "rep\\d+"),
    replicate = as.numeric(stringr::str_replace(replicate, "rep", "")),
    # Crear IDs de muestra que coincidan con los nombres de columna de la matriz de conteo
    # La matriz de conteo probablemente usa IDs GSM o nombres de muestra del título
    sample_id_clean = geo_accession
  ) %>%
  dplyr::mutate(
    condition = factor(condition, levels = c("Control", "CBD")),
    cell_line = factor(cell_line)
  )

cat("  Número de muestras:", nrow(sample_sheet), "\n")
cat("  Condiciones:", paste(levels(sample_sheet$condition), collapse = ", "), "\n")
cat("  Líneas celulares:", paste(levels(sample_sheet$cell_line), collapse = ", "), "\n\n")


# 4. Cargar matriz de conteos de genes crudos ----------------------------------------------------------------
cat("Cargando matriz de conteos de genes crudos...\n")
counts_file <- file.path(paths$data_dir, paste0(geo_id, "_Raw_gene_counts_matrix.txt.gz"))

if (!file.exists(counts_file)) {
  stop(paste("Archivo de conteos no encontrado:", counts_file))
}

# Leer archivo de conteos comprimido
# El formato del archivo puede variar - intentar diferentes enfoques
counts_df <- tryCatch({
  readr::read_tsv(counts_file, comment = "#", show_col_types = FALSE)
}, error = function(e) {
  # Alternativa: leer como tabla
  read.table(gzfile(counts_file), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
})

# Identificar columna de ID de gen (podría ser Geneid, gene_id, o la primera columna)
gene_col <- NULL
if ("gene_id" %in% colnames(counts_df)) {
  gene_col <- "gene_id"
} else if ("gene_id" %in% colnames(counts_df)) {
  gene_col <- "gene_id"
} else {
  gene_col <- colnames(counts_df)[1]
}

# Preparar matriz de conteo
# Eliminar columnas de anotación si están presentes
annotation_cols <- c("Chr", "Start", "End", "Length", "Strand", "Geneid", "gene_id")
cols_to_remove <- intersect(annotation_cols, colnames(counts_df))
cols_to_remove <- setdiff(cols_to_remove, gene_col)

if (length(cols_to_remove) > 0) {
  counts_mat <- counts_df %>%
    dplyr::select(-dplyr::all_of(cols_to_remove)) %>%
    tibble::column_to_rownames(gene_col) %>%
    as.matrix()
} else {
  counts_mat <- counts_df %>%
    tibble::column_to_rownames(gene_col) %>%
    as.matrix()
}

# Asegurar que los conteos sean numéricos
counts_mat <- apply(counts_mat, 2, as.numeric)
rownames(counts_mat) <- rownames(counts_df)

cat("  Dimensiones de la matriz de conteo:", dim(counts_mat), "\n")
cat("  Nombres de muestras en la matriz:", paste(head(colnames(counts_mat), 5), collapse = ", "), "...\n\n")

# Emparejar nombres de muestras entre matriz de conteo y datos de fenotipo
cat("Emparejando nombres de muestras...\n")
sample_cols <- colnames(counts_mat)

# Diagnóstico: mostrar lo que tenemos
cat("  Columnas de matriz de conteo (", length(sample_cols), "):", paste(sample_cols, collapse = ", "), "\n")
cat("  Sample_id de fenotipo (", length(sample_sheet$sample_id), "):", paste(sample_sheet$sample_id, collapse = ", "), "\n")
cat("  Geo_accession de fenotipo (", length(sample_sheet$geo_accession), "):", paste(sample_sheet$geo_accession, collapse = ", "), "\n")

# Función auxiliar para normalizar nombres de muestras
# Convierte "H_0 rep1" -> "h_0_1", "HE_40 rep2" -> "he_40_2", etc.
normalize_sample_name <- function(name) {
  # Convertir a minúsculas
  name_lower <- tolower(name)
  # Reemplazar " rep" o "rep" con "_"
  name_lower <- gsub("\\s*rep\\s*", "_", name_lower)
  # Eliminar espacios extra y normalizar guiones bajos
  name_lower <- gsub("\\s+", "_", name_lower)
  name_lower <- gsub("_+", "_", name_lower)
  # Eliminar guiones bajos iniciales/finales
  name_lower <- gsub("^_|_$", "", name_lower)
  return(name_lower)
}

# Crear versiones normalizadas para emparejamiento
sample_sheet <- sample_sheet %>%
  dplyr::mutate(
    sample_id_normalized = normalize_sample_name(sample_id),
    geo_accession_normalized = tolower(geo_accession)
  )

sample_cols_normalized <- tolower(sample_cols)

# Estrategia 1: Coincidencia directa con sample_id normalizado
cat("  Estrategia 1: Emparejamiento por sample_id normalizado...\n")
matched_indices <- match(sample_cols_normalized, sample_sheet$sample_id_normalized)
matched_samples <- sample_cols[!is.na(matched_indices)]

if (length(matched_samples) > 0) {
  cat("    ✓ Emparejadas", length(matched_samples), "muestras por sample_id normalizado\n")
  sample_sheet$sample_id_clean[matched_indices[!is.na(matched_indices)]] <- sample_cols[!is.na(matched_indices)]
} else {
  # Estrategia 2: Coincidencia directa con geo_accession
  cat("  Estrategia 2: Emparejamiento por geo_accession...\n")
  matched_indices <- match(sample_cols_normalized, sample_sheet$geo_accession_normalized)
  matched_samples <- sample_cols[!is.na(matched_indices)]
  
  if (length(matched_samples) > 0) {
    cat("    ✓ Emparejadas", length(matched_samples), "muestras por geo_accession\n")
    sample_sheet$sample_id_clean[matched_indices[!is.na(matched_indices)]] <- sample_cols[!is.na(matched_indices)]
  } else {
    # Estrategia 3: Coincidencia parcial - verificar si los nombres de la matriz contienen partes de sample_id
    cat("  Estrategia 3: Emparejamiento parcial por patrón...\n")
    matched_samples <- character(0)
    matched_pheno_indices <- integer(0)
    
    for (i in seq_along(sample_cols)) {
      col_name <- sample_cols_normalized[i]
      # Intentar coincidir con sample_id normalizado
      for (j in seq_along(sample_sheet$sample_id_normalized)) {
        pheno_name <- sample_sheet$sample_id_normalized[j]
        # Verificar si uno contiene al otro (emparejamiento flexible)
        if (grepl(gsub("_", ".", col_name), pheno_name, ignore.case = TRUE) ||
            grepl(gsub("_", ".", pheno_name), col_name, ignore.case = TRUE) ||
            # Verificar patrón como "h_0" coincidiendo con "H_0 rep1"
            grepl(paste0("^", gsub("_\\d+$", "", col_name)), pheno_name, ignore.case = TRUE)) {
          if (!(j %in% matched_pheno_indices)) {
            matched_samples <- c(matched_samples, sample_cols[i])
            matched_pheno_indices <- c(matched_pheno_indices, j)
            sample_sheet$sample_id_clean[j] <- sample_cols[i]
            break
          }
        }
      }
    }
    
    if (length(matched_samples) > 0) {
      cat("    ✓ Emparejadas", length(matched_samples), "muestras por patrón parcial\n")
    } else {
      # Estrategia 4: Intentar emparejar extrayendo patrones clave (H_0, H_40, HE_0, HE_40)
      cat("  Estrategia 4: Emparejamiento por patrones extraídos...\n")
      matched_samples <- character(0)
      matched_pheno_indices <- integer(0)
      
      # Extraer patrón de nombres de matriz de conteo (ej: "h_0_1" -> "h_0")
      extract_pattern <- function(name) {
        # Eliminar números finales y guiones bajos
        pattern <- gsub("_\\d+$", "", tolower(name))
        return(pattern)
      }
      
      # Extraer patrón de nombres de fenotipo (ej: "H_0 rep1" -> "h_0")
      extract_pheno_pattern <- function(name) {
        pattern <- tolower(name)
        pattern <- gsub("\\s*rep\\d+", "", pattern)
        pattern <- gsub("\\s+", "_", pattern)
        return(pattern)
      }
      
      count_patterns <- sapply(sample_cols, extract_pattern)
      pheno_patterns <- sapply(sample_sheet$sample_id, extract_pheno_pattern)
      
      for (i in seq_along(count_patterns)) {
        match_idx <- which(pheno_patterns == count_patterns[i])
        if (length(match_idx) > 0 && !(match_idx[1] %in% matched_pheno_indices)) {
          matched_samples <- c(matched_samples, sample_cols[i])
          matched_pheno_indices <- c(matched_pheno_indices, match_idx[1])
          sample_sheet$sample_id_clean[match_idx[1]] <- sample_cols[i]
        }
      }
      
      if (length(matched_samples) > 0) {
        cat("    ✓ Emparejadas", length(matched_samples), "muestras por patrones extraídos\n")
      }
    }
  }
}

# Chequeo final: si aún no hay coincidencias, usar asignación secuencial como respaldo
if (length(matched_samples) == 0 || length(matched_samples) < nrow(sample_sheet)) {
  cat("  Estrategia 5: Respaldo - asignación secuencial\n")
  unmatched_pheno <- which(is.na(sample_sheet$sample_id_clean) | sample_sheet$sample_id_clean == sample_sheet$geo_accession)
  unmatched_counts <- setdiff(sample_cols, matched_samples)
  
  if (length(unmatched_pheno) > 0 && length(unmatched_counts) > 0) {
    n_assign <- min(length(unmatched_pheno), length(unmatched_counts))
    for (i in 1:n_assign) {
      sample_sheet$sample_id_clean[unmatched_pheno[i]] <- unmatched_counts[i]
      matched_samples <- c(matched_samples, unmatched_counts[i])
    }
    cat("    ⚠ Asignadas", n_assign, "muestras secuencialmente (¡verificar manualmente!)\n")
  }
}

# Asegurar que todas las muestras tengan un ID limpio
if (any(is.na(sample_sheet$sample_id_clean))) {
  sample_sheet$sample_id_clean[is.na(sample_sheet$sample_id_clean)] <- sample_sheet$geo_accession[is.na(sample_sheet$sample_id_clean)]
}

# Mostrar resultados de emparejamiento
cat("\n  Resultados de emparejamiento:\n")
matching_df <- data.frame(
  Phenotype_ID = sample_sheet$sample_id,
  Count_Matrix_ID = sample_sheet$sample_id_clean,
  Condition = sample_sheet$condition,
  Cell_Line = sample_sheet$cell_line,
  Matched = sample_sheet$sample_id_clean %in% sample_cols
)
print(matching_df)

# Filtrar y reordenar
matched_samples <- unique(matched_samples[matched_samples %in% sample_cols])
sample_sheet <- sample_sheet %>%
  dplyr::filter(sample_id_clean %in% matched_samples | sample_id_clean %in% sample_cols) %>%
  dplyr::arrange(match(sample_id_clean, sample_cols))

# Asegurar que la matriz de conteos tenga las muestras emparejadas en el orden correcto
available_samples <- intersect(sample_sheet$sample_id_clean, sample_cols)
counts_mat <- counts_mat[, available_samples, drop = FALSE]
sample_sheet <- sample_sheet %>%
  dplyr::filter(sample_id_clean %in% available_samples) %>%
  dplyr::arrange(match(sample_id_clean, colnames(counts_mat)))

# Verificación final
if (ncol(counts_mat) != nrow(sample_sheet)) {
  warning("Desajuste: la matriz de conteo tiene ", ncol(counts_mat), " columnas pero sample_sheet tiene ", nrow(sample_sheet), " filas")
}

cat("\n  ✓ Muestras emparejadas finales:", length(available_samples), "\n")
cat("  ✓ Columnas de matriz de conteo:", ncol(counts_mat), "\n")
cat("  ✓ Filas de sample_sheet:", nrow(sample_sheet), "\n\n")


# 5. Crear conjunto de datos DESeq2 y pre-filtrado ------------------------------------------------
cat("Creando conjunto de datos DESeq2...\n")

# Crear metadatos de muestras para DESeq2
col_data <- sample_sheet %>%
  dplyr::select(sample_id = sample_id_clean, condition, cell_line, replicate) %>%
  as.data.frame()
rownames(col_data) <- col_data$sample_id

# Crear conjunto de datos DESeq2
# Diseño: tener en cuenta línea celular y condición
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = col_data,
  design = ~ cell_line + condition
)

# Pre-filtrado: mantener genes con suficiente expresión
# IMPORTANTE PARA ESTUDIOS DE CÁNCER: Usar filtrado menos agresivo para preservar
# genes sub-expresados que pueden ser biológicamente relevantes (ej. supresores de tumores)
# Estrategia: al menos 5 conteos en al menos 25% de las muestras (menos estricto que lo típico)
# Esto preserva genes que podrían estar regulados a la baja en cáncer mientras elimina
# genes con esencialmente ninguna señal (ceros técnicos)
n_samples <- ncol(dds)
min_samples_with_counts <- max(2, ceiling(n_samples * 0.25))  # Al menos 25% o mínimo 2 (menos estricto)
count_threshold <- 5  # Umbral más bajo (5 en lugar de 10) para preservar genes sub-expresados

keep <- rowSums(DESeq2::counts(dds) >= count_threshold) >= min_samples_with_counts
dds <- dds[keep, ]

cat("  Criterios de pre-filtrado (ajustados para estudio de cáncer - preserva genes sub-expresados):\n")
cat("    - Umbral de conteo:", count_threshold, "conteos por muestra (bajado de 10 para preservar genes regulados a la baja)\n")
cat("    - Muestras mínimas:", min_samples_with_counts, "de", n_samples, "muestras (", 
    round(min_samples_with_counts/n_samples*100, 1), "% - menos estricto para capturar genes sub-expresados)\n")
cat("    - Razón: En cáncer, los genes sub-expresados (ej. supresores de tumores) son biológicamente relevantes\n")
cat("  Genes después del filtrado:", nrow(dds), "\n")
cat("  Muestras:", ncol(dds), "\n\n")


# 6. Control de calidad y normalización --------------------------------------------------------
cat("Realizando control de calidad y normalización...\n")
cat("  IMPORTANTE: Se calculan múltiples normalizaciones, pero solo los factores de tamaño de DESeq2\n")
cat("  se usan en el análisis de expresión diferencial. Las otras son solo para QC.\n\n")

# Factores de tamaño (DESeq2)
# Estimaciones: factores de normalización para cada muestra (tiene en cuenta diferencias en tamaño de biblioteca)
# 
# MEDIDAS ESTADÍSTICAS QUE USA:
# - Median ratio method: Calcula la mediana de las razones (ratios) entre cada muestra y una 
#   muestra de referencia (pseudo-reference) construida como la mediana geométrica de todos los genes
# - Fórmula: sizeFactor_i = median( counts_ij / geometric_mean(counts_j) ) para todos los genes j
# - Excluye genes con conteos extremos (outliers) usando el método de mediana absoluta (MAD)
# - Resultado: Factor de normalización por muestra (típicamente entre 0.5 y 2.0)
# 
# Acceder con: sizeFactors(dds) o colData(dds)$sizeFactor
dds <- DESeq2::estimateSizeFactors(dds)
cat("  Factores de tamaño DESeq2:\n")
print(sizeFactors(dds))
cat("\n")

# Definir paleta de colores personalizada para QC (orden de prioridad)
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

sample_colors_qc <- rep(qc_colors, length.out = ncol(dds))
names(sample_colors_qc) <- colnames(dds)

# Estimación de dispersión (DESeq2)
# Estimaciones: dispersión por gen, dispersión ajustada y dispersión a priori
# 
# MEDIDAS ESTADÍSTICAS QUE USA:
# 1. Gene-wise dispersion (dispGeneEst):
#    - Método de momentos: Var = μ + α*μ² donde α es la dispersión
#    - Calcula dispersión inicial para cada gen usando la fórmula: α = (Var - μ) / μ²
#    - Usa estimación de máxima verosimilitud (MLE) para genes con suficientes conteos
# 
# 2. Fitted dispersion (dispFit):
#    - Regresión local (loess): Ajusta una curva suave de dispersión vs. expresión media
#    - Modelo: log(α) ~ f(log(μ)) donde μ es la expresión media del gen
#    - Predice dispersión esperada basada en el nivel de expresión
# 
# 3. Prior dispersion (dispersion):
#    - Shrinkage hacia el valor ajustado: Combina dispersión gen-específica y ajustada
#    - Fórmula: α_final = w*α_gene + (1-w)*α_fitted, donde w es un peso basado en la precisión
#    - Usa distribución a priori (prior) para estabilizar estimaciones de genes con pocos datos
# 
# Acceder con: dispersions(dds), mcols(dds)$dispGeneEst, mcols(dds)$dispFit, mcols(dds)$dispersion
dds <- DESeq2::estimateDispersions(dds)
cat("  Estimación de dispersión completada\n")
disp_range <- range(dispersions(dds)[is.finite(dispersions(dds))], na.rm = TRUE)
if (all(is.finite(disp_range))) {
  cat("  Rango de dispersión por gen:", round(disp_range, 4), "\n\n")
} else {
  cat("  Rango de dispersión por gen: No se puede calcular (valores no finitos)\n\n")
}

# Transformación de estabilización de varianza (VST)
# Transforma conteos a escala tipo log2 con varianza estabilizada
# 
# MEDIDAS ESTADÍSTICAS QUE USA:
# - Transformación basada en la dispersión estimada: Usa los valores de dispersión (α) calculados
#   previamente para estabilizar la varianza
# - Fórmula aproximada: VST(x) = ∫[1/sqrt(μ + α*μ²)] dμ desde 0 hasta x
#   donde μ es la expresión media y α es la dispersión
# - Para valores grandes: Se aproxima a log2(x) (comportamiento logarítmico)
# - Para valores pequeños: Mantiene la escala lineal (evita problemas con ceros)
# - Resultado: Datos transformados donde la varianza es aproximadamente constante (homocedasticidad)
# - Ventaja sobre log2: No necesita pseudocounts y la varianza es más estable
# 
# Parámetro 'blind = FALSE': Usa la información de las condiciones experimentales en la 
# transformación (mejor para análisis exploratorio cuando ya conocemos los grupos)
# 
# Acceder con: assay(vsd) - devuelve conteos transformados
vsd <- DESeq2::vst(dds, blind = FALSE)
cat("  Transformación VST completada\n")
cat("  Dimensiones de matriz VST:", dim(assay(vsd)), "\n\n")

# Normalización TMM (edgeR) para comparación
# Estimaciones: factores de normalización usando método Trimmed Mean of M-values
# 
# MEDIDAS ESTADÍSTICAS QUE USA:
# 1. M-values (log fold changes): M = log2(counts_sample / counts_reference)
#    - Calcula el log2 de la razón entre cada muestra y una muestra de referencia
# 
# 2. A-values (mean expression): A = (log2(counts_sample) + log2(counts_reference)) / 2
#    - Calcula el promedio de expresión en escala logarítmica
# 
# 3. Trimmed Mean (media recortada):
#    - Elimina el 30% superior e inferior de los M-values (por defecto trim=0.3)
#    - Calcula la media de los M-values restantes (40% central)
#    - Esto elimina genes diferencialmente expresados y genes con baja expresión
# 
# 4. Normalization factor: 
#    - Factor = 2^(trimmed_mean_M)
#    - Ajusta el tamaño de biblioteca de cada muestra relativo a la referencia
# 
# 5. Reference sample:
#    - Usa la muestra con el percentil 75 más cercano a la mediana de los conteos totales
#    - O puede usar una muestra específica como referencia
# 
# Ventaja: Robusto a genes diferencialmente expresados (al eliminar extremos)
# 
# Acceder con: dds_edgeR$samples$norm.factors
dds_edgeR <- edgeR::DGEList(counts = counts_mat[rownames(dds), colnames(dds)])
dds_edgeR <- edgeR::calcNormFactors(dds_edgeR, method = "TMM")
cat("  Factores de normalización TMM de edgeR:\n")
print(dds_edgeR$samples$norm.factors)
cat("\n")

# Calcular CPM para QC
# Calcula: Counts Per Million (normalizado por tamaño de biblioteca y factores TMM)
# 
# MEDIDAS ESTADÍSTICAS QUE USA:
# - Fórmula: CPM = (counts / library_size) * 1,000,000
#   donde library_size = suma de conteos de todos los genes en la muestra
# 
# - Con factores TMM: CPM = (counts / (library_size * TMM_factor)) * 1,000,000
#   - Ajusta por diferencias en tamaño de biblioteca Y por sesgos de composición (TMM)
# 
# - Interpretación: Número de conteos por millón de reads mapeados
#   - Permite comparar expresión entre muestras con diferentes tamaños de biblioteca
#   - Ejemplo: 10 CPM significa 10 conteos por cada millón de reads totales
# 
# - Uso: Control de calidad, visualización, filtrado (ej: genes con CPM < 1 en todas las muestras)
# 
# Parámetro 'log = FALSE': Devuelve CPM en escala lineal (no logarítmica)
# 
# Acceder con: cpm_mat - matriz de valores CPM
cpm_mat <- edgeR::cpm(dds_edgeR, log = FALSE)
cat("  Matriz CPM calculada\n")
cat("  Dimensiones de matriz CPM:", dim(cpm_mat), "\n")
cpm_range <- range(cpm_mat[is.finite(cpm_mat)], na.rm = TRUE)
if (all(is.finite(cpm_range))) {
  cat("  Rango CPM:", round(cpm_range, 2), "\n\n")
} else {
  cat("  Rango CPM: No se puede calcular (valores no finitos)\n\n")
}

# Guardar resumen de normalización para inspección
normalization_summary <- data.frame(
  Sample = colnames(dds),
  DESeq2_SizeFactor = sizeFactors(dds),
  edgeR_TMM_Factor = dds_edgeR$samples$norm.factors,
  Library_Size = colSums(counts(dds)),
  Mean_CPM = colMeans(cpm_mat)
)
readr::write_tsv(normalization_summary, 
                 file.path(paths$qc_dir, "normalization_factors_summary.tsv"))
cat("  Resumen de normalización guardado en:", 
    file.path(paths$qc_dir, "normalization_factors_summary.tsv"), "\n\n")

# QC plots directory
fs::dir_create(paths$qc_dir)

# 6.1. Sample count distribution (boxplot) with colors
# Prepare data for boxplot with colors based on condition and cell line
# Ensure counts_mat columns match col_data rows
counts_log2 <- log2(counts_mat[, rownames(col_data), drop = FALSE] + 1)

pdf(file.path(paths$qc_dir, "count_distribution_boxplot.pdf"), width = 14, height = 7)
par(mar = c(10, 4, 4, 6))  # Increase margins for sample names and legend
boxplot(counts_log2, 
        main = "",
        xlab = "", 
        ylab = "Log2(counts + 1)",
        las = 2, 
        cex.axis = 0.65,
        col = sample_colors_qc[colnames(counts_log2)],
        border = "black",
        names = colnames(counts_log2))

# Add legend
legend("topright", 
       legend = c("Control - MHCC97H", "Control - HepG2", "CBD - MHCC97H", "CBD - HepG2"),
       fill = c("#6BAED6", "#2E86AB", "#D67BA8", "#A23B72"),
       cex = 0.8,
       title = "Condition - Cell Line")

# Add grid for better readability
grid(nx = NA, ny = NULL, col = "gray90", lty = "dotted")

dev.off()

# Calculate and report mean values for boxplot data
cat("Boxplot data statistics (log2(counts + 1)):\n")
cat("  ===========================================\n")
boxplot_stats <- data.frame(
  Sample = colnames(counts_log2),
  Mean = colMeans(counts_log2),
  Median = apply(counts_log2, 2, median),
  SD = apply(counts_log2, 2, sd),
  Min = apply(counts_log2, 2, min),
  Max = apply(counts_log2, 2, max),
  Q1 = apply(counts_log2, 2, quantile, probs = 0.25),
  Q3 = apply(counts_log2, 2, quantile, probs = 0.75)
)

# Add condition and cell line info
boxplot_stats <- boxplot_stats %>%
  dplyr::left_join(
    col_data %>% 
      tibble::rownames_to_column("Sample") %>%
      dplyr::select(Sample, Condition = condition, CellLine = cell_line),
    by = "Sample"
  ) %>%
  dplyr::select(Sample, Condition, CellLine, Mean, Median, SD, Min, Max, Q1, Q3)

cat("  Summary statistics by sample:\n")
print(boxplot_stats, row.names = FALSE)
cat("\n")

# Overall statistics
cat("  Overall statistics (all samples combined):\n")
counts_log2_finite <- counts_log2[is.finite(counts_log2)]
if (length(counts_log2_finite) > 0) {
  cat("    Overall mean:", round(mean(counts_log2_finite), 4), "\n")
  cat("    Overall median:", round(median(counts_log2_finite), 4), "\n")
  cat("    Overall SD:", round(sd(counts_log2_finite), 4), "\n")
} else {
  cat("    Warning: No finite values in counts_log2\n")
}
cat("\n")

# Statistics by condition
cat("  Statistics by condition:\n")
for (cond in levels(col_data$condition)) {
  cond_samples <- rownames(col_data)[col_data$condition == cond]
  cond_data <- counts_log2[, cond_samples, drop = FALSE]
  cond_data_finite <- cond_data[is.finite(cond_data)]
  if (length(cond_data_finite) > 0) {
    cat("    ", cond, ":\n")
    cat("      Mean:", round(mean(cond_data_finite), 4), "\n")
    cat("      Median:", round(median(cond_data_finite), 4), "\n")
    cat("      SD:", round(sd(cond_data_finite), 4), "\n")
  } else {
    cat("    ", cond, ": No finite values\n")
  }
}
cat("\n")
      dplyr::filter(!is.na(stat) & !is.na(entrez_final) & entrez_final != "NA" & entrez_final != "") %>%
      dplyr::mutate(entrez_final = as.character(entrez_final)) %>%
      dplyr::arrange(desc(stat)) %>%
      dplyr::select(entrez_final, stat) %>%
      dplyr::distinct(entrez_final, .keep_all = TRUE) %>%
      tibble::deframe()
    
    if (length(ranked_genes) > 0) {
      cat("    Running GSEA (this may take 5-10 minutes)...\n")
      fgsea_res <- fgsea::fgsea(
        pathways = hallmark_list,
        stats = ranked_genes,
        minSize = 15,
        maxSize = 500,
        nperm = 10000,
        nproc = 1  # Use single core to avoid issues
      ) %>%
        dplyr::arrange(padj)
      cat("    ✅ GSEA completed\n")
      
      readr::write_tsv(fgsea_res, 
                      file.path(paths$enrich_dir, "FGSEA_hallmark.tsv"))
      cat("    Saved GSEA Hallmark results:", sum(fgsea_res$padj < 0.05, na.rm = TRUE), 
          "significant pathways\n")
      
      # Plot top pathways
      topPathwaysUp <- fgsea_res %>% 
        dplyr::filter(ES > 0, padj < 0.05) %>% 
        head(10) %>% 
        dplyr::pull(pathway)
      topPathwaysDown <- fgsea_res %>% 
        dplyr::filter(ES < 0, padj < 0.05) %>% 
        head(10) %>% 
        dplyr::pull(pathway)
      
      if (length(topPathwaysUp) > 0) {
        pdf(file.path(paths$enrich_dir, "FGSEA_hallmark_topPathways_up.pdf"), 
            width = 14, height = 10)
        fgsea::plotGseaTable(hallmark_list[topPathwaysUp], ranked_genes, fgsea_res, 
                            gseaParam = 1)
        dev.off()
        cat("    Saved GSEA plot for", length(topPathwaysUp), "up-regulated pathways\n")
      }
      
      if (length(topPathwaysDown) > 0) {
        pdf(file.path(paths$enrich_dir, "FGSEA_hallmark_topPathways_down.pdf"), 
            width = 14, height = 10)
        fgsea::plotGseaTable(hallmark_list[topPathwaysDown], ranked_genes, fgsea_res, 
                            gseaParam = 1)
        dev.off()
        cat("    Saved GSEA plot for", length(topPathwaysDown), "down-regulated pathways\n")
      }
      
      # Enrichment plot for top pathway
      if (nrow(fgsea_res) > 0 && sum(fgsea_res$padj < 0.05, na.rm = TRUE) > 0) {
        top_pathway <- fgsea_res$pathway[which.min(fgsea_res$padj)]
        tryCatch({
          pdf(file.path(paths$enrich_dir, "FGSEA_hallmark_enrichment_plot.pdf"), 
              width = 10, height = 8)
          g_enrich <- enrichplot::gseaplot2(fgsea_res, 
                                            geneSetID = top_pathway,
                                            title = top_pathway)
          print(g_enrich)
          dev.off()
        }, error = function(e) {
          cat("    Warning: Could not create enrichment plot:", e$message, "\n")
        })
      }
    } else {
      cat("    Warning: No ranked genes available for GSEA\n")
    }
  }, error = function(e) {
    cat("    Error in GSEA:", e$message, "\n")
  })
  
} else {
  warning("Could not perform enrichment analysis: insufficient genes with Entrez IDs")
  cat("  sig_entrez:", length(sig_entrez), "\n")
  cat("  gene_universe:", length(gene_universe), "\n")
}

cat("  Enrichment results saved to:", paths$enrich_dir, "\n\n")


# 8.10. Prepare files for DAVID Functional Annotation Tool -------------------------
cat("  8.10. Preparing gene lists for DAVID Functional Annotation Tool...\n")

david_dir <- file.path(paths$results_dir, "david")
fs::dir_create(david_dir)

# Function to prepare DAVID input files
prepare_david_files <- function(deg_data, prefix, description) {
  # Filter genes with valid IDs
  deg_clean <- deg_data %>%
    dplyr::filter(!is.na(entrez_final) & entrez_final != "" & entrez_final != "NA")
  
  if (nrow(deg_clean) == 0) {
    cat("    No genes with Entrez IDs for", description, "\n")
    return(NULL)
  }
  
  # 1. Simple gene list (one ID per line) - Entrez IDs
  entrez_list <- deg_clean %>%
    dplyr::filter(!is.na(entrez_final)) %>%
    dplyr::pull(entrez_final) %>%
    unique()
  
  if (length(entrez_list) > 0) {
    readr::write_lines(entrez_list, 
                       file.path(david_dir, paste0(prefix, "_EntrezID_list.txt")))
    cat("    Saved", length(entrez_list), "Entrez IDs to", 
        paste0(prefix, "_EntrezID_list.txt"), "\n")
  }
  
  # 2. Gene list with Gene Symbols
  symbol_list <- deg_clean %>%
    dplyr::filter(!is.na(symbol) & symbol != "" & symbol != "NA") %>%
    dplyr::pull(symbol) %>%
    unique()
  
  if (length(symbol_list) > 0) {
    readr::write_lines(symbol_list, 
                       file.path(david_dir, paste0(prefix, "_GeneSymbol_list.txt")))
    cat("    Saved", length(symbol_list), "Gene Symbols to", 
        paste0(prefix, "_GeneSymbol_list.txt"), "\n")
  }
  
  # 3. Gene list with Ensembl IDs
  ensembl_list <- deg_clean %>%
    dplyr::filter(!is.na(ensgene) & ensgene != "" & ensgene != "NA") %>%
    dplyr::pull(ensgene) %>%
    unique()
  
  if (length(ensembl_list) > 0) {
    readr::write_lines(ensembl_list, 
                       file.path(david_dir, paste0(prefix, "_EnsemblID_list.txt")))
    cat("    Saved", length(ensembl_list), "Ensembl IDs to", 
        paste0(prefix, "_EnsemblID_list.txt"), "\n")
  }
  
  # 4. DAVID format with fold change (for functional annotation clustering)
  # Format: GeneID\tFoldChange
  david_fc <- deg_clean %>%
    dplyr::select(entrez_final, log2FoldChange) %>%
    dplyr::filter(!is.na(entrez_final) & !is.na(log2FoldChange)) %>%
    dplyr::distinct(entrez_final, .keep_all = TRUE) %>%
    dplyr::arrange(desc(abs(log2FoldChange)))
  
  if (nrow(david_fc) > 0) {
    readr::write_tsv(david_fc, 
                     file.path(david_dir, paste0(prefix, "_DAVID_FC_format.txt")),
                     col_names = FALSE)
    cat("    Saved", nrow(david_fc), "genes with fold changes to", 
        paste0(prefix, "_DAVID_FC_format.txt"), "\n")
  }
  
  # 5. Complete table with all IDs for reference
  david_complete <- deg_clean %>%
    dplyr::select(entrez_final, symbol, ensgene, log2FoldChange, padj, description) %>%
    dplyr::arrange(padj)
  
  readr::write_tsv(david_complete, 
                   file.path(david_dir, paste0(prefix, "_complete_table.tsv")))
  cat("    Saved complete table to", paste0(prefix, "_complete_table.tsv"), "\n")
  
  return(list(
    entrez_count = length(entrez_list),
    symbol_count = length(symbol_list),
    ensembl_count = length(ensembl_list)
  ))
}

# Prepare files for all DEGs
if (nrow(deg_tbl) > 0) {
  cat("    Preparing files for all DEGs...\n")
  prepare_david_files(deg_tbl, "DAVID_all_DEGs", "all DEGs")
}

# Prepare files for up-regulated genes
if (nrow(deg_up) > 0) {
  cat("    Preparing files for up-regulated DEGs...\n")
  prepare_david_files(deg_up, "DAVID_upregulated", "up-regulated DEGs")
}

# Prepare files for down-regulated genes
if (nrow(deg_down) > 0) {
  cat("    Preparing files for down-regulated DEGs...\n")
  prepare_david_files(deg_down, "DAVID_downregulated", "down-regulated DEGs")
}

# Create a README file with instructions for using DAVID
david_readme <- c(
  "DAVID Functional Annotation Tool - Input Files",
  "==============================================",
  "",
  "This directory contains gene lists prepared for DAVID analysis.",
  "",
  "FILE FORMATS:",
  "-------------",
  "",
  "1. *_EntrezID_list.txt:",
  "   - Simple list of Entrez Gene IDs (one per line)",
  "   - Use this in DAVID: Select 'ENTREZ_GENE_ID' as identifier",
  "",
  "2. *_GeneSymbol_list.txt:",
  "   - Simple list of Official Gene Symbols (one per line)",
  "   - Use this in DAVID: Select 'OFFICIAL_GENE_SYMBOL' as identifier",
  "",
  "3. *_EnsemblID_list.txt:",
  "   - Simple list of Ensembl Gene IDs (one per line)",
  "   - Use this in DAVID: Select 'ENSEMBL_GENE_ID' as identifier",
  "",
  "4. *_DAVID_FC_format.txt:",
  "   - Tab-separated file: EntrezID\\tFoldChange",
  "   - Use for Functional Annotation Clustering with fold change weighting",
  "",
  "5. *_complete_table.tsv:",
  "   - Complete table with all annotations for reference",
  "",
  "HOW TO USE DAVID:",
  "-----------------",
  "",
  "1. Go to: https://david.ncifcrf.gov/",
  "2. Click 'Start Analysis'",
  "3. Upload your gene list file (choose appropriate format)",
  "4. Select the correct identifier type:",
  "   - For *_EntrezID_list.txt: Select 'ENTREZ_GENE_ID'",
  "   - For *_GeneSymbol_list.txt: Select 'OFFICIAL_GENE_SYMBOL'",
  "   - For *_EnsemblID_list.txt: Select 'ENSEMBL_GENE_ID'",
  "5. Select species: Homo sapiens",
  "6. Click 'Submit List'",
  "7. Select annotation categories:",
  "   - GOTERM_BP_DIRECT (GO Biological Process)",
  "   - GOTERM_CC_DIRECT (GO Cellular Component)",
  "   - GOTERM_MF_DIRECT (GO Molecular Function)",
  "   - KEGG_PATHWAY",
  "   - BIOCARTA",
  "   - REACTOME_PATHWAY",
  "   - etc.",
  "8. Click 'Functional Annotation Clustering' for grouped results",
  "9. Click 'Functional Annotation Chart' for detailed results",
  "",
  "RECOMMENDED SETTINGS:",
  "--------------------",
  "- EASE Score (p-value threshold): 0.1 (default)",
  "- Count: 2 (minimum genes per term)",
  "- Classification Stringency: Medium",
  "",
  "FILES AVAILABLE:",
  "----------------",
  paste0("- DAVID_all_DEGs_*: All significant DEGs (FDR ≤ 0.05, |log2FC| ≥ 1)"),
  paste0("- DAVID_upregulated_*: Up-regulated DEGs (CBD > Control)"),
  paste0("- DAVID_downregulated_*: Down-regulated DEGs (CBD < Control)"),
  "",
  "NOTE:",
  "-----",
  "DAVID provides complementary analysis to clusterProfiler.",
  "Compare results from both tools for comprehensive interpretation."
)

readr::write_lines(david_readme, file.path(david_dir, "README_DAVID.txt"))
cat("    Created README_DAVID.txt with instructions\n")

cat("  DAVID input files saved to:", david_dir, "\n\n")


# 9. Analysis of key genes and signatures -----------------------------------------------------
cat("Analyzing key genes and pathway signatures...\n")

# 9.1. Key genes from the paper
key_genes <- c("GSDME", "CASP3", "ATF4", "DDIT3", "IGFBP1", "TRPV3")
# DDIT3 is CHOP

# Find these genes in the results
key_genes_data <- res_tbl %>%
  dplyr::filter(symbol %in% key_genes) %>%
  dplyr::select(symbol, log2FoldChange, padj, baseMean)

if (nrow(key_genes_data) > 0) {
  readr::write_tsv(key_genes_data, 
                  file.path(paths$deseq_dir, "key_genes_expression.tsv"))
  
  # Plot key genes expression
  key_genes_expr <- SummarizedExperiment::assay(vsd)[
    key_genes_data$symbol[key_genes_data$symbol %in% rownames(vsd)], 
  ]
  
  if (nrow(key_genes_expr) > 0) {
    pdf(file.path(paths$deseq_dir, "key_genes_heatmap.pdf"), width = 8, height = 6)
    pheatmap::pheatmap(
      key_genes_expr,
      scale = "row",
      annotation_col = col_data[, c("condition", "cell_line"), drop = FALSE],
      main = ""
    )
    dev.off()
  }
  entrez_list <- deg_clean %>%
    dplyr::filter(!is.na(entrez_final)) %>%
    dplyr::pull(entrez_final) %>%
    unique()
  
  if (length(entrez_list) > 0) {
    readr::write_lines(entrez_list, 
                       file.path(david_dir, paste0(prefix, "_EntrezID_list.txt")))
    cat("    Saved", length(entrez_list), "Entrez IDs to", 
        paste0(prefix, "_EntrezID_list.txt"), "\n")
  }
  
  # 2. Gene list with Gene Symbols
  symbol_list <- deg_clean %>%
    dplyr::filter(!is.na(symbol) & symbol != "" & symbol != "NA") %>%
    dplyr::pull(symbol) %>%
    unique()
  
  if (length(symbol_list) > 0) {
    readr::write_lines(symbol_list, 
                       file.path(david_dir, paste0(prefix, "_GeneSymbol_list.txt")))
    cat("    Saved", length(symbol_list), "Gene Symbols to", 
        paste0(prefix, "_GeneSymbol_list.txt"), "\n")
  }
  
  # 3. Gene list with Ensembl IDs
  ensembl_list <- deg_clean %>%
    dplyr::filter(!is.na(ensgene) & ensgene != "" & ensgene != "NA") %>%
    dplyr::pull(ensgene) %>%
    unique()
  
  if (length(ensembl_list) > 0) {
    readr::write_lines(ensembl_list, 
                       file.path(david_dir, paste0(prefix, "_EnsemblID_list.txt")))
    cat("    Saved", length(ensembl_list), "Ensembl IDs to", 
        paste0(prefix, "_EnsemblID_list.txt"), "\n")
  }
  
  # 4. DAVID format with fold change (for functional annotation clustering)
  # Format: GeneID\tFoldChange
  david_fc <- deg_clean %>%
    dplyr::select(entrez_final, log2FoldChange) %>%
    dplyr::filter(!is.na(entrez_final) & !is.na(log2FoldChange)) %>%
    dplyr::distinct(entrez_final, .keep_all = TRUE) %>%
    dplyr::arrange(desc(abs(log2FoldChange)))
  
  if (nrow(david_fc) > 0) {
    readr::write_tsv(david_fc, 
                     file.path(david_dir, paste0(prefix, "_DAVID_FC_format.txt")),
                     col_names = FALSE)
    cat("    Saved", nrow(david_fc), "genes with fold changes to", 
        paste0(prefix, "_DAVID_FC_format.txt"), "\n")
  }
  
  # 5. Complete table with all IDs for reference
  david_complete <- deg_clean %>%
    dplyr::select(entrez_final, symbol, ensgene, log2FoldChange, padj, description) %>%
    dplyr::arrange(padj)
  
  readr::write_tsv(david_complete, 
                   file.path(david_dir, paste0(prefix, "_complete_table.tsv")))
  cat("    Saved complete table to", paste0(prefix, "_complete_table.tsv"), "\n")
  
  return(list(
    entrez_count = length(entrez_list),
    symbol_count = length(symbol_list),
    ensembl_count = length(ensembl_list)
  ))
}

# Preparar archivos para todos los DEGs
if (nrow(deg_tbl) > 0) {
  cat("    Preparando archivos para todos los DEGs...\n")
  prepare_david_files(deg_tbl, "DAVID_all_DEGs", "all DEGs")
}

# Preparar archivos para genes regulados al alza
if (nrow(deg_up) > 0) {
  cat("    Preparando archivos para DEGs regulados al alza...\n")
  prepare_david_files(deg_up, "DAVID_upregulated", "up-regulated DEGs")
}

# Preparar archivos para genes regulados a la baja
if (nrow(deg_down) > 0) {
  cat("    Preparando archivos para DEGs regulados a la baja...\n")
  prepare_david_files(deg_down, "DAVID_downregulated", "down-regulated DEGs")
}

# Crear un archivo README con instrucciones para usar DAVID
david_readme <- c(
  "Herramienta de Anotación Funcional DAVID - Archivos de Entrada",
  "==============================================================",
  "",
  "Este directorio contiene listas de genes preparadas para el análisis DAVID.",
  "",
  "FORMATOS DE ARCHIVO:",
  "--------------------",
  "",
  "1. *_EntrezID_list.txt:",
  "   - Lista simple de IDs de genes Entrez (uno por línea)",
  "   - Usar en DAVID: Seleccionar 'ENTREZ_GENE_ID' como identificador",
  "",
  "2. *_GeneSymbol_list.txt:",
  "   - Lista simple de Símbolos de Genes Oficiales (uno por línea)",
  "   - Usar en DAVID: Seleccionar 'OFFICIAL_GENE_SYMBOL' como identificador",
  "",
  "3. *_EnsemblID_list.txt:",
  "   - Lista simple de IDs de genes Ensembl (uno por línea)",
  "   - Usar en DAVID: Seleccionar 'ENSEMBL_GENE_ID' como identificador",
  "",
  "4. *_DAVID_FC_format.txt:",
  "   - Archivo separado por tabulaciones: EntrezID\\tFoldChange",
  "   - Usar para Agrupación de Anotación Funcional con ponderación por cambio de expresión",
  "",
  "5. *_complete_table.tsv:",
  "   - Tabla completa con todas las anotaciones para referencia",
  "",
  "CÓMO USAR DAVID:",
  "-----------------",
  "",
  "1. Ir a: https://david.ncifcrf.gov/",
  "2. Hacer clic en 'Start Analysis'",
  "3. Subir el archivo de lista de genes (elegir el formato apropiado)",
  "4. Seleccionar el tipo de identificador correcto:",
  "   - Para *_EntrezID_list.txt: Seleccionar 'ENTREZ_GENE_ID'",
  "   - Para *_GeneSymbol_list.txt: Seleccionar 'OFFICIAL_GENE_SYMBOL'",
  "   - Para *_EnsemblID_list.txt: Seleccionar 'ENSEMBL_GENE_ID'",
  "5. Seleccionar especie: Homo sapiens",
  "6. Hacer clic en 'Submit List'",
  "7. Seleccionar categorías de anotación:",
  "   - GOTERM_BP_DIRECT (Proceso Biológico GO)",
  "   - GOTERM_CC_DIRECT (Componente Celular GO)",
  "   - GOTERM_MF_DIRECT (Función Molecular GO)",
  "   - KEGG_PATHWAY",
  "   - BIOCARTA",
  "   - REACTOME_PATHWAY",
  "   - etc.",
  "8. Hacer clic en 'Functional Annotation Clustering' para resultados agrupados",
  "9. Hacer clic en 'Functional Annotation Chart' para resultados detallados",
  "",
  "CONFIGURACIÓN RECOMENDADA:",
  "--------------------------",
  "- Puntuación EASE (umbral de p-valor): 0.1 (por defecto)",
  "- Conteo: 2 (mínimo de genes por término)",
  "- Rigurosidad de Clasificación: Media",
  "",
  "ARCHIVOS DISPONIBLES:",
  "--------------------",
  paste0("- DAVID_all_DEGs_*: Todos los DEGs significativos (FDR ≤ 0.05, |log2FC| ≥ 1)"),
  paste0("- DAVID_upregulated_*: DEGs regulados al alza (CBD > Control)"),
  paste0("- DAVID_downregulated_*: DEGs regulados a la baja (CBD < Control)"),
  "",
  "NOTA:",
  "-----",
  "DAVID proporciona un análisis complementario a clusterProfiler.",
  "Comparar los resultados de ambas herramientas para una interpretación completa."
)

readr::write_lines(david_readme, file.path(david_dir, "README_DAVID.txt"))
cat("    Creado README_DAVID.txt con instrucciones\n")

cat("  Archivos de entrada DAVID guardados en:", david_dir, "\n\n")


# 9. Análisis de genes clave y firmas -----------------------------------------------------
cat("Analizando genes clave y firmas de vías...\n")

# 9.1. Genes clave del estudio
key_genes <- c("GSDME", "CASP3", "ATF4", "DDIT3", "IGFBP1", "TRPV3")
# DDIT3 es CHOP

# Encontrar estos genes en los resultados
key_genes_data <- res_tbl %>%
  dplyr::filter(symbol %in% key_genes) %>%
  dplyr::select(symbol, log2FoldChange, padj, baseMean)

if (nrow(key_genes_data) > 0) {
  readr::write_tsv(key_genes_data, 
                  file.path(paths$deseq_dir, "key_genes_expression.tsv"))
  
  # Graficar la expresión de genes clave
  key_genes_expr <- SummarizedExperiment::assay(vsd)[
    key_genes_data$symbol[key_genes_data$symbol %in% rownames(vsd)], 
  ]
  
  if (nrow(key_genes_expr) > 0) {
    pdf(file.path(paths$deseq_dir, "key_genes_heatmap.pdf"), width = 8, height = 6)
    pheatmap::pheatmap(
      key_genes_expr,
      scale = "row",
      annotation_col = col_data[, c("condition", "cell_line"), drop = FALSE],
      main = ""
    )
    dev.off()
  }
}

# 9.2. Firmas de vías (Piroptosis, ISR, Glicólisis)
# Inicializar una lista para almacenar los conjuntos de genes clave
key_genes_list <- list()

# 9.2.1. Piroptosis
pyroptosis_genes <- c("GSDME", "GSDMD", "CASP1", "CASP3", "CASP4", "CASP5", 
                      "NLRP3", "PYCARD", "IL1B", "IL18")
key_genes_list[["Pyroptosis"]] <- pyroptosis_genes

# 9.2.2. Respuesta Integrada al Estrés (ISR)
isr_genes <- c("EIF2AK1", "EIF2AK2", "EIF2AK3", "EIF2AK4", "ATF4", "DDIT3", "PPP1R15A", "ASNS")
key_genes_list[["ISR"]] <- isr_genes

# 9.2.3. Glicólisis (Efecto Warburg)
glycolysis_genes <- c("HK2", "PFKP", "ALDOA", "GAPDH", "PGK1", "ENO1", "PKM", "LDHA", "SLC2A1")
key_genes_list[["Glycolysis"]] <- glycolysis_genes

# Analizar cada conjunto de genes clave
for (pathway_name in names(key_genes_list)) {
  genes <- key_genes_list[[pathway_name]]
  
  # Filtrar genes que están en nuestros datos
  valid_genes <- intersect(genes, rownames(vsd))
  
  if (length(valid_genes) >= 3) {
    cat("  Analizando firma de", pathway_name, "(", length(valid_genes), "genes encontrados)...\n")
    
    # 1. Heatmap de estos genes
    pathway_expr <- SummarizedExperiment::assay(vsd)[valid_genes, , drop = FALSE]
    
    # Escalar por fila para visualización
    pdf(file.path(paths$enrich_dir, paste0("heatmap_", pathway_name, ".pdf")), 
        width = 8, height = max(6, length(valid_genes) * 0.2 + 2))
    
    pheatmap::pheatmap(
      pathway_expr,
      scale = "row",
      clustering_distance_rows = "correlation",
      clustering_distance_cols = "correlation",
      annotation_col = col_data[, c("condition", "cell_line"), drop = FALSE],
      main = paste("Firma:", pathway_name),
      fontsize_row = 8
    )
    dev.off()
    
    # 2. Boxplot de puntuación de firma (promedio de expresión escalada Z)
    # Calcular Z-scores para cada gen
    z_scores <- t(scale(t(pathway_expr)))
    
    # Calcular puntuación promedio por muestra
    signature_scores <- colMeans(z_scores, na.rm = TRUE)
    
    # Crear data frame para graficar
    plot_data <- data.frame(
      Sample = names(signature_scores),
      Score = signature_scores,
      Condition = col_data$condition,
      CellLine = col_data$cell_line
    )
    
    # Boxplot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Condition, y = Score, fill = Condition)) +
      ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      ggplot2::geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
      ggplot2::facet_wrap(~CellLine) +
      ggplot2::scale_fill_manual(values = c("Control" = "gray70", "CBD" = "forestgreen")) +
      ggplot2::labs(
        title = paste("Puntuación de Firma:", pathway_name),
        subtitle = paste(length(valid_genes), "genes"),
        y = "Puntuación Z-score Promedio",
        x = "Condición"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none")
    
    ggplot2::ggsave(file.path(paths$enrich_dir, paste0("boxplot_score_", pathway_name, ".pdf")), 
           p, width = 6, height = 5)
    
    # Prueba estadística para diferencia de puntuación
    # ANOVA de dos vías para probar efectos de Condición y Línea Celular
    tryCatch({
      model <- lm(Score ~ Condition * CellLine, data = plot_data)
      anova_res <- anova(model)
      
      cat("    Resultados ANOVA para puntuación de", pathway_name, ":\n")
      print(anova_res)
      
      # Guardar resultados estadísticos
      capture.output(print(anova_res), 
                     file = file.path(paths$enrich_dir, paste0("stats_", pathway_name, ".txt")))
    }, error = function(e) {
      cat("    No se pudo realizar prueba estadística para", pathway_name, ":", e$message, "\n")
    })
    
  } else {
    cat("  Saltando", pathway_name, "- insuficientes genes encontrados (", length(valid_genes), ")\n")
  }
}

cat("  Análisis de genes clave y firmas completado.\n\n")


# 10. Guardar información de sesión -----------------------------------------------------------
cat("Guardando información de la sesión...\n")
writeLines(capture.output(sessionInfo()), 
          file.path(paths$results_dir, "session_info.txt"))

cat("\n================================================================\n")
cat("  ANÁLISIS COMPLETADO EXITOSAMENTE\n")
cat("================================================================\n")
cat("Resumen de resultados:\n")
cat("  - Directorio de salida: ", paths$results_dir, "\n")
cat("  - Resultados DESeq2: ", paths$deseq_dir, "\n")
cat("  - Gráficos QC: ", paths$qc_dir, "\n")
cat("  - Resultados de enriquecimiento: ", paths$enrich_dir, "\n")

# Fin del pipeline --------------------------------------------------------------------------------
