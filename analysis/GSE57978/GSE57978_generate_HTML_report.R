# ==============================================================================
# Script para generar reporte HTML de resumen de QC y DEGs
# Dataset: GSE57978
# ==============================================================================

cat("\n", rep("=", 70), "\n")
cat("GENERANDO REPORTE HTML DE RESUMEN\n")
cat(rep("=", 70), "\n\n")

# Obtener directorio del script
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

if (is.null(script_dir) || length(script_dir) == 0 || script_dir == "") {
    current_dir <- getwd()
    if (basename(current_dir) == "GSE57978") {
        script_dir <- current_dir
    } else {
        script_dir <- file.path(current_dir, "ANALISIS", "GSE57978")
    }
}

output_dir <- script_dir
qc_dir <- file.path(output_dir, "QC")
degs_dir <- file.path(output_dir, "DEGs")

# Cargar datos necesarios
cat("Cargando datos...\n")

# Cargar QC metrics
qc_metrics_file <- file.path(qc_dir, "QC_metrics_summary.csv")
if (file.exists(qc_metrics_file)) {
    qc_metrics <- read.csv(qc_metrics_file, stringsAsFactors = FALSE)
    cat("  ✓ QC metrics cargados\n")
} else {
    cat("  ⚠ No se encontró QC_metrics_summary.csv\n")
    qc_metrics <- NULL
}

# Cargar DEGs summary
degs_summary_file <- file.path(degs_dir, "DEGs_summary.csv")
if (file.exists(degs_summary_file)) {
    degs_summary <- read.csv(degs_summary_file, stringsAsFactors = FALSE)
    cat("  ✓ DEGs summary cargado\n")
} else {
    cat("  ⚠ No se encontró DEGs_summary.csv\n")
    degs_summary <- NULL
}

# Cargar sample info
sample_info_file <- file.path(degs_dir, "sample_info.csv")
if (file.exists(sample_info_file)) {
    sample_info <- read.csv(sample_info_file, stringsAsFactors = FALSE)
    cat("  ✓ Sample info cargado\n")
} else {
    cat("  ⚠ No se encontró sample_info.csv\n")
    sample_info <- NULL
}

# Cargar DEGs
degs_file <- file.path(degs_dir, "all_DEGs.tsv")
if (file.exists(degs_file)) {
    degs <- read.table(degs_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    cat("  ✓ DEGs cargados\n")
} else {
    cat("  ⚠ No se encontró all_DEGs.tsv\n")
    degs <- NULL
}

# Cargar up y down genes
up_genes_file <- file.path(degs_dir, "up_genes_complete.tsv")
down_genes_file <- file.path(degs_dir, "down_genes_complete.tsv")

if (file.exists(up_genes_file)) {
    up_genes <- read.table(up_genes_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    cat("  ✓ Up-regulated genes cargados\n")
} else {
    up_genes <- NULL
}

if (file.exists(down_genes_file)) {
    down_genes <- read.table(down_genes_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    cat("  ✓ Down-regulated genes cargados\n")
} else {
    down_genes <- NULL
}

# Cargar filtering summary
filtering_summary_file <- file.path(output_dir, "GSE57978_filtering_summary.csv")
if (file.exists(filtering_summary_file)) {
    filtering_summary <- read.csv(filtering_summary_file, stringsAsFactors = FALSE)
    cat("  ✓ Filtering summary cargado\n")
} else {
    filtering_summary <- NULL
}

cat("\nGenerando contenido HTML...\n")

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
        
        .image-container {
            text-align: center;
            margin: 20px 0;
        }
        
        .image-container img {
            max-width: 100%;
            height: auto;
            border-radius: 8px;
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
if (!is.null(sample_info)) {
    html_content <- paste0(html_content, '
        <h2>🧬 Información de Muestras</h2>
        <div class="section">
            <p><strong>Total de muestras:</strong> ', nrow(sample_info), '</p>
            <table>
                <tr>
                    <th>Muestra</th>
                    <th>GSM ID</th>
                    <th>Línea Celular</th>
                    <th>Tratamiento</th>
                </tr>')
    
    for (i in 1:nrow(sample_info)) {
        html_content <- paste0(html_content, '<tr>
            <td>', sample_info$Sample[i], '</td>
            <td>', sample_info$GSM_ID[i], '</td>
            <td>', sample_info$Cell_Line[i], '</td>
            <td>', sample_info$Treatment[i], '</td>
        </tr>')
    }
    
    html_content <- paste0(html_content, '</table>
        </div>')
}

# Sección de QC
html_content <- paste0(html_content, '
        <h2>🔍 Análisis de Calidad (QC)</h2>')

if (!is.null(qc_metrics)) {
    html_content <- paste0(html_content, '
        <div class="summary-grid">
            <div class="info-box">
                <div class="metric-label">Mediana de Intensidad (Log2)</div>
                <div class="metric-value">', round(mean(qc_metrics$Median_Intensity, na.rm = TRUE), 2), '</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Media de Intensidad (Log2)</div>
                <div class="metric-value">', round(mean(qc_metrics$Mean_Intensity, na.rm = TRUE), 2), '</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Desviación Estándar Promedio</div>
                <div class="metric-value">', round(mean(qc_metrics$SD_Intensity, na.rm = TRUE), 2), '</div>
            </div>
            <div class="info-box">
                <div class="metric-label">IQR Promedio</div>
                <div class="metric-value">', round(mean(qc_metrics$IQR_Intensity, na.rm = TRUE), 2), '</div>
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
    
    for (i in 1:nrow(qc_metrics)) {
        html_content <- paste0(html_content, '<tr>
            <td>', qc_metrics$Sample[i], '</td>
            <td>', qc_metrics$Treatment[i], '</td>
            <td>', round(qc_metrics$Median_Intensity[i], 3), '</td>
            <td>', round(qc_metrics$Mean_Intensity[i], 3), '</td>
            <td>', round(qc_metrics$SD_Intensity[i], 3), '</td>
            <td>', round(qc_metrics$IQR_Intensity[i], 3), '</td>
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
if (!is.null(filtering_summary)) {
    html_content <- paste0(html_content, '
        <h2>🔬 Filtrado de Genes</h2>
        <div class="section">
            <table>
                <tr>
                    <th>Criterio</th>
                    <th>Número</th>
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
}

# Sección de DEGs
html_content <- paste0(html_content, '
        <h2>📈 Análisis de Expresión Diferencial (DEGs)</h2>')

if (!is.null(degs_summary)) {
    html_content <- paste0(html_content, '
        <div class="summary-grid">
            <div class="info-box">
                <div class="metric-label">Total Genes Analizados</div>
                <div class="metric-value">', degs_summary$Numero[degs_summary$Categoria == "Total genes analizados"], '</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Total DEGs</div>
                <div class="metric-value">', degs_summary$Numero[degs_summary$Categoria == "DEGs totales"], '</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Up-regulated</div>
                <div class="metric-value">', degs_summary$Numero[degs_summary$Categoria == "Up-regulated (CBD > Vehicle)"], '</div>
            </div>
            <div class="info-box">
                <div class="metric-label">Down-regulated</div>
                <div class="metric-value">', degs_summary$Numero[degs_summary$Categoria == "Down-regulated (CBD < Vehicle)"], '</div>
            </div>
        </div>
        
        <div class="section">
            <table>
                <tr>
                    <th>Categoría</th>
                    <th>Número</th>
                    <th>Porcentaje</th>
                </tr>')
    
    for (i in 1:nrow(degs_summary)) {
        html_content <- paste0(html_content, '<tr>
            <td>', degs_summary$Categoria[i], '</td>
            <td>', degs_summary$Numero[i], '</td>
            <td>', degs_summary$Porcentaje[i], '%</td>
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

if (!is.null(degs) && nrow(degs) > 0) {
    html_content <- paste0(html_content, '
            <a href="DEGs/heatmap_top_genes.png" class="link-button" target="_blank">🔥 Heatmap Top Genes</a>')
}

html_content <- paste0(html_content, '
        </p>')

# Mostrar top DEGs
if (!is.null(degs) && nrow(degs) > 0) {
    # Ordenar por adj.P.Val
    degs_sorted <- degs[order(degs$adj.P.Val), ]
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
        
        html_content <- paste0(html_content, '<tr>
            <td><strong>', ifelse("SYMBOL" %in% colnames(top_degs), top_degs$SYMBOL[i], "N/A"), '</strong></td>
            <td>', ifelse("ENTREZID" %in% colnames(top_degs), top_degs$ENTREZID[i], "N/A"), '</td>
            <td>', round(top_degs$logFC[i], 3), '</td>
            <td>', format(top_degs$P.Value[i], scientific = TRUE, digits = 3), '</td>
            <td>', format(top_degs$adj.P.Val[i], scientific = TRUE, digits = 3), '</td>
            <td style="color: ', regulation_color, '; font-weight: bold;">', regulation, '</td>
        </tr>')
    }
    
    html_content <- paste0(html_content, '</table>
        </div>')
    
    # Separar up y down
    if (!is.null(up_genes) && nrow(up_genes) > 0) {
        html_content <- paste0(html_content, '
        <h3>Genes Up-regulated (CBD > Vehicle)</h3>
        <div class="success">
            <p><strong>Total:</strong> ', nrow(up_genes), ' genes</p>
        </div>')
    }
    
    if (!is.null(down_genes) && nrow(down_genes) > 0) {
        html_content <- paste0(html_content, '
        <h3>Genes Down-regulated (CBD < Vehicle)</h3>
        <div class="success">
            <p><strong>Total:</strong> ', nrow(down_genes), ' genes</p>
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

if (!is.null(up_genes) && nrow(up_genes) > 0) {
    html_content <- paste0(html_content, '
            <a href="DEGs/up_genes.txt" class="link-button">📄 Genes Up-regulated</a>
            <a href="DEGs/up_genes_complete.tsv" class="link-button">📄 Up-regulated Completo (TSV)</a>')
}

if (!is.null(down_genes) && nrow(down_genes) > 0) {
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

