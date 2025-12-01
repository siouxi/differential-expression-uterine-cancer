# 🧬 Análisis de Expresión Diferencial en Cáncer Uterino

Este proyecto contiene análisis de expresión diferencial de genes (DEGs) a partir de múltiples datasets de GEO (Gene Expression Omnibus) relacionados con cáncer uterino y otros tipos de cáncer.

## 📋 Tabla de Contenidos

- [Estructura del Proyecto](#estructura-del-proyecto)
- [Datasets Analizados](#datasets-analizados)
- [Archivos y Directorios](#archivos-y-directorios)
- [Scripts de Python](#scripts-de-python)
- [Análisis Realizados](#análisis-realizados)
- [Diagrama de Venn Interactivo](#diagrama-de-venn-interactivo)
- [Requisitos](#requisitos)
- [Uso](#uso)

## 📁 Estructura del Proyecto

```
differential-expression-uterine-cancer/
├── analysis/           # Scripts R para análisis por dataset
├── DATA/              # Datos crudos descargados
├── ESTUDIOS/          # Documentación de estudios
├── results/           # Resultados de análisis DESeq2
├── standard/          # Archivos estandarizados con IDs normalizados
├── venn/              # Datos para diagramas de Venn
├── *.py               # Scripts Python de procesamiento
└── venn_diagram_interactive.html  # Visualización interactiva
```

## 🔬 Datasets Analizados

El proyecto analiza los siguientes datasets de GEO:

| Dataset | Descripción | Líneas Celulares/Condiciones |
|---------|-------------|------------------------------|
| **GSE21656** | Expresión diferencial en cáncer | Múltiples condiciones |
| **GSE102787** | Análisis de expresión génica | Múltiples condiciones |
| **GSE131565** | Comparación CBD vs Control | Tratamiento con CBD |
| **GSE173201** | Expresión diferencial | Múltiples condiciones |
| **GSE179661** | Análisis en líneas celulares | HepG2, MHCC97H |
| **GSE197561** | Expresión diferencial | Múltiples condiciones |
| **GSE223827** | Análisis de expresión | Múltiples condiciones |
| **GSE285498** | Expresión diferencial | Múltiples condiciones |

## 📂 Archivos y Directorios

### Directorios Principales

#### `analysis/`
Contiene scripts R para el análisis de expresión diferencial usando DESeq2. Cada subdirectorio corresponde a un dataset específico (GSE) y contiene:
- Scripts R de pipeline de análisis
- Gráficos de resultados (PCA, heatmaps, volcano plots)
- Análisis de enriquecimiento funcional

#### `results/`
Almacena los resultados crudos de los análisis DESeq2:
- Archivos TSV/CSV con genes diferencialmente expresados
- Estadísticas de expresión (log2FoldChange, p-values, FDR)
- Anotaciones de genes
- Archivos de normalización

#### `standard/`
Contiene archivos **estandarizados** con identificadores de genes normalizados. Estos archivos tienen el sufijo `_standard.csv` y contienen:
- **Ensembl_ID**: Identificador Ensembl del gen
- **Entrez_ID**: Identificador Entrez del gen
- **Gene_Symbol**: Símbolo oficial del gen
- **log2FoldChange**: Cambio de expresión (log2)

**Propósito**: Estos archivos permiten comparar genes entre diferentes estudios usando identificadores consistentes.

#### `venn/`
Directorio para almacenar datos procesados para diagramas de Venn (actualmente vacío, los datos se cargan directamente desde `standard/`).

### Archivos HTML

#### `venn_diagram_interactive.html`
**Visualización interactiva de diagramas de Venn** para comparar genes diferencialmente expresados entre múltiples datasets.

**Características:**
- Carga hasta 6 datasets simultáneamente
- Selección de tipo de ID (Ensembl, Entrez, Gene Symbol)
- Filtrado por regulación (todos, up-regulados, down-regulados)
- Visualización de intersecciones entre datasets
- Exportación de resultados a CSV
- Leyenda de colores para cada dataset

**Uso:**
1. Abrir el archivo en un navegador web
2. Cargar archivos `_standard.csv` desde el directorio `standard/`
3. Asignar nombres personalizados a cada dataset
4. Seleccionar tipo de ID y filtros
5. Generar el diagrama
6. Hacer clic en intersecciones para ver genes específicos

## 🐍 Scripts de Python

### Scripts de Estandarización (MyGene)

Estos scripts convierten identificadores de genes a un formato estándar usando la API de MyGene.info:

#### `standardize_gse21656_mygene.py`
- **Entrada**: `standard/GSE21656/GSE21656_DEGS_significativos.csv`
- **Proceso**: Mapea Gene Symbols → Ensembl ID, Entrez ID, Gene Symbol
- **Salida**: `GSE21656_DEGS_significativos_standard.csv`

#### `standardize_gse102787_mygene.py`
- **Entrada**: `standard/GSE102787/GSE102787_DEGS_Significativos.csv`
- **Proceso**: Mapea Gene Symbols → Ensembl ID, Entrez ID, Gene Symbol
- **Salida**: `GSE102787_DEGS_Significativos_standard.csv`

#### `standardize_gse131565_mygene.py`
- **Entrada**: `standard/GSE131565/DESeq2_CBD_vs_Control.csv`
- **Proceso**: Mapea Entrez IDs → Ensembl ID, Entrez ID, Gene Symbol
- **Salida**: `DESeq2_CBD_vs_Control_standard.csv`

#### `standardize_gse173201_mygene.py`
- **Entrada**: Archivos CSV en `standard/GSE173201/`
- **Proceso**: Mapea IDs numéricos (Entrez) → Ensembl ID, Entrez ID, Gene Symbol
- **Salida**: Archivos `*_standard.csv`

#### `standardize_gse179661_mygene.py`
- **Entrada**: 
  - `standard/GSE179661/DESeq2_DEG_FDR0.05_log2FC2_HepG2.csv`
  - `standard/GSE179661/DESeq2_DEG_FDR0.05_log2FC2_MHCC97H.csv`
- **Proceso**: Mapea Entrez IDs → Ensembl ID, Entrez ID, Gene Symbol
- **Características especiales**: 
  - Usa MyGene.info como fuente principal
  - Usa archivo local de anotación como respaldo
  - Reporta estadísticas detalladas de mapeo
- **Salida**: Archivos `*_standard.csv`

#### `standardize_gse197561_mygene.py`
- **Entrada**: Archivos en `standard/GSE197561/`
- **Proceso**: Estandarización de IDs
- **Salida**: Archivos `*_standard.csv`

#### `standardize_gse223827_mygene.py`
- **Entrada**: Múltiples archivos CSV en `standard/GSE223827/`
- **Proceso**: Mapea IDs → Ensembl ID, Entrez ID, Gene Symbol
- **Salida**: Archivos `*_standard.csv`

#### `standardize_gse285498_mygene.py`
- **Entrada**: Archivos en `standard/GSE285498/`
- **Proceso**: Estandarización de IDs
- **Salida**: Archivos `*_standard.csv`

### Scripts de Procesamiento Inicial

#### `process_gse131565.py`
Procesa archivos originales del dataset GSE131565 y los prepara para estandarización.

#### `process_gse173201.py`
Procesa archivos originales del dataset GSE173201.

#### `process_gse179661.py`
Procesa archivos originales del dataset GSE179661 (HepG2 y MHCC97H).

#### `process_gse223827.py`
Procesa archivos originales del dataset GSE223827.

#### `process_gse285498.py`
Procesa archivos originales del dataset GSE285498.

#### `process_excel_files.py`
Procesa archivos Excel y los convierte a formato CSV.

### Scripts de Utilidad

#### `fix_entrez_decimals.py`
**Propósito**: Corrige IDs de Entrez que tienen formato decimal (ej: `3.0` → `3`)

**Uso:**
```python
# Convierte Entrez_ID de formato float a integer
# Ejemplo: 12345.0 → 12345
```

#### `check_ids_and_symbols.py`
**Propósito**: Verifica la validez de IDs de genes y símbolos

**Funciones:**
- Verifica IDs en archivo local de anotación
- Consulta MyGene.info para validar símbolos
- Genera reporte en `check_results.txt`

#### `inspect_ids.py`
**Propósito**: Inspecciona IDs específicos de genes usando MyGene.info

**Uso:**
```python
# Verifica una lista de IDs de Entrez
# Genera inspection_results.txt con detalles
```

#### `inspect_gse173201.py`
Inspecciona estructura de archivos del dataset GSE173201.

#### `inspect_gse179661.py`
Inspecciona estructura de archivos del dataset GSE179661.

#### `inspect_gse223827.py`
Inspecciona estructura de archivos del dataset GSE223827.

#### `inspect_gse285498.py`
Inspecciona estructura de archivos del dataset GSE285498.

#### `inspect_excel.py`
Inspecciona contenido de archivos Excel.

### Scripts de Descarga

#### `gse_downloader.py`
Descarga datasets de GEO usando GEOparse.

#### `download_gse.R`
Script R para descargar datasets de GEO.

#### `extract_pdf_text.py`
Extrae texto de archivos PDF (para documentación de estudios).

## 📊 Análisis Realizados

### 1. Análisis de Expresión Diferencial (DESeq2)

Cada dataset se analiza usando DESeq2 con los siguientes criterios:

- **FDR (False Discovery Rate)**: < 0.05
- **log2FoldChange**: ≥ 2 o ≤ -2 (en la mayoría de los casos)
- **Normalización**: DESeq2 normalization

**Resultados incluyen:**
- Lista de genes diferencialmente expresados (DEGs)
- Estadísticas (p-value, adjusted p-value, log2FC)
- Gráficos de volcano plot
- Heatmaps de expresión
- Análisis PCA

### 2. Diagramas de Venn

**Propósito**: Identificar genes comunes y únicos entre diferentes estudios.

**Tipos de análisis:**
- **Intersección total**: Genes presentes en todos los datasets
- **Intersecciones parciales**: Genes compartidos entre subconjuntos de datasets
- **Genes únicos**: Genes presentes solo en un dataset

**Filtros disponibles:**
- Por tipo de ID (Ensembl, Entrez, Gene Symbol)
- Por regulación (up, down, todos)

### 3. DEGs (Genes Diferencialmente Expresados)

**Definición**: Genes cuya expresión cambia significativamente entre condiciones.

**Clasificación:**
- **Up-regulados**: log2FoldChange > 0 (mayor expresión en condición tratada)
- **Down-regulados**: log2FoldChange < 0 (menor expresión en condición tratada)

**Archivos de DEGs:**
- Formato: CSV con columnas estandarizadas
- Ubicación: `standard/GSE_ID/*_standard.csv`
- Contenido: Ensembl_ID, Entrez_ID, Gene_Symbol, log2FoldChange

## 🎨 Diagrama de Venn Interactivo

El archivo `venn_diagram_interactive.html` proporciona una herramienta visual para:

### Funcionalidades

1. **Carga de Datasets**
   - Hasta 6 datasets simultáneamente
   - Soporte para archivos CSV estandarizados
   - Nombres personalizables para cada dataset

2. **Configuración de Análisis**
   - **Tipo de ID**: Ensembl_ID, Entrez_ID, o Gene_Symbol
   - **Filtro de regulación**:
     - Todos los genes
     - Solo up-regulados (log2FC > 0)
     - Solo down-regulados (log2FC < 0)

3. **Visualización**
   - Diagramas de Venn para 2, 3, o 4 datasets
   - Vista de lista para 5-6 datasets
   - Etiquetas con número de genes en cada intersección
   - Leyenda de colores por dataset

4. **Interactividad**
   - Click en intersecciones para ver lista de genes
   - Tabla con detalles de cada gen
   - Indicadores de regulación (up/down)
   - Exportación a CSV de genes seleccionados

### Ejemplo de Uso

```
1. Abrir venn_diagram_interactive.html en navegador
2. Cargar archivos:
   - Dataset 1: standard/GSE179661/DESeq2_DEG_FDR0.05_log2FC2_HepG2_standard.csv
   - Dataset 2: standard/GSE179661/DESeq2_DEG_FDR0.05_log2FC2_MHCC97H_standard.csv
3. Nombrar datasets: "HepG2", "MHCC97H"
4. Seleccionar: Tipo ID = "Ensembl_ID", Regulación = "Todos"
5. Click en "Generar Diagrama"
6. Explorar intersecciones y exportar resultados
```

## 📦 Requisitos

### Python
```bash
pip install pandas
pip install mygene
pip install openpyxl  # Para archivos Excel
pip install PyPDF2    # Para extract_pdf_text.py
```

### R
```r
install.packages("DESeq2")
install.packages("GEOquery")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("EnhancedVolcano")
```

## 🚀 Uso

### 1. Estandarizar Nuevos Datasets

```bash
# Ejemplo para GSE179661
python standardize_gse179661_mygene.py
```

### 2. Generar Diagramas de Venn

```bash
# Abrir en navegador
start venn_diagram_interactive.html
```

### 3. Ejecutar Análisis DESeq2

```r
# En R
source("analysis/GSE179661/GSE179661_HepG2_pipeline.R")
```

## 📝 Notas Importantes

### Sobre los IDs de Genes

- **Ensembl_ID**: Formato `ENSG00000000003` - Más estable, recomendado para análisis
- **Entrez_ID**: Formato numérico `3` - Ampliamente usado en bases de datos
- **Gene_Symbol**: Formato `A2M` - Más legible pero puede cambiar

### Sobre los Archivos Estandarizados

Los archivos `*_standard.csv` son el resultado de:
1. Extracción de DEGs de análisis DESeq2
2. Mapeo de IDs usando MyGene.info
3. Validación con archivos de anotación locales
4. Normalización de formato

**Ventajas:**
- Comparabilidad entre estudios
- Múltiples tipos de ID disponibles
- Validación de IDs
- Formato consistente

### Sobre el Log2FoldChange

- **Positivo**: Gen up-regulado (más expresión en tratamiento)
- **Negativo**: Gen down-regulado (menos expresión en tratamiento)
- **Magnitud**: Indica cuánto cambió la expresión
  - log2FC = 1 → 2x más expresión
  - log2FC = 2 → 4x más expresión
  - log2FC = -1 → 2x menos expresión

## 📧 Contacto

Para preguntas o sugerencias sobre este proyecto, por favor contactar al investigador principal.

## 📄 Licencia

Ver archivo `LICENSE` para detalles.

---

**Última actualización**: Diciembre 2025
