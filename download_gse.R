# download_gse.R
library(GEOquery)

# Función para descargar GSEs y guardar matrices y metadata
download_gse_list <- function(gse_list, base_dir = "differential-expression-uterine-cancer/DATA") {
  
  for (gse_id in gse_list) {
    cat("\n⬇️ Descargando", gse_id, "...\n")
    
    # Carpeta destino
    dest_dir <- file.path(base_dir, gse_id)
    if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
    
    # Descargar GSE (Series Matrix + metadata)
    gse <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE, destdir = dest_dir)
    if (length(gse) == 0) {
      warning(paste("No se pudo descargar", gse_id))
      next
    }
    
    gse_data <- gse[[1]]
    
    # Matriz de expresión
    exprs_mat <- exprs(gse_data)
    exprs_file <- file.path(dest_dir, paste0(gse_id, "_exprs_matrix.csv"))
    write.csv(exprs_mat, file = exprs_file, row.names = TRUE)
    cat("💾 Matriz de expresión guardada en", exprs_file, "\n")
    
    # Metadata de muestras
    pheno <- pData(gse_data)
    pheno_file <- file.path(dest_dir, paste0(gse_id, "_pheno_data.csv"))
    write.csv(pheno, file = pheno_file, row.names = TRUE)
    cat("💾 Metadata guardada en", pheno_file, "\n")
  }
  
  cat("\n✅ Descarga completada.\n")
}
