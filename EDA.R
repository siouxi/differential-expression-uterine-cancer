if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")

library(GEOquery)

gse_id <- "GSE285498"

dest_dir <- file.path("DATA", gse_id)

gse <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE, destdir = dest_dir)

gse_data <- gse[[1]]


# Guardar metadata de muestras
pheno <- pData(gse_data)
pheno_file <- file.path(dest_dir, paste0(gse_id, "_pheno_data.csv"))
write.csv(pheno, file = "xd.csv", row.names = FALSE)
