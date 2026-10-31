# Checkpoints de la lección "De las cuentas de RNA-seq a la exploración de datos"
# (vignettes/intro_RNAseq_es.Rmd). Reproduce los mismos pasos de la lección para
# que las lecciones siguientes puedan empezar con `dds` y `zscore_mat` ya listos.
# Ejecutar desde la raíz del paquete: source("data-raw/lesson2_checkpoints.R")
# Si cambian los pasos de la lección, actualizar este script y volver a ejecutarlo.

library(DESeq2)
filtered_txi <- GeneExprLupusMx::download_data("filtered_txi") # ver ?download_data
load("data/metadata_cleaned.rda")

# Metadatos: voluntaria y niveles de referencia
metadata_cleaned$Volunteer_ID <- factor(sub("_.*", "", metadata_cleaned$sample_ID))
metadata_cleaned$Group     <- relevel(metadata_cleaned$Group, ref = "Ctrl")
metadata_cleaned$Cell_type <- relevel(metadata_cleaned$Cell_type, ref = "monocyte")
stopifnot(all(colnames(filtered_txi$counts) == metadata_cleaned$sample_ID))

# DESeqDataSet, filtrado de genes poco expresados y normalización
dds <- DESeqDataSetFromTximport(txi     = filtered_txi,
                                colData = metadata_cleaned,
                                design  = ~ Cell_type + Group)
dds <- dds[rowSums(counts(dds) >= 10) >= 10, ]
dds <- estimateSizeFactors(dds)

# VST y z-scores por gen
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)
vsd_mat <- vsd_mat[apply(vsd_mat, 1, var) > 0, ]
zscore_mat <- t(scale(t(vsd_mat)))

stopifnot(identical(dim(dds), c(13596L, 99L)), ncol(zscore_mat) == 99)

# dds y zscore_mat no van en data/ por su tamaño: se suben como archivos del
# release v0.1.0 de GitHub y se obtienen con download_data(). data-raw/release/
# está en .gitignore para no subir los archivos al repositorio.
dir.create("data-raw/release", showWarnings = FALSE)
saveRDS(dds,        "data-raw/release/dds.rds",        compress = "xz")
saveRDS(zscore_mat, "data-raw/release/zscore_mat.rds", compress = "xz")
# Luego, en la terminal:
# gh release upload v0.1.0 data-raw/release/dds.rds data-raw/release/zscore_mat.rds --clobber
