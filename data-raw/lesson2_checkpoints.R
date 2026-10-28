# Checkpoints de la lección "De las cuentas de RNA-seq a la exploración de datos"
# (vignettes/intro_RNAseq_es.Rmd). Reproduce los mismos pasos de la lección para
# que las lecciones siguientes puedan empezar con `dds` y `zscore_mat` ya listos.
# Ejecutar desde la raíz del paquete: source("data-raw/lesson2_checkpoints.R")
# Si cambian los pasos de la lección, actualizar este script y volver a ejecutarlo.

library(DESeq2)
load("data/filtered_txi.rda")
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

usethis::use_data(dds, zscore_mat, overwrite = TRUE, compress = "xz")
