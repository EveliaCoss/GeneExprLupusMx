# Checkpoint de la lección "Genes diferencialmente expresados entre pacientes y
# controles" (vignettes/l3_DEGs_es.Rmd). Reproduce el análisis de expresión
# diferencial (SLE vs Ctrl, diseño ~ Group) dentro de cada tipo celular.
# Ejecutar desde la raíz del paquete: source("data-raw/lesson3_checkpoints.R")
# Si cambian los pasos de la lección, actualizar este script y volver a ejecutarlo.

library(DESeq2)
filtered_txi <- GeneExprLupusMx::download_data("filtered_txi") # ver ?download_data
load("data/metadata_cleaned.rda")

metadata_cleaned$Group <- relevel(metadata_cleaned$Group, ref = "Ctrl")
full_dds <- DESeqDataSetFromTximport(txi     = filtered_txi,
                                     colData = metadata_cleaned,
                                     design  = ~ Group)

run_de <- function(cell) {
  dds_cell <- full_dds[, full_dds$Cell_type == cell]
  dds_cell$Cell_type <- droplevels(dds_cell$Cell_type)
  dds_cell <- DESeq(dds_cell, quiet = TRUE)
  res <- as.data.frame(results(dds_cell, name = "Group_SLE_vs_Ctrl"))
  data.frame(gene = rownames(res), res, row.names = NULL)
}

cell_types <- c("monocyte", "moDC", "tolDC")
de_results <- lapply(cell_types, run_de)
names(de_results) <- cell_types

# Verificación: número de DEG de la Figura 2 (|LFC| >= 1 y padj < 0.05)
n_deg <- sapply(de_results, function(df) {
  df <- df[!is.na(df$padj) & df$padj < 0.05, ]
  c(up = sum(df$log2FoldChange >= 1), down = sum(df$log2FoldChange <= -1))
})
stopifnot(identical(unname(n_deg), matrix(c(37L, 0L, 4L, 40L, 4L, 12L), nrow = 2)))

usethis::use_data(de_results, overwrite = TRUE, compress = "xz")
