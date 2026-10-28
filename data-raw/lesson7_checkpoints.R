# Checkpoints de la lección "¿Se diferencian distinto las células de las pacientes?"
# (vignettes/l7_DEGs_between_groups_es.Rmd; Figura 7A y 7B del artículo).
#
# 1. Expresión diferencial con un modelo de interacción en las 99 muestras,
#    igual que en el script original (BQuestion2_v2.qmd, STEP 6,
#    "between-group comparison of contrasts"):
#      diseño ~ Group + Cell_type + Group:Cell_type
#    Los coeficientes de interacción GroupSLE.Cell_typemoDC y
#    GroupSLE.Cell_typetolDC miden cuánto cambia, en las pacientes respecto a
#    los controles, la diferencia célula dendrítica vs monocito.
# 2. Enriquecimiento funcional (g:Profiler; GO:BP, KEGG y REAC) de los genes con
#    padj < 0.05 en cada interacción, separados por dirección (signo del LFC).
#
# Requiere conexión a internet (g:Profiler).
# Ejecutar desde la raíz del paquete: source("data-raw/lesson7_checkpoints.R")

library(DESeq2)
library(gprofiler2)
load("data/filtered_txi.rda")
load("data/metadata_cleaned.rda")

# --- 1. Modelo de interacción ---
dds <- DESeqDataSetFromTximport(txi     = filtered_txi,
                                colData = metadata_cleaned,
                                design  = ~ Group + Cell_type + Group:Cell_type)
dds$Group     <- relevel(dds$Group, ref = "Ctrl")
dds$Cell_type <- relevel(dds$Cell_type, ref = "monocyte")
dds <- DESeq(dds)
print(resultsNames(dds))

de_interaction_results <- lapply(c(moDC = "moDC", tolDC = "tolDC"), function(cell) {
  res <- as.data.frame(results(dds, alpha = 0.05, name = paste0("GroupSLE.Cell_type", cell)))
  data.frame(gene = rownames(res), res, row.names = NULL)
})

# DEG: padj < 0.05 (sin umbral de LFC, como en el script original)
deg_lists <- list()
for (cell in c("moDC", "tolDC")) {
  df <- de_interaction_results[[cell]]
  df <- df[!is.na(df$padj) & df$padj < 0.05, ]
  deg_lists[[paste0("SLEvsCtrl_UP_", cell)]]   <- df$gene[df$log2FoldChange > 0]
  deg_lists[[paste0("SLEvsCtrl_DOWN_", cell)]] <- df$gene[df$log2FoldChange < 0]
}
print(lengths(deg_lists))

# --- 2. Enriquecimiento funcional con g:Profiler ---
# Versión archivada de g:Profiler (junio de 2025)
set_base_url("https://biit.cs.ut.ee/gprofiler_archive3/e113_eg59_p19")

gost_res <- gost(deg_lists, organism = "hsapiens", correction_method = "fdr",
                 user_threshold = 0.05, sources = c("GO:BP", "KEGG", "REAC"), evcodes = TRUE)
r <- gost_res$result
go_interaction_results <- data.frame(
  query             = r$query,
  source            = r$source,
  term_id           = r$term_id,
  term_name         = stringi::stri_trans_general(r$term_name, "Latin-ASCII"),
  term_size         = r$term_size,
  query_size        = r$query_size,
  intersection_size = r$intersection_size,
  p_value           = r$p_value,
  parents           = vapply(r$parents, paste, character(1), collapse = ","),
  intersection      = r$intersection)
attr(go_interaction_results, "gprofiler_version") <- gost_res$meta$version
print(table(go_interaction_results$query, go_interaction_results$source))

usethis::use_data(de_interaction_results, go_interaction_results, overwrite = TRUE, compress = "xz")
