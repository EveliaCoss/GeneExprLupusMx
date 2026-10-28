# Checkpoints de la lección "Diferenciación de monocitos a células dendríticas"
# (vignettes/l5_DEGs_monocytes_diff_es.Rmd; Figuras 4 y 5 del artículo).
#
# 1. Expresión diferencial de cada tipo de célula dendrítica (moDC, tolDC)
#    contra monocitos, dentro de cada grupo, igual que en los scripts originales
#    (Step4B_DE_Type_Treatment_control.qmd y Step4C_DE_Type_Treatment_SLE.qmd):
#      - Controles: diseño ~ Cell_type.
#      - Pacientes: diseño ~ SLEDAI_category + Cell_type + Dose_category
#        (ajusta por actividad de la enfermedad y dosis de glucocorticoides).
# 2. Enriquecimiento funcional (g:Profiler; GO:BP, KEGG y REAC) de los DEG
#    (|log2FoldChange| >= 1 y padj < 0.05), separados por dirección.
#
# Requiere conexión a internet (g:Profiler).
# Ejecutar desde la raíz del paquete: source("data-raw/lesson5_checkpoints.R")

library(DESeq2)
library(gprofiler2)
load("data/filtered_txi.rda")
load("data/metadata_cleaned.rda")

# --- 1. Expresión diferencial: célula dendrítica vs monocito ---
full_dds <- DESeqDataSetFromTximport(txi     = filtered_txi,
                                     colData = metadata_cleaned,
                                     design  = ~ Cell_type)

run_celltype_de <- function(group, design) {
  dds_group <- full_dds[, full_dds$Group == group]
  dds_group$Cell_type       <- relevel(dds_group$Cell_type, ref = "monocyte")
  dds_group$SLEDAI_category <- relevel(dds_group$SLEDAI_category, ref = "Zero")
  dds_group$Dose_category   <- relevel(dds_group$Dose_category, ref = "0")
  design(dds_group) <- design
  dds_group <- DESeq(dds_group, quiet = TRUE)
  lapply(c(moDC = "moDC", tolDC = "tolDC"), function(cell) {
    res <- as.data.frame(results(dds_group, alpha = 0.05,
                                 name = paste0("Cell_type_", cell, "_vs_monocyte")))
    data.frame(gene = rownames(res), res, row.names = NULL)
  })
}

ctrl <- run_celltype_de("Ctrl", ~ Cell_type)
sle  <- run_celltype_de("SLE",  ~ SLEDAI_category + Cell_type + Dose_category)
de_celltype_results <- list(Ctrl_moDC = ctrl$moDC, Ctrl_tolDC = ctrl$tolDC,
                            SLE_moDC  = sle$moDC,  SLE_tolDC  = sle$tolDC)

# DEG: |log2FoldChange| >= 1 y padj < 0.05
deg_genes <- function(name, direction) {
  df <- de_celltype_results[[name]]
  df <- df[!is.na(df$padj) & df$padj < 0.05, ]
  if (direction == "UP") df$gene[df$log2FoldChange >= 1] else df$gene[df$log2FoldChange <= -1]
}
print(sapply(names(de_celltype_results), function(n) c(up = length(deg_genes(n, "UP")),
                                                       down = length(deg_genes(n, "DOWN")))))
# Artículo (Figura 4B; DESeq2 en R 4.0.2): Ctrl moDC 2688/2947, Ctrl tolDC 2736/2851,
# SLE moDC 2945/3652, SLE tolDC 2804/3592. Versiones recientes de DESeq2 pueden
# diferir en unos pocos genes.

# --- 2. Enriquecimiento funcional con g:Profiler ---
# Versión archivada de g:Profiler (2025), usada en el artículo
set_base_url("https://biit.cs.ut.ee/gprofiler_archive3/e113_eg59_p19")

deg_lists <- list()
for (group in c("Ctrl", "SLE")) {
  for (cell in c("moDC", "tolDC")) {
    for (direction in c("UP", "DOWN")) {
      deg_lists[[paste(group, direction, cell, sep = "_")]] <- deg_genes(paste(group, cell, sep = "_"), direction)
    }
  }
}

# Una consulta por lista: las listas son grandes (~3,000 genes) y el servidor
# archivado puede no responder si se envían todas juntas.
gost_res <- lapply(names(deg_lists), function(query) {
  message("g:Profiler: ", query)
  gost(deg_lists[query], organism = "hsapiens", correction_method = "fdr",
       user_threshold = 0.05, sources = c("GO:BP", "KEGG", "REAC"), evcodes = TRUE)
})

go_celltype_results <- do.call(rbind, lapply(gost_res, function(g) {
  r <- g$result
  data.frame(query             = r$query,
             source            = r$source,
             term_id           = r$term_id,
             term_name         = stringi::stri_trans_general(r$term_name, "Latin-ASCII"),
             term_size         = r$term_size,
             query_size        = r$query_size,
             intersection_size = r$intersection_size,
             p_value           = r$p_value,
             parents           = vapply(r$parents, paste, character(1), collapse = ","),
             intersection      = r$intersection)
}))
attr(go_celltype_results, "gprofiler_version") <- gost_res[[1]]$meta$version
print(table(go_celltype_results$query))
# Artículo: Ctrl UP moDC 1288, DOWN moDC 1827, UP tolDC 979, DOWN tolDC 2001;
# SLE UP moDC 1016, DOWN moDC 2079, UP tolDC 963, DOWN tolDC 2157.

usethis::use_data(de_celltype_results, go_celltype_results, overwrite = TRUE, compress = "xz")
