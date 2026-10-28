# Checkpoint de la lección de funciones biológicas (vignettes/l4_GOterms_es.Rmd).
# Enriquecimiento funcional con g:Profiler de los DEG de cada tipo celular.
# g:Profiler actualiza sus bases de datos; para reproducir el artículo se usa la
# versión archivada e113_eg59_p19 (junio de 2025). Requiere conexión a internet.
# Ejecutar desde la raíz del paquete: source("data-raw/lesson4_checkpoints.R")

library(gprofiler2)
load("data/de_results.rda")

# DEG: |log2FoldChange| >= 1 y padj < 0.05 (igual que en la lección 3)
deg_genes <- function(cell, direction) {
  df <- de_results[[cell]]
  df <- df[!is.na(df$padj) & df$padj < 0.05, ]
  if (direction == "UP") df$gene[df$log2FoldChange >= 1] else df$gene[df$log2FoldChange <= -1]
}
deg_lists <- list(SLE_UP_monocyte = deg_genes("monocyte", "UP"),
                  SLE_UP_moDC     = deg_genes("moDC", "UP"),
                  SLE_DOWN_moDC   = deg_genes("moDC", "DOWN"),
                  SLE_UP_tolDC    = deg_genes("tolDC", "UP"),
                  SLE_DOWN_tolDC  = deg_genes("tolDC", "DOWN"))

set_base_url("https://biit.cs.ut.ee/gprofiler_archive3/e113_eg59_p19")
gost_res <- gost(deg_lists, organism = "hsapiens", correction_method = "fdr",
                 user_threshold = 0.05, evcodes = TRUE)

go_results <- data.frame(
  query             = gost_res$result$query,
  source            = gost_res$result$source,
  term_id           = gost_res$result$term_id,
  term_name         = stringi::stri_trans_general(gost_res$result$term_name, "Latin-ASCII"),
  term_size         = gost_res$result$term_size,
  query_size        = gost_res$result$query_size,
  intersection_size = gost_res$result$intersection_size,
  p_value           = gost_res$result$p_value,
  parents           = vapply(gost_res$result$parents, paste, character(1), collapse = ","),
  intersection      = gost_res$result$intersection
)
attr(go_results, "gprofiler_version") <- gost_res$meta$version

# Verificación: número de términos del artículo
n_terms <- table(go_results$query)
stopifnot(identical(as.integer(n_terms[c("SLE_UP_monocyte", "SLE_UP_moDC", "SLE_DOWN_moDC",
                                         "SLE_UP_tolDC", "SLE_DOWN_tolDC")]),
                    c(322L, 9L, 291L, 282L, 39L)))

usethis::use_data(go_results, overwrite = TRUE, compress = "xz")
