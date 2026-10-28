# Listas de genes curadas para la lección de diferenciación de monocitos a
# células dendríticas (vignettes/l5_DEGs_monocytes_diff_es.Rmd; Figuras 4C y 5B).
# Fuentes (en data-raw/), copiadas del proyecto original bulkRNA_lupus2024:
#   - Genes_moDC_monocyte_cleaned.csv: genes descritos en la diferenciación de
#     monocito a moDC, por proceso (Figura 4C).
#   - Genes_dendritic_GO.csv: genes asociados a la diferenciación de células
#     dendríticas, por proceso (Figura 5B).
# Ejecutar desde la raíz del paquete: source("data-raw/lesson5_annotations.R")

# --- Genes de la diferenciación monocito -> moDC (Figura 4C) ---
moDC_genes <- read.csv("data-raw/Genes_moDC_monocyte_cleaned.csv")
moDC_genes[] <- lapply(moDC_genes, trimws)
gbp2_artifact <- "GBP\u00a02.00" # "GBP 2.00" con espacio no separable: artefacto de Excel
moDC_genes$gene[moDC_genes$gene == gbp2_artifact] <- "GBP2"
moDC_genes$Gene_name[moDC_genes$Gene_name == gbp2_artifact] <- "GBP2"
moDC_genes$Gene_name <- gsub("γ", "gamma", moDC_genes$Gene_name) # FcγRII -> FcgammaRII
moDC_genes$Gene_name <- gsub("β", "beta", moDC_genes$Gene_name)  # IL-1β -> IL-1beta
moDC_genes$process[moDC_genes$process == "Maduration"] <- "Maturation"
names(moDC_genes)[names(moDC_genes) == "Associated.pathways"] <- "Associated_pathways"
moDC_genes[moDC_genes == ""] <- NA

# --- Genes de la diferenciación de células dendríticas (Figura 5B) ---
dc_genes <- read.csv("data-raw/Genes_dendritic_GO.csv")
dc_genes[] <- lapply(dc_genes, trimws)

ascii <- function(df) all(grepl("^[ -~]*$", na.omit(unlist(df))))
stopifnot(nrow(moDC_genes) == 90, nrow(dc_genes) == 94, ascii(moDC_genes), ascii(dc_genes),
          !anyDuplicated(moDC_genes$gene), !anyDuplicated(dc_genes$gene))

usethis::use_data(moDC_genes, dc_genes, overwrite = TRUE, compress = "xz")
