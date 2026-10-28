# Firma transcripcional de células dendríticas tolerogénicas (tolDC) para la
# lección vignettes/l6_DEGs_tolDCs_es.Rmd (Figura 6 del artículo).
# Fuente (en data-raw/), copiada del proyecto original bulkRNA_lupus2024:
#   - Genes_tolDC_siganture_cleaned.csv: genes activados o reprimidos en tolDC
#     según Robertson et al. (2021), Front Immunol 12:733231
#     (https://doi.org/10.3389/fimmu.2021.733231).
# Ejecutar desde la raíz del paquete: source("data-raw/lesson6_annotations.R")

raw <- read.csv("data-raw/Genes_tolDC_siganture_cleaned.csv")

tolDC_genes <- data.frame(
  gene            = trimws(raw$gene),
  tolDC_signature = trimws(raw$tolDC_signature),
  biotype         = trimws(raw$biotype),
  link            = trimws(raw$link)
) # las columnas X e info del archivo original están vacías

stopifnot(nrow(tolDC_genes) == 86, !anyDuplicated(tolDC_genes$gene),
          all(grepl("^(Activated|Repressed) - ", tolDC_genes$tolDC_signature)),
          all(grepl("^[ -~]*$", unlist(tolDC_genes))))

usethis::use_data(tolDC_genes, overwrite = TRUE, compress = "xz")
