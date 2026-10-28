# Anotaciones de genes para la lección de funciones biológicas (Figura 3).
# Fuentes (en data-raw/):
#   - biotype.RData: biotipo de cada gen (objeto biotype_summary).
#   - SLE_associated_genes.csv: genes asociados a SLE reunidos de la literatura.
#   - Yu_etal_2025_ISG.xlsx: lista de genes estimulados por interferón (ISG),
#     Yu et al. 2025, Table S3 (https://doi.org/10.6084/m9.figshare.28353289).
# Ejecutar desde la raíz del paquete: source("data-raw/lesson4_annotations.R")

# --- Biotipos ---
load("data-raw/biotype.RData") # biotype_summary, lncRNAs_names, PCG_names
biotype_summary <- biotype_summary[, c("gene", "gene_biotype", "Group_biotype")]
rownames(biotype_summary) <- NULL
stopifnot(!anyDuplicated(biotype_summary$gene),
          setequal(PCG_names, biotype_summary$gene[biotype_summary$Group_biotype == "PCG"]),
          setequal(lncRNAs_names, biotype_summary$gene[biotype_summary$Group_biotype == "lncRNA"]))

# --- Genes asociados a SLE ---
sle_raw <- read.csv("data-raw/SLE_associated_genes.csv", check.names = FALSE)
to_ascii <- function(x) stringi::stri_trans_general(trimws(x), "Latin-ASCII") # Fernández -> Fernandez
sle_genes <- data.frame(
  Gene_name           = trimws(sle_raw$Gene_name),
  Associated_pathways = to_ascii(sle_raw$`Associated pathways`),
  Reference           = to_ascii(sle_raw$Reference),
  Link                = trimws(sle_raw$Link)
)
sle_genes$Associated_pathways[sle_genes$Associated_pathways == ""] <- NA
stopifnot(nrow(sle_genes) == 219, all(grepl("^[ -~]*$", unlist(sle_genes[!is.na(sle_genes)]))))

# --- Genes estimulados por interferón (ISG) ---
isg_genes <- data.frame(genes = readxl::read_excel("data-raw/Yu_etal_2025_ISG.xlsx")[[1]])
stopifnot(nrow(isg_genes) == 423, !anyDuplicated(isg_genes$genes))

usethis::use_data(biotype_summary, sle_genes, isg_genes, overwrite = TRUE, compress = "xz")
