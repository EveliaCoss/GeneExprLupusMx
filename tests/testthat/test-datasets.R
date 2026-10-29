# Pruebas de estructura y consistencia de los datos del paquete

de_cols <- c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
go_cols <- c("query", "source", "term_id", "term_name", "term_size", "query_size",
             "intersection_size", "p_value", "parents", "intersection")
sample_ids <- as.character(metadata_cleaned$sample_ID)

expect_de_table <- function(d) {
  expect_s3_class(d, "data.frame")
  expect_named(d, de_cols)
  expect_identical(d$gene, rownames(filtered_txi$counts))
  expect_true(all(d$pvalue >= 0 & d$pvalue <= 1, na.rm = TRUE))
  expect_true(all(d$padj >= 0 & d$padj <= 1, na.rm = TRUE))
  # padj (Benjamini-Hochberg) nunca es menor que el p-valor
  expect_true(all(d$padj >= d$pvalue, na.rm = TRUE))
  expect_true(all(d$baseMean >= 0))
}

expect_go_table <- function(d) {
  expect_s3_class(d, "data.frame")
  expect_named(d, go_cols)
  expect_true(all(d$p_value > 0 & d$p_value <= 1))
  expect_true(all(d$intersection_size <= d$term_size))
  expect_true(all(d$intersection_size <= d$query_size))
  expect_false(anyNA(d$term_id))
}

# ---- Cuentas -----------------------------------------------------------------

test_that("filtered_txi has the tximport structure", {
  expect_type(filtered_txi, "list")
  expect_named(filtered_txi, c("abundance", "counts", "length", "countsFromAbundance"))
  expect_identical(filtered_txi$countsFromAbundance, "no")
  for (m in filtered_txi[c("abundance", "counts", "length")]) {
    expect_true(is.matrix(m))
    expect_identical(dimnames(m), dimnames(filtered_txi$counts))
    expect_false(anyNA(m))
  }
})

test_that("filtered_txi has 21,718 genes and the 99 samples of metadata_cleaned", {
  expect_identical(dim(filtered_txi$counts), c(21718L, 99L))
  expect_identical(colnames(filtered_txi$counts), sample_ids)
  expect_equal(anyDuplicated(rownames(filtered_txi$counts)), 0)
  expect_true(all(filtered_txi$counts >= 0))
  expect_true(all(filtered_txi$abundance >= 0))
  expect_true(all(filtered_txi$length > 0))
})

test_that("dds matches the counts and metadata", {
  skip_if_not_installed("DESeq2")
  expect_s4_class(dds, "DESeqDataSet")
  expect_identical(dim(dds), c(13596L, 99L))
  expect_identical(colnames(dds), sample_ids)
  expect_true(all(rownames(dds) %in% rownames(filtered_txi$counts)))
  expect_identical(deparse(DESeq2::design(dds)), "~Cell_type + Group")
  counts <- DESeq2::counts(dds)
  expect_true(is.integer(counts))
  expect_true(all(counts >= 0))
  expect_true(all(c("Group", "Cell_type", "Volunteer_ID") %in% names(dds@colData)))
})

test_that("zscore_mat has one Z-score per gene and sample of dds", {
  expect_true(is.matrix(zscore_mat))
  expect_identical(dim(zscore_mat), c(13596L, 99L))
  expect_identical(colnames(zscore_mat), sample_ids)
  expect_false(anyNA(zscore_mat))
  # cada gen tiene media 0 y desviación estándar 1
  expect_lt(max(abs(rowMeans(zscore_mat))), 1e-8)
  expect_equal(unname(apply(zscore_mat, 1, stats::sd)), rep(1, nrow(zscore_mat)),
               tolerance = 1e-8)
  skip_if_not_installed("DESeq2")
  expect_identical(rownames(zscore_mat), rownames(dds))
})

# ---- Metadatos ---------------------------------------------------------------

test_that("metadata_cleaned has one row per sample", {
  expect_s3_class(metadata_cleaned, "data.frame")
  expect_identical(dim(metadata_cleaned), c(99L, 14L))
  expect_true(all(c("sample_ID", "Group", "Cell_type", "SLEDAI_category") %in%
                    colnames(metadata_cleaned)))
  expect_equal(anyDuplicated(sample_ids), 0)
  expect_identical(levels(metadata_cleaned$Group), c("Ctrl", "SLE"))
  expect_setequal(levels(metadata_cleaned$Cell_type), c("monocyte", "moDC", "tolDC"))
})

test_that("metadata_cleaned has 10 controls and 23 patients per cell type", {
  tab <- table(metadata_cleaned$Group, metadata_cleaned$Cell_type)
  expect_true(all(tab["Ctrl", ] == 10))
  expect_true(all(tab["SLE", ] == 23))
})

test_that("sample_ID encodes the volunteer and the cell type", {
  volunteer <- sub("_.*", "", sample_ids)
  suffix <- sub(".*_", "", sample_ids)
  expect_true(all(table(volunteer) == 3))
  cell_by_suffix <- c("0" = "monocyte", "1" = "moDC", "3" = "tolDC")
  expect_identical(unname(cell_by_suffix[suffix]), as.character(metadata_cleaned$Cell_type))
})

test_that("volunteer_info has one row per volunteer", {
  expect_s3_class(volunteer_info, "data.frame")
  expect_identical(nrow(volunteer_info), 33L)
  expect_equal(anyDuplicated(volunteer_info$Volunteer_ID), 0)
  expect_identical(as.vector(table(volunteer_info$Group)), c(10L, 23L))
  expect_s3_class(volunteer_info$Time_of_recruitment, "Date")
  expect_true(all(volunteer_info$Age_years > 0))
})

test_that("volunteer_info links to metadata_cleaned with the same group", {
  volunteer <- sub("_.*", "", sample_ids)
  expect_setequal(volunteer_info$Volunteer_ID, volunteer)
  group_meta <- tapply(as.character(metadata_cleaned$Group), volunteer, unique)
  group_info <- stats::setNames(as.character(volunteer_info$Group), volunteer_info$Volunteer_ID)
  expect_identical(as.vector(group_meta[names(group_info)]), unname(group_info))
})

test_that("clinical variables are only filled in for patients", {
  ctrl <- volunteer_info$Group == "Ctrl"
  expect_true(all(is.na(volunteer_info$SLEDAI[ctrl])))
  expect_false(anyNA(volunteer_info$SLEDAI[!ctrl]))
  expect_true(all(volunteer_info$SLEDAI[!ctrl] >= 0))
  expect_true(all(is.na(volunteer_info$Nephritis[ctrl])))
  expect_false(anyNA(volunteer_info$Nephritis[!ctrl]))
})

# ---- Resultados de expresión diferencial -------------------------------------

test_that("de_results has SLE vs Ctrl per cell type", {
  expect_type(de_results, "list")
  expect_named(de_results, c("monocyte", "moDC", "tolDC"))
  for (d in de_results) expect_de_table(d)
})

test_that("de_celltype_results has the differentiation contrasts", {
  expect_type(de_celltype_results, "list")
  expect_named(de_celltype_results, c("Ctrl_moDC", "Ctrl_tolDC", "SLE_moDC", "SLE_tolDC"))
  for (d in de_celltype_results) expect_de_table(d)
})

test_that("de_interaction_results has the group x cell type interaction", {
  expect_type(de_interaction_results, "list")
  expect_named(de_interaction_results, c("moDC", "tolDC"))
  for (d in de_interaction_results) expect_de_table(d)
})

# ---- Resultados de enriquecimiento funcional ---------------------------------

test_that("go_results has the g:Profiler results for SLE vs Ctrl", {
  expect_go_table(go_results)
  expect_setequal(unique(go_results$query),
                  c("SLE_UP_monocyte", "SLE_UP_moDC", "SLE_DOWN_moDC",
                    "SLE_UP_tolDC", "SLE_DOWN_tolDC"))
  expect_true("GO:BP" %in% go_results$source)
  # cada término lista tantos genes como su intersección
  expect_identical(lengths(strsplit(go_results$intersection, ",")),
                   go_results$intersection_size)
})

test_that("go_celltype_results has UP and DOWN for each group and cell type", {
  expect_go_table(go_celltype_results)
  expected <- as.vector(outer(c("Ctrl", "SLE"), c("UP", "DOWN"), paste, sep = "_"))
  expected <- as.vector(outer(expected, c("moDC", "tolDC"), paste, sep = "_"))
  expect_setequal(unique(go_celltype_results$query), expected)
  expect_setequal(unique(go_celltype_results$source), c("GO:BP", "KEGG", "REAC"))
})

test_that("go_interaction_results has UP and DOWN for moDC and tolDC", {
  expect_go_table(go_interaction_results)
  expect_setequal(unique(go_interaction_results$query),
                  c("SLEvsCtrl_UP_moDC", "SLEvsCtrl_DOWN_moDC",
                    "SLEvsCtrl_UP_tolDC", "SLEvsCtrl_DOWN_tolDC"))
  expect_setequal(unique(go_interaction_results$source), c("GO:BP", "KEGG", "REAC"))
})

# ---- Anotaciones y listas de genes -------------------------------------------

test_that("biotype_summary annotates genes of filtered_txi", {
  expect_s3_class(biotype_summary, "data.frame")
  expect_named(biotype_summary, c("gene", "gene_biotype", "Group_biotype"))
  expect_equal(anyDuplicated(biotype_summary$gene), 0)
  expect_true(all(biotype_summary$gene %in% rownames(filtered_txi$counts)))
  expect_setequal(unique(biotype_summary$Group_biotype),
                  c("PCG", "lncRNA", "Pseudogene", "small ncRNAs", "TEC"))
})

test_that("isg_genes is a list of unique genes", {
  expect_s3_class(isg_genes, "data.frame")
  expect_named(isg_genes, "genes")
  expect_identical(nrow(isg_genes), 423L)
  expect_false(anyNA(isg_genes$genes))
  expect_equal(anyDuplicated(isg_genes$genes), 0)
})

test_that("sle_genes lists genes with their reference", {
  expect_s3_class(sle_genes, "data.frame")
  expect_named(sle_genes, c("Gene_name", "Associated_pathways", "Reference", "Link"))
  expect_false(anyNA(sle_genes$Gene_name))
  expect_false(anyNA(sle_genes$Reference))
  expect_true(all(grepl("^https?://", stats::na.omit(sle_genes$Link))))
})

test_that("dc_genes, moDC_genes and tolDC_genes have unique genes", {
  gene_lists <- list(dc_genes = dc_genes, moDC_genes = moDC_genes, tolDC_genes = tolDC_genes)
  for (d in gene_lists) {
    expect_s3_class(d, "data.frame")
    expect_true("gene" %in% names(d))
    expect_false(anyNA(d$gene))
    expect_equal(anyDuplicated(d$gene), 0)
  }
  expect_named(dc_genes, c("gene", "process", "GOterm", "info"))
  expect_setequal(unique(moDC_genes$Expression), c("Upregulated", "Downregulated"))
  expect_true(all(grepl("^(Activated|Repressed) - ", tolDC_genes$tolDC_signature)))
  expect_true(all(grepl("^https?://", c(moDC_genes$link, tolDC_genes$link))))
})
