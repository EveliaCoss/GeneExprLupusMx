# Expectativas propias para las tablas de resultados del paquete, siguiendo
# https://testthat.r-lib.org/articles/custom-expectation.html: cada una registra
# un solo éxito o un solo fallo, con un mensaje que dice qué falló

# Tabla de expresión diferencial (DESeq2) con un renglón por gen de `genes`
expect_de_table <- function(object, genes) {
  act <- quasi_label(rlang::enquo(object))
  d <- act$val
  cols <- c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")

  if (!is.data.frame(d)) {
    fail(c(sprintf("Expected %s to be a data frame.", act$lab),
           sprintf("Actual class: %s.", toString(class(d)))))
  } else if (!identical(names(d), cols)) {
    fail(c(sprintf("Expected %s to have columns %s.", act$lab, toString(cols)),
           sprintf("Actual columns: %s.", toString(names(d)))))
  } else if (!identical(d$gene, genes)) {
    fail(sprintf("Expected %s to have the genes of filtered_txi, in the same order.", act$lab))
  } else if (any(d$pvalue < 0 | d$pvalue > 1 | d$padj < 0 | d$padj > 1, na.rm = TRUE)) {
    fail(sprintf("Expected `pvalue` and `padj` of %s to be between 0 and 1.", act$lab))
  } else if (any(d$padj < d$pvalue, na.rm = TRUE)) {
    # padj (Benjamini-Hochberg) nunca es menor que el p-valor
    fail(sprintf("Expected `padj` of %s to be at least `pvalue`.", act$lab))
  } else if (any(d$baseMean < 0)) {
    fail(sprintf("Expected `baseMean` of %s to be non-negative.", act$lab))
  } else {
    pass()
  }

  invisible(act$val)
}

# Tabla de enriquecimiento funcional (g:Profiler)
expect_go_table <- function(object) {
  act <- quasi_label(rlang::enquo(object))
  d <- act$val
  cols <- c("query", "source", "term_id", "term_name", "term_size", "query_size",
            "intersection_size", "p_value", "parents", "intersection")

  if (!is.data.frame(d)) {
    fail(c(sprintf("Expected %s to be a data frame.", act$lab),
           sprintf("Actual class: %s.", toString(class(d)))))
  } else if (!identical(names(d), cols)) {
    fail(c(sprintf("Expected %s to have columns %s.", act$lab, toString(cols)),
           sprintf("Actual columns: %s.", toString(names(d)))))
  } else if (anyNA(d$term_id)) {
    fail(sprintf("Expected `term_id` of %s to have no missing values.", act$lab))
  } else if (!all(d$p_value > 0 & d$p_value <= 1)) {
    fail(sprintf("Expected `p_value` of %s to be in (0, 1].", act$lab))
  } else if (any(d$intersection_size > pmin(d$term_size, d$query_size))) {
    fail(sprintf("Expected `intersection_size` of %s to be at most `term_size` and `query_size`.",
                 act$lab))
  } else {
    pass()
  }

  invisible(act$val)
}
