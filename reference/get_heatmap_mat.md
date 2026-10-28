# Generate a heatmap-ready matrix from a DESeq2 object

Applies the variance stabilizing transformation (VST) to a
`DESeqDataSet`, removes genes with zero variance, converts each gene to
Z-scores across samples and orders the columns.

## Usage

``` r
get_heatmap_mat(dds_object, sample_order = colnames(dds_object))
```

## Arguments

- dds_object:

  A `DESeqDataSet` (for example, the samples of a single cell type).

- sample_order:

  Character vector with the sample IDs in the order in which the columns
  should appear. By default, the order of `dds_object`.

## Value

A numeric matrix of Z-scores (genes in rows, samples in columns).

## Details

Requires the Bioconductor packages DESeq2 and SummarizedExperiment.

## Examples

``` r
if (requireNamespace("DESeq2", quietly = TRUE)) {
  mo <- dds[, dds$Cell_type == "monocyte"]
  DESeq2::design(mo) <- ~ Group # a single cell type: compare groups only
  mo_mat <- get_heatmap_mat(mo)
  dim(mo_mat)
}
#> [1] 13569    33
```
