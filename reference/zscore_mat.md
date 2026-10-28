# Z-score matrix checkpoint from the RNA-seq lesson

Gene-wise z-scores of the variance-stabilized expression (VST) computed
in the lesson "De las cuentas de RNA-seq a la exploración de datos"
([`vignette("intro_RNAseq_es")`](https://eveliacoss.github.io/GeneExprLupusMx/articles/intro_RNAseq_es.md)).
Useful for heatmaps and clustering.

## Usage

``` r
zscore_mat
```

## Format

A numeric matrix with 13,596 genes (rows) and 99 samples (columns).

## Source

Generated from
[`dds`](https://eveliacoss.github.io/GeneExprLupusMx/reference/dds.md).

## Details

Steps: `vst(dds, blind = FALSE)`, removal of genes with zero variance,
and `t(scale(t(...)))` so that each gene has mean 0 and standard
deviation 1 across samples.

The script that generates it is `data-raw/lesson2_checkpoints.R`.

## See also

[`dds`](https://eveliacoss.github.io/GeneExprLupusMx/reference/dds.md)

## Examples

``` r
dim(zscore_mat)
#> [1] 13596    99
zscore_mat[1:4, 1:4]
#>             QR011_0      QR011_1    QR011_3    QR013_0
#> A1BG      1.0171423  0.152116477  2.5101854  0.3869970
#> A1BG-AS1 -0.3773136 -0.258875884  1.4792124 -0.4765105
#> A2M      -1.3482322  1.116066300  0.2047104 -1.2848021
#> A2M-AS1  -0.2822032  0.007389773 -0.5396974  0.7246234
```
