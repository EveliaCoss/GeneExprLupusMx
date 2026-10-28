# Differential expression results: between-group comparison of contrasts

DESeq2 interaction results that test whether the change from monocytes
to in vitro differentiated dendritic cells (moDC or tolDC) is different
in SLE patients than in healthy controls. Used in the lesson about the
between-group comparison of contrasts
([`vignette("l7_DEGs_between_groups_es")`](https://eveliacoss.github.io/GeneExprLupusMx/articles/l7_DEGs_between_groups_es.md);
Figure 7A and 7B of the article).

## Usage

``` r
de_interaction_results
```

## Format

A named list with two data frames (`moDC` and `tolDC`), each with 21,718
rows (genes) and 7 columns: `gene`, `baseMean`, `log2FoldChange`,
`lfcSE`, `stat`, `pvalue` and `padj` (see
[`de_results`](https://eveliacoss.github.io/GeneExprLupusMx/reference/de_results.md)).

## Source

Generated from
[`filtered_txi`](https://eveliacoss.github.io/GeneExprLupusMx/reference/filtered_txi.md)
and
[`metadata_cleaned`](https://eveliacoss.github.io/GeneExprLupusMx/reference/metadata_cleaned.md).

## Details

Design (all genes of `filtered_txi`, the 99 samples, reference levels
`Ctrl` and `monocyte`): `~ Group + Cell_type + Group:Cell_type`. Results
were extracted with
`results(alpha = 0.05, name = "GroupSLE.Cell_type<moDC|tolDC>")` (Wald
test, Benjamini-Hochberg adjustment). The `log2FoldChange` is a
difference of differences: (dendritic cell vs monocyte in SLE) minus
(dendritic cell vs monocyte in Ctrl).

With `padj < 0.05` (up / down): moDC 1 / 28 and tolDC 3 / 38. The script
that generates this object is `data-raw/lesson7_checkpoints.R`.

## See also

[`de_celltype_results`](https://eveliacoss.github.io/GeneExprLupusMx/reference/de_celltype_results.md)
(dendritic cells vs monocytes within each group),
[`go_interaction_results`](https://eveliacoss.github.io/GeneExprLupusMx/reference/go_interaction_results.md).

## Examples

``` r
names(de_interaction_results)
#> [1] "moDC"  "tolDC"
head(de_interaction_results$tolDC[order(de_interaction_results$tolDC$padj), ])
#>         gene   baseMean log2FoldChange     lfcSE      stat       pvalue
#> 7705   ISG15 1571.72540     -2.3242527 0.3964015 -5.863381 4.535367e-09
#> 7348    IFI6 5171.46042     -2.1636287 0.4035187 -5.361905 8.234902e-08
#> 7345   IFI35 1270.40241     -1.1028139 0.2149937 -5.129518 2.904851e-07
#> 1690    BST2 6388.65236     -0.8485252 0.1706817 -4.971388 6.647534e-07
#> 16889 SH3TC2   18.24524     -4.6162394 0.9279896 -4.974452 6.543263e-07
#> 5079    ETV7   57.87252     -3.7694878 0.7671254 -4.913783 8.933559e-07
#>               padj
#> 7705  6.634789e-05
#> 7348  6.023419e-04
#> 7345  1.416502e-03
#> 1690  1.944936e-03
#> 16889 1.944936e-03
#> 5079  2.178151e-03
```
