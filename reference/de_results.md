# Differential expression results checkpoint (SLE vs Ctrl by cell type)

DESeq2 results of the comparison between SLE patients and healthy
controls within each cell type, generated in the lesson "Genes
diferencialmente expresados entre pacientes y controles"
([`vignette("l3_DEGs_es")`](https://eveliacoss.github.io/GeneExprLupusMx/articles/l3_DEGs_es.md)).
Starting point of the functional enrichment lesson.

## Usage

``` r
de_results
```

## Format

A named list with three data frames (`monocyte`, `moDC`, `tolDC`), each
with 21,718 rows (genes) and 7 columns:

- gene:

  Gene symbol.

- baseMean:

  Mean normalized expression across samples.

- log2FoldChange:

  Log2 fold change, SLE vs Ctrl.

- lfcSE:

  Standard error of the log2 fold change.

- stat:

  Wald statistic.

- pvalue:

  Wald test p-value.

- padj:

  Benjamini-Hochberg adjusted p-value (`NA` for genes removed by
  independent filtering).

## Source

Generated from
[`filtered_txi`](https://eveliacoss.github.io/GeneExprLupusMx/reference/filtered_txi.md)
and
[`metadata_cleaned`](https://eveliacoss.github.io/GeneExprLupusMx/reference/metadata_cleaned.md).

## Details

For each cell type:
[`DESeqDataSetFromTximport()`](https://rdrr.io/pkg/DESeq2/man/DESeqDataSet.html)
with all genes of `filtered_txi`, design `~ Group` (reference `Ctrl`),
[`DESeq()`](https://rdrr.io/pkg/DESeq2/man/DESeq.html) and
`results(name = "Group_SLE_vs_Ctrl")` (Wald test, Benjamini-Hochberg
adjustment).

With `|log2FoldChange| >= 1` and `padj < 0.05`: 37 up / 0 down in
monocytes, 4 up / 40 down in moDC and 4 up / 12 down in tolDC (94
distinct genes). The script that generates it is
`data-raw/lesson3_checkpoints.R`.

## Examples

``` r
names(de_results)
#> [1] "monocyte" "moDC"     "tolDC"   
head(de_results$monocyte[order(de_results$monocyte$padj), ])
#>          gene  baseMean log2FoldChange     lfcSE     stat       pvalue
#> 7341    IFI27  745.8090       6.505656 0.9658170 6.735910 1.629076e-11
#> 13612    OTOF  725.5703       5.089491 0.8190471 6.213918 5.167965e-10
#> 16938 SIGLEC1 8360.4454       3.110930 0.5216396 5.963754 2.465075e-09
#> 20398   USP18  758.4596       3.123161 0.5411373 5.771475 7.858055e-09
#> 6301     GMPR  403.5613       1.645094 0.3014701 5.456907 4.844996e-08
#> 7347   IFI44L 6547.0501       2.824219 0.5508756 5.126783 2.947351e-07
#>               padj
#> 7341  1.137584e-07
#> 13612 1.804395e-06
#> 16938 5.737873e-06
#> 20398 1.371820e-05
#> 6301  6.766521e-05
#> 7347  2.725426e-04
```
