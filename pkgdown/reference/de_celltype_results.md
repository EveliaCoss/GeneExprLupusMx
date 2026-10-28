# Differential expression results: dendritic cells vs monocytes

DESeq2 results comparing in vitro differentiated dendritic cells (moDC
or tolDC) against monocytes, within each group (healthy controls and SLE
patients). Used in the lesson about monocyte differentiation
([`vignette("l5_DEGs_monocytes_diff_es")`](https://eveliacoss.github.io/GeneExprLupusMx/articles/l5_DEGs_monocytes_diff_es.md);
Figures 4 and 5 of the article).

## Uso

``` r
de_celltype_results
```

## Formato

A named list with four data frames (`Ctrl_moDC`, `Ctrl_tolDC`,
`SLE_moDC`, `SLE_tolDC`), each with 21,718 rows (genes) and 7 columns:
`gene`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue` and
`padj` (see
[`de_results`](https://eveliacoss.github.io/GeneExprLupusMx/reference/de_results.md)).

## Fuente

Generated from
[`filtered_txi`](https://eveliacoss.github.io/GeneExprLupusMx/reference/filtered_txi.md)
and
[`metadata_cleaned`](https://eveliacoss.github.io/GeneExprLupusMx/reference/metadata_cleaned.md).

## Detalles

Designs (all genes of `filtered_txi`, reference level `monocyte`):

- Controls: `~ Cell_type`.

- SLE patients: `~ SLEDAI_category + Cell_type + Dose_category`
  (adjusted for disease activity and glucocorticoid dose).

Results were extracted with
`results(alpha = 0.05, name = "Cell_type_<moDC|tolDC>_vs_monocyte")`
(Wald test, Benjamini-Hochberg adjustment).

With `|log2FoldChange| >= 1` and `padj < 0.05` (up / down): Ctrl moDC
2688 / 2947, Ctrl tolDC 2731 / 2849, SLE moDC 2945 / 3652 and SLE tolDC
2804 / 3596. The article (DESeq2 on R 4.0.2) reports the same numbers
for moDC and differs by a few genes in tolDC (2736 / 2851 and 2804 /
3592). The script that generates this object is
`data-raw/lesson5_checkpoints.R`.

## Ver también

[`de_results`](https://eveliacoss.github.io/GeneExprLupusMx/reference/de_results.md)
(SLE vs Ctrl within each cell type).

## Ejemplos

``` r
names(de_celltype_results)
#> [1] "Ctrl_moDC"  "Ctrl_tolDC" "SLE_moDC"   "SLE_tolDC" 
head(de_celltype_results$Ctrl_moDC[order(de_celltype_results$Ctrl_moDC$padj), ])
#>         gene  baseMean log2FoldChange     lfcSE     stat        pvalue
#> 11881   MRC1 68435.137       9.146538 0.2613979 34.99086 3.098523e-268
#> 17371 SLC7A8  4486.754       8.216655 0.2413625 34.04280 5.187183e-254
#> 6460   GPNMB 64903.360       9.443920 0.2816243 33.53376 1.553122e-246
#> 19706  TREM2  6019.286      10.614697 0.3247763 32.68310 2.714885e-234
#> 4641    ECM1  2623.108       7.358515 0.2325977 31.63623 1.172972e-219
#> 890     APOE 70848.544      11.967432 0.3869356 30.92874 4.907850e-210
#>                padj
#> 11881 4.765219e-264
#> 17371 3.988684e-250
#> 6460  7.961823e-243
#> 19706 1.043806e-230
#> 4641  3.607826e-216
#> 890   1.257964e-206
```
