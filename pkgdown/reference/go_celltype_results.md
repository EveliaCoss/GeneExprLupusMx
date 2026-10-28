# Functional enrichment checkpoint: dendritic cells vs monocytes (g:Profiler)

Functional enrichment of the differentially expressed genes (DEG)
between dendritic cells and monocytes
([`de_celltype_results`](https://eveliacoss.github.io/GeneExprLupusMx/reference/de_celltype_results.md)),
obtained with gprofiler2 in the lesson about monocyte differentiation
([`vignette("l5_DEGs_monocytes_diff_es")`](https://eveliacoss.github.io/GeneExprLupusMx/articles/l5_DEGs_monocytes_diff_es.md)).

## Uso

``` r
go_celltype_results
```

## Formato

A data frame with one row per significant term and the same 10 columns
as
[`go_results`](https://eveliacoss.github.io/GeneExprLupusMx/reference/go_results.md):
`query`, `source`, `term_id`, `term_name`, `term_size`, `query_size`,
`intersection_size`, `p_value`, `parents` and `intersection`.

## Fuente

g:Profiler (<https://biit.cs.ut.ee/gprofiler>).

## Detalles

DEG were defined with `|log2FoldChange| >= 1` and `padj < 0.05` and
queried as eight lists (`<Ctrl|SLE>_<UP|DOWN>_<moDC|tolDC>`), one at a
time, with
`gost(organism = "hsapiens", correction_method = "fdr", user_threshold = 0.05, sources = c("GO:BP", "KEGG", "REAC"), evcodes = TRUE)`.

The archived g:Profiler version `e113_eg59_p19` (2025) was used; it is
stored in `attr(go_celltype_results, "gprofiler_version")`. The script
that generates it is `data-raw/lesson5_checkpoints.R`.

## Ver también

[`de_celltype_results`](https://eveliacoss.github.io/GeneExprLupusMx/reference/de_celltype_results.md),
[`go_results`](https://eveliacoss.github.io/GeneExprLupusMx/reference/go_results.md)

## Ejemplos

``` r
table(go_celltype_results$query)
#> 
#>  Ctrl_DOWN_moDC Ctrl_DOWN_tolDC    Ctrl_UP_moDC   Ctrl_UP_tolDC   SLE_DOWN_moDC 
#>            1823            2000            1289             974            2077 
#>  SLE_DOWN_tolDC     SLE_UP_moDC    SLE_UP_tolDC 
#>            2155            1016             963 
```
