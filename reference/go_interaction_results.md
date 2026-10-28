# Functional enrichment checkpoint: between-group comparison of contrasts (g:Profiler)

Functional enrichment of the genes whose change from monocytes to
dendritic cells differs between SLE patients and controls
([`de_interaction_results`](https://eveliacoss.github.io/GeneExprLupusMx/reference/de_interaction_results.md)),
obtained with gprofiler2 in
[`vignette("l7_DEGs_between_groups_es")`](https://eveliacoss.github.io/GeneExprLupusMx/articles/l7_DEGs_between_groups_es.md).

## Usage

``` r
go_interaction_results
```

## Format

A data frame with one row per significant term and the same 10 columns
as
[`go_results`](https://eveliacoss.github.io/GeneExprLupusMx/reference/go_results.md):
`query`, `source`, `term_id`, `term_name`, `term_size`, `query_size`,
`intersection_size`, `p_value`, `parents` and `intersection`.

## Source

g:Profiler (<https://biit.cs.ut.ee/gprofiler>).

## Details

Genes with `padj < 0.05` were split by the sign of `log2FoldChange` and
queried as four lists (`SLEvsCtrl_<UP|DOWN>_<moDC|tolDC>`) with
`gost(organism = "hsapiens", correction_method = "fdr", user_threshold = 0.05, sources = c("GO:BP", "KEGG", "REAC"), evcodes = TRUE)`.

The archived g:Profiler version `e113_eg59_p19` (2025) was used; it is
stored in `attr(go_interaction_results, "gprofiler_version")`. The
script that generates it is `data-raw/lesson7_checkpoints.R`.

## See also

[`de_interaction_results`](https://eveliacoss.github.io/GeneExprLupusMx/reference/de_interaction_results.md),
[`go_celltype_results`](https://eveliacoss.github.io/GeneExprLupusMx/reference/go_celltype_results.md)

## Examples

``` r
table(go_interaction_results$query, go_interaction_results$source)
#>                       
#>                        GO:BP KEGG REAC
#>   SLEvsCtrl_DOWN_moDC    228    5   12
#>   SLEvsCtrl_DOWN_tolDC   190    3   14
#>   SLEvsCtrl_UP_moDC      115    0   18
#>   SLEvsCtrl_UP_tolDC      94    2    9
```
