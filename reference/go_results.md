# Functional enrichment results checkpoint (g:Profiler)

Functional enrichment of the differentially expressed genes (DEG) of
each cell type, obtained with gprofiler2 in the lesson "Funciones
biológicas de los genes diferencialmente expresados"
([`vignette("l4_GOterms_es")`](https://eveliacoss.github.io/GeneExprLupusMx/articles/l4_GOterms_es.md)).

## Usage

``` r
go_results
```

## Format

A data frame with 943 rows (significant terms) and 10 columns:

- query:

  Gene list (cell type and direction of change).

- source:

  Database: `GO:BP`, `GO:CC`, `GO:MF`, `REAC`, `KEGG`, `WP`, `TF`,
  `MIRNA`, `HP` or `CORUM`.

- term_id:

  Term identifier.

- term_name:

  Term name.

- term_size:

  Number of genes annotated to the term.

- query_size:

  Number of genes of the list with annotation in the database.

- intersection_size:

  Number of genes of the list annotated to the term.

- p_value:

  Adjusted p-value (FDR) of the hypergeometric test.

- parents:

  Comma-separated parent terms in the database hierarchy.

- intersection:

  Comma-separated genes of the list annotated to the term.

## Source

g:Profiler (<https://biit.cs.ut.ee/gprofiler>).

## Details

DEG were defined in
[`de_results`](https://eveliacoss.github.io/GeneExprLupusMx/reference/de_results.md)
with `|log2FoldChange| >= 1` and `padj < 0.05`, and queried as five
lists (`SLE_UP_monocyte`, `SLE_UP_moDC`, `SLE_DOWN_moDC`,
`SLE_UP_tolDC`, `SLE_DOWN_tolDC`) with
`gost(organism = "hsapiens", correction_method = "fdr", user_threshold = 0.05, evcodes = TRUE)`.

g:Profiler updates its databases periodically. To reproduce the article,
the archived version `e113_eg59_p19` (June 2025) was used; it is stored
in `attr(go_results, "gprofiler_version")`. The script that generates it
is `data-raw/lesson4_checkpoints.R`.

## Examples

``` r
table(go_results$query)
#> 
#>   SLE_DOWN_moDC  SLE_DOWN_tolDC     SLE_UP_moDC SLE_UP_monocyte    SLE_UP_tolDC 
#>             291              39               9             322             282 
table(go_results$source)
#> 
#> CORUM GO:BP GO:CC GO:MF    HP  KEGG MIRNA  REAC    TF    WP 
#>    24   556    15    80    65    26     1    47    78    51 
```
