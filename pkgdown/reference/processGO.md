# Process GO Enrichment Result from a Single File

Load and analyze a single `.RData` file containing GO enrichment results
from either UP or DOWN regulated genes. Includes validation, semantic
similarity reduction, merging metadata, and plotting.

## Uso

``` r
processGO(f, direction = c("UP", "DOWN"), outdir, figdir)
```

## Argumentos

- f:

  Character. Path to the `.RData` file (must include full or relative
  path). The file must contain exactly one data frame with the columns
  `go_id`, `term` and `p.adjust`.

- direction:

  Character. Either `"UP"` or `"DOWN"` to specify the type of
  enrichment.

- outdir:

  Character. Directory to save `.RData` output files.

- figdir:

  Character. Directory to save PNG plots.

## Valor

Invisibly, `NULL`. Saves processed files and figures in `outdir` and
`figdir`.

## Detalles

Requires the Bioconductor packages rrvgo and org.Hs.eg.db. Datasets with
fewer than 40 unique GO terms are skipped.

## Ver también

[`Reduce_GOTerm`](https://eveliacoss.github.io/GeneExprLupusMx/reference/Reduce_GOTerm.md)
for batch processing of multiple enrichment files.

## Ejemplos

``` r
if (FALSE) { # \dontrun{
# Analyze a specific UP-regulated dataset
processGO(
  f = "data/enrichment/UP_GO_MySample1.RData",
  direction = "UP",
  outdir = "results/GO_terms/",
  figdir = "figures/GO_terms/"
)
} # }
```
