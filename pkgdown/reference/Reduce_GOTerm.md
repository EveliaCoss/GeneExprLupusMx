# Batch Reduce and Visualize GO Terms Across Multiple Datasets

Applies semantic similarity reduction to GO enrichment results from a
batch of `.RData` files (both UP and DOWN regulated genes). For each
dataset, it invokes
[`processGO`](https://eveliacoss.github.io/GeneExprLupusMx/reference/processGO.md)
to reduce redundancy, merge metadata, and generate visualizations.

## Uso

``` r
Reduce_GOTerm(indir, outdir, figdir)
```

## Argumentos

- indir:

  Character. Input directory containing GO enrichment `.RData` files.
  Expected filenames follow the pattern `UP_GO_<name>.RData` and
  `DOWN_GO_<name>.RData`.

- outdir:

  Character. Output directory to save reduced GO term `.RData` files.

- figdir:

  Character. Directory to save output plots (treemaps and bar plots).

## Valor

Invisibly, `NULL`. It generates and saves `.RData` and `.png` files in
the specified output directories.

## Detalles

This wrapper function automates analysis across multiple files following
the naming pattern `UP_GO_<name>.RData` or `DOWN_GO_<name>.RData`,
typically produced by GO enrichment pipelines.

For each file, the function:

- Loads the GO enrichment result.

- Filters invalid, duplicated, or missing GO IDs.

- Skips datasets with fewer than 40 unique GO terms.

- Computes semantic similarity with `rrvgo::calculateSimMatrix()` using
  the "Rel" method.

- Reduces terms via `rrvgo::reduceSimMatrix()` with a similarity
  threshold of 0.7.

- Saves both the reduced term object and merged gene-GO mapping.

- Generates treemap and bar plots showing parent GO categories.

Internal errors in similarity computation are caught with
[`tryCatch()`](https://rdrr.io/r/base/conditions.html) to allow
uninterrupted batch processing across datasets.

Requires the Bioconductor packages rrvgo and org.Hs.eg.db.

## Ver también

[`processGO`](https://eveliacoss.github.io/GeneExprLupusMx/reference/processGO.md)
for analyzing individual files.

## Ejemplos

``` r
if (FALSE) { # \dontrun{
# Batch process all GO enrichment results in the directory
Reduce_GOTerm(
  indir = "data/enrichment/",
  outdir = "results/GO_reduced/",
  figdir = "figures/GO_terms/"
)
} # }
```
