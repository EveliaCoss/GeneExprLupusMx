# Genes associated with systemic lupus erythematosus (SLE)

Genes reported as associated with SLE in the literature, used to
identify SLE-related genes among the differentially expressed genes.

## Usage

``` r
sle_genes
```

## Format

A data frame with 219 rows (213 unique genes; some genes are reported by
more than one study) and 4 columns:

- Gene_name:

  Gene symbol.

- Associated_pathways:

  Pathway or function associated with the gene (`NA` if not reported).

- Reference:

  Study that reports the association.

- Link:

  Link to the study.

## Source

See the `Link` column.

## Details

Compiled from Shen et al. (2022), Deng and Tsao (2017), Bing et al.
(2016) and Ortiz-Fernandez et al. (2023). The script that generates it
is `data-raw/lesson4_annotations.R`.

## Examples

``` r
table(sle_genes$Reference)
#> 
#> Bing et al, 2016, PLOS ONE        Deng and Tsao, 2017 
#>                          8                        108 
#> Ortiz-Fernandez et al 2023          Shen et al., 2022 
#>                         93                         10 
length(unique(sle_genes$Gene_name))
#> [1] 213
```
