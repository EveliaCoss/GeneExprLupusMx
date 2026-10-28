# Interferon-stimulated genes (ISGs)

List of 423 interferon-stimulated genes, used to identify differentially
expressed ISGs.

## Usage

``` r
isg_genes
```

## Format

A data frame with 423 rows and 1 column:

- genes:

  Gene symbol.

## Source

Yu et al. (2025), Table S3: "The list of interferon-stimulated genes"
([doi:10.6084/m9.figshare.28353289](https://doi.org/10.6084/m9.figshare.28353289)
).

## Details

The script that generates it is `data-raw/lesson4_annotations.R`.

## Examples

``` r
head(isg_genes)
#>        genes
#> 1      ABTB2
#> 2 AC025171.1
#> 3      ADAM9
#> 4       ADAR
#> 5       AIDA
#> 6     AKR7A2
"IFI27" %in% isg_genes$genes
#> [1] TRUE
```
