# Transcriptional signature of tolerogenic dendritic cells (tolDCs)

Genes reported as activated (up-regulated) or repressed (down-regulated)
in tolerogenic dendritic cells by Robertson et al. (2021). Used to
evaluate whether the tolDCs of the study acquire the expected
tolerogenic signature (Figure 6 of the article;
`vignette("l6_DEGs_tolDCs_es")`).

## Usage

``` r
tolDC_genes
```

## Format

A data frame with 86 rows and 4 columns:

- gene:

  Gene symbol.

- tolDC_signature:

  Category reported by Robertson et al.: direction (`Activated` or
  `Repressed`) and set in which the gene was found: tolDCs generated
  with IL-10 (`IL10`), with IL-10 plus dexamethasone (`IL10+DXM`) or the
  consensus tolDC signature (`tolDC_signature`). Genes found in several
  sets list all of them, separated by commas.

- biotype:

  Gene biotype (`PCG` or `lncRNA`).

- link:

  Link to the source publication.

## Source

Robertson H, et al. (2021). Transcriptomic analysis identifies a
tolerogenic dendritic cell signature. Front Immunol 12:733231.
[doi:10.3389/fimmu.2021.733231](https://doi.org/10.3389/fimmu.2021.733231)

## Details

The script that generates it is `data-raw/lesson6_annotations.R`.

## Examples

``` r
table(tolDC_genes$tolDC_signature)
#> 
#>                                                    Activated - IL10 
#>                                                                   9 
#>                                                Activated - IL10+DXM 
#>                                                                   8 
#>                              Activated - IL10+DXM, Activated - IL10 
#>                                                                   1 
#>                   Activated - IL10+DXM, Activated - tolDC_signature 
#>                                                                   1 
#>                                         Activated - tolDC_signature 
#>                                                                  23 
#>                                                    Repressed - IL10 
#>                                                                   7 
#>                                                Repressed - IL10+DXM 
#>                                                                   8 
#>                                         Repressed - tolDC_signature 
#>                                                                  25 
#>                   Repressed - tolDC_signature, Repressed - IL10+DXM 
#>                                                                   3 
#> Repressed - tolDC_signature, Repressed - IL10+DXM, Repressed - IL10 
#>                                                                   1 
```
