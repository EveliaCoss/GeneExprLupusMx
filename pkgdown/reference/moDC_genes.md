# Genes described in the differentiation of monocytes into moDCs

Curated list of genes reported in the literature as up- or
down-regulated when monocytes differentiate in vitro into
monocyte-derived dendritic cells (moDCs), grouped by biological process.
Used to check the expected transcriptional changes (Figure 4C of the
article).

## Uso

``` r
moDC_genes
```

## Formato

A data frame with 90 rows and 7 columns:

- gene:

  Gene symbol.

- Gene_name:

  Name used in the source publication (e.g., `"IL-1beta"`).

- process:

  Biological process (e.g., `"Maturation"`, `"Cytokines"`).

- Associated_pathways:

  Pathway or description reported in the source.

- Expression:

  Reported direction: `Upregulated` or `Downregulated`.

- link:

  Link to the source publication.

- contrast:

  Comparison in which the change was reported.

## Fuente

See the `link` column.

## Detalles

The script that generates it is `data-raw/lesson5_data.R`. Text was
converted to ASCII (Greek letters spelled out) and the spreadsheet
artifact `"GBP 2.00"` was corrected to `"GBP2"`. Eleven entries use
aliases or old names (`CD31`, `CD54`, `IL8`, `IRF7A`, `ERF2`, `CD32`,
`TNFR2`, `IP3KB`, `HSP27`, `GRP78`) or genes absent from the data
(`CXCL8`); they do not match the gene symbols of `filtered_txi` and, as
in the original analysis, are not shown.

## Ejemplos

``` r
table(moDC_genes$process)
#> 
#>                              Adenosine receptors 
#>                                                2 
#>              Antigen processing and presentation 
#>                                                4 
#>                                     Ca2+ binding 
#>                                                2 
#>                       Cell adhesion and motility 
#>                                                6 
#>                                       Chemokines 
#>                                                4 
#>                                        Cytokines 
#>                                               11 
#>                                   Growth factors 
#>                                                2 
#>                                Immune regulation 
#>                                                2 
#>                                 Lipid metabolism 
#>                                                6 
#>                                       Maturation 
#>                                                7 
#>                  Pattern and antigen recognition 
#>                                                3 
#>                                       Regulatory 
#>                                                6 
#>               Signal transduction/Growth control 
#>                                                4 
#>                              Signaling molecules 
#>                                                6 
#> TP53 Regulates Transcription of Cell Death Genes 
#>                                                1 
#>                                     Th2 response 
#>                                               10 
#>                            Transcription factors 
#>                                                3 
#>                         Transcription regulation 
#>                                                7 
#>                                       chaperones 
#>                                                2 
#>                                 oxidative stress 
#>                                                2 
```
