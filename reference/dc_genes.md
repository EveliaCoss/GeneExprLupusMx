# Genes associated with dendritic cell differentiation and function

Curated list of genes associated with Gene Ontology terms related to
dendritic cells, grouped by biological process. Used to explore
canonical dendritic cell differentiation pathways (Figure 5B of the
article).

## Usage

``` r
dc_genes
```

## Format

A data frame with 94 rows and 4 columns:

- gene:

  Gene symbol.

- process:

  Biological process (e.g., `"DC differentiation"`, `"Migration"`,
  `"TAM signaling"`).

- GOterm:

  Group of Gene Ontology terms used to select the gene.

- info:

  Gene Ontology term associated with the gene.

## Source

Gene Ontology annotations (<https://geneontology.org>).

## Details

The script that generates it is `data-raw/lesson5_data.R`.

## Examples

``` r
table(dc_genes$process)
#> 
#>                  Adhesion molecules Antigen processing and presentation 
#>                                   1                                   7 
#>                           Apoptosis                 Chemokine receptors 
#>                                   6                                   6 
#>                 Cytokine production                  DC differentiation 
#>                                   4                                  13 
#>                         DC-specific                      Growth factors 
#>                                   3                                   1 
#>                  Immune checkpoints                   Immune modulators 
#>                                   2                                   5 
#>                    Lipid metabolism                           Migration 
#>                                   2                                   7 
#>                 Pattern recognition                  Protein processing 
#>                                   7                                   2 
#>                   Signaling kinases                       TAM signaling 
#>                                   5                                  13 
#>               Transcription factors            Transport and Metabolism 
#>                                   9                                   1 
```
