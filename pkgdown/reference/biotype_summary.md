# Gene biotypes

Biotype of each gene in `filtered_txi`, used to annotate differentially
expressed genes as protein-coding genes (PCG) or long non-coding RNAs
(lncRNA).

## Uso

``` r
biotype_summary
```

## Formato

A data frame with 21,699 rows and 3 columns:

- gene:

  Gene symbol.

- gene_biotype:

  Ensembl biotype (e.g., `"protein_coding"`, `"lncRNA"`).

- Group_biotype:

  Grouped biotype: `PCG`, `lncRNA`, `Pseudogene`, `small ncRNAs` or
  `TEC`.

## Fuente

Ensembl gene annotation.

## Detalles

The script that generates it is `data-raw/lesson4_annotations.R`.

## Ejemplos

``` r
table(biotype_summary$Group_biotype)
#> 
#>          PCG   Pseudogene          TEC       lncRNA small ncRNAs 
#>        17673          483            3         1379         2161 
```
