# Gene-level RNA-seq quantification (tximport)

Gene-level RNA-seq quantification of the 99 samples of the study (33
volunteers x 3 cell types: monocytes, moDCs and tolDCs). Transcripts
were quantified with Salmon and summarized to genes with tximport. Genes
are identified by their gene symbol; duplicated entries in gene biotypes
were removed before creating this object.

## Uso

``` r
filtered_txi
```

## Formato

A list with the structure returned by
[`tximport::tximport()`](https://rdrr.io/pkg/tximport/man/tximport.html):

- abundance:

  Numeric matrix (21,718 genes x 99 samples) of gene abundances in
  transcripts per million (TPM).

- counts:

  Numeric matrix (21,718 genes x 99 samples) of estimated read counts.
  Values are not integers because Salmon estimates them.

- length:

  Numeric matrix (21,718 genes x 99 samples) with the average effective
  transcript length of each gene in each sample; DESeq2 uses it as a
  normalization offset.

- countsFromAbundance:

  Character, `"no"`: counts were not derived from abundances.

## Fuente

Salmon quantifications of the study, imported with tximport.

## Detalles

Row names are gene symbols (e.g., `"A1BG"`) and column names are sample
identifiers (e.g., `"QR011_1"`), in the same order as the rows of
[`metadata_cleaned`](https://eveliacoss.github.io/GeneExprLupusMx/reference/metadata_cleaned.md).
Use it with
[`DESeq2::DESeqDataSetFromTximport()`](https://rdrr.io/pkg/DESeq2/man/DESeqDataSet.html).

## Ver también

[`metadata_cleaned`](https://eveliacoss.github.io/GeneExprLupusMx/reference/metadata_cleaned.md),
[`biotype_summary`](https://eveliacoss.github.io/GeneExprLupusMx/reference/biotype_summary.md)

## Ejemplos

``` r
names(filtered_txi)
#> [1] "abundance"           "counts"              "length"             
#> [4] "countsFromAbundance"
dim(filtered_txi$counts)
#> [1] 21718    99
filtered_txi$counts[1:4, 1:4]
#>          QR011_0   QR011_1  QR011_3 QR013_0
#> A1BG     179.000   134.000  401.000 163.000
#> A1BG-AS1  69.977    76.949  205.607  76.942
#> A1CF       0.000     0.000    0.000   0.000
#> A2M       16.000 27777.937 3700.000  29.584
```
