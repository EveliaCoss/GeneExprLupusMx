# DESeqDataSet checkpoint from the RNA-seq lesson

`DESeqDataSet` built in the lesson "De las cuentas de RNA-seq a la
exploración de datos"
([`vignette("intro_RNAseq_es")`](https://eveliacoss.github.io/GeneExprLupusMx/articles/intro_RNAseq_es.md)),
ready to be used as the starting point of later lessons.

## Usage

``` r
dds
```

## Format

A `DESeqDataSet` with 13,596 genes and 99 samples. Its `colData` is
`metadata_cleaned` plus a `Volunteer_ID` column.

## Source

Generated from `filtered_txi` and `metadata_cleaned`.

## Details

Steps:
[`DESeqDataSetFromTximport()`](https://rdrr.io/pkg/DESeq2/man/DESeqDataSet.html)
from `filtered_txi` and `metadata_cleaned` with design
`~ Cell_type + Group` (reference levels `monocyte` and `Ctrl`), removal
of genes with fewer than 10 reads in at least 10 samples, and
`estimateSizeFactors()`.

Requires the DESeq2 package (Bioconductor) to be used. The script that
generates it is `data-raw/lesson2_checkpoints.R`.

## See also

[`zscore_mat`](https://eveliacoss.github.io/GeneExprLupusMx/reference/zscore_mat.md),
[`filtered_txi`](https://eveliacoss.github.io/GeneExprLupusMx/reference/filtered_txi.md),
[`metadata_cleaned`](https://eveliacoss.github.io/GeneExprLupusMx/reference/metadata_cleaned.md)

## Examples

``` r
if (requireNamespace("DESeq2", quietly = TRUE)) {
  library(DESeq2)
  dds
  design(dds)
}
#> Loading required package: S4Vectors
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: GenomicRanges
#> Loading required package: Seqinfo
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
#> ~Cell_type + Group
```
