# GeneExprLupusMx

The goal of GeneExprLupusMx is to provide a bilingual, reproducible R
data package that integrates transcriptomic profiles and clinical
metadata from Mexican women with systemic lupus erythematosus and
healthy controls. By combining gene expression data from monocytes,
monocyte-derived dendritic cells (moDCs), and tolerogenic dendritic
cells (tolDCs), the package enables:

- **Educational use:** serving as a teaching resource in bioinformatics,
  immunology, and biomedical sciences.

- **Comparative analysis:** supporting exploration of gene regulation
  networks and differential expression studies.

- \*\*Reproducible workflows:\*+ offering auxiliary functions and
  scripts to replicate figures and analyses from the original study.

- **Global integration:** making locally generated datasets visible and
  usable in international research contexts.

- **Scientific outreach:** promoting open science and the dissemination
  of biomedical data from Latin America.

## Installation

You can install the development version of GeneExprLupusMx from
[GitHub](https://github.com/) with:

``` r

# install.packages("pak")
pak::pak("EveliaCoss/GeneExprLupusMx")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r

library(GeneExprLupusMx)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r

summary(metadata_cleaned)
#>    sample_ID   Group     GC      Dose    Dose_category      Age    
#>  QR011_0: 1   Ctrl:30   No:57   0  :57   0:57          19     : 9  
#>  QR011_1: 1   SLE :69   Si:42   2.5: 6   1:42          27     : 9  
#>  QR011_3: 1                     3  : 3                 35     : 9  
#>  QR013_0: 1                     5  :30                 23     : 6  
#>  QR013_1: 1                     8  : 3                 25     : 6  
#>  QR013_3: 1                                            32     : 6  
#>  (Other):93                                            (Other):54  
#>     Cell_type   Group_age  Nephritis Other_tx  Sledai_score
#>  moDC    :33   less30:39   No :60    No :12   0      :18   
#>  monocyte:33   more30:60   Si : 9    Si :57   2      : 6   
#>  tolDC   :33               NAs:30    NAs:30   4      : 6   
#>                                               6      : 6   
#>                                               14     : 6   
#>                                               (Other):27   
#>                                               NAs    :30   
#>            Sledai_group SLEDAI_Md    SLEDAI_category
#>  High            :12    No :30    Activity   :39    
#>  Low             :21    Si :36    No_activity:12    
#>  Moderate        : 6    Sí : 3    Zero       :48    
#>  Very_high       :15    NAs:30                      
#>  Without_activity:15                                
#>  NAs             :30                                
#> 
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.
[`devtools::build_readme()`](https://devtools.r-lib.org/reference/build_readme.html)
is handy for this.
