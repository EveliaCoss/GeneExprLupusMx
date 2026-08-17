# Cleaned metadata with patient and cell type annotations

This dataset includes patient-level information along with cell type
annotations for monocytes (`mo`), monocyte-derived dendritic cells
(`moDC`), and tolerogenic dendritic cells (`tolDC`).

## Uso

``` r
metadata_cleaned
```

## Formato

A data frame with multiple rows (samples) and 13 columns:

- sample_ID:

  Unique identifier for each sample.

- Group:

  Disease group (e.g., SLE patient or control).

- GC:

  Glucocorticoid treatment status.

- Dose:

  Numeric dose of treatment.

- Dose_category:

  Categorical dose classification.

- Age:

  Age of the patient.

- Cell_type:

  Annotated cell type (mo, moDC, tolDC).

- Group_age:

  Combined group and age category.

- Nephritis:

  Presence of nephritis (yes/no).

- Other_tx:

  Other treatments received.

- Sledai_score:

  SLEDAI disease activity score.

- Sledai_group:

  Grouped SLEDAI classification.

- SLEDAI_Md:

  Median SLEDAI score.

- SLEDAI_category:

  Categorical SLEDAI classification.

## Fuente

Patient metadata curated and cleaned for downstream analysis.

## Detalles

This dataset is useful for linking transcriptomic data with clinical and
cell type annotations in systemic lupus erythematosus (SLE).

## Ejemplos

``` r
data(metadata_cleaned)
head(metadata_cleaned)
#>         sample_ID Group GC Dose Dose_category Age Cell_type Group_age Nephritis
#> QR011_0   QR011_0  Ctrl No    0             0  25  monocyte    less30      <NA>
#> QR011_1   QR011_1  Ctrl No    0             0  25      moDC    less30      <NA>
#> QR011_3   QR011_3  Ctrl No    0             0  25     tolDC    less30      <NA>
#> QR013_0   QR013_0  Ctrl No    0             0  23  monocyte    less30      <NA>
#> QR013_1   QR013_1  Ctrl No    0             0  23      moDC    less30      <NA>
#> QR013_3   QR013_3  Ctrl No    0             0  23     tolDC    less30      <NA>
#>         Other_tx Sledai_score Sledai_group SLEDAI_Md SLEDAI_category
#> QR011_0     <NA>         <NA>         <NA>      <NA>            Zero
#> QR011_1     <NA>         <NA>         <NA>      <NA>            Zero
#> QR011_3     <NA>         <NA>         <NA>      <NA>            Zero
#> QR013_0     <NA>         <NA>         <NA>      <NA>            Zero
#> QR013_1     <NA>         <NA>         <NA>      <NA>            Zero
#> QR013_3     <NA>         <NA>         <NA>      <NA>            Zero
table(metadata_cleaned$Cell_type)
#> 
#>     moDC monocyte    tolDC 
#>       33       33       33 
```
