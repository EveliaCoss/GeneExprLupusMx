# Sample metadata: group, cell type and clinical information

Metadata of the 99 RNA-seq samples: 33 volunteers (10 healthy controls
and 23 women with systemic lupus erythematosus, SLE), each contributing
three cell types: monocytes (`monocyte`), monocyte-derived dendritic
cells (`moDC`) and tolerogenic dendritic cells (`tolDC`). Clinical
variables are repeated in the three samples of each volunteer.

## Usage

``` r
metadata_cleaned
```

## Format

A data frame with 99 rows (samples) and 14 columns. All columns are
factors.

- sample_ID:

  Sample identifier: volunteer ID plus cell type suffix (`_0` =
  monocyte, `_1` = moDC, `_3` = tolDC; e.g., `"QR011_1"`). Also used as
  row names.

- Group:

  `Ctrl` (healthy control) or `SLE` (patient).

- GC:

  Glucocorticoid treatment: `No` or `Si`.

- Dose:

  Glucocorticoid dose in mg (`0`, `2.5`, `3`, `5` or `8`); `0` for
  controls and untreated patients.

- Dose_category:

  Glucocorticoid treatment coded as `0` (no) or `1` (yes); equivalent to
  `GC`.

- Age:

  Age in years.

- Cell_type:

  Cell type: `monocyte`, `moDC` or `tolDC`.

- Group_age:

  Age group: `less30` (younger than 30 years) or `more30` (30 years or
  older).

- Nephritis:

  Lupus nephritis: `No` or `Si`; `NA` in controls.

- Other_tx:

  Treatment other than glucocorticoids (e.g., antimalarials,
  immunosuppressants): `No` or `Si`; `NA` in controls.

- Sledai_score:

  SLEDAI disease activity score (0 to 62); `NA` in controls.

- Sledai_group:

  SLEDAI activity group: `Without_activity` (score 0), `Low` (1-4),
  `Moderate` (6), `High` (10-19) or `Very_high` (27 or more); `NA` in
  controls.

- SLEDAI_Md:

  Yes/no indicator from the original clinical records (`No`, `Si` or
  `Sí`); `NA` in controls.

- SLEDAI_category:

  Disease activity category: `Zero` (SLEDAI score 0, and all controls),
  `No_activity` (score 1-3) or `Activity` (score 4 or more).

## Source

Clinical records of the study, curated for the analysis.

## Details

- Rows are in the same order as the columns of the count matrices of
  [`filtered_txi`](https://eveliacoss.github.io/GeneExprLupusMx/reference/filtered_txi.md).

- `Age`, `Dose` and `Sledai_score` are stored as factors; convert them
  with `as.numeric(as.character(x))` before doing calculations.

- The volunteer identifier (e.g., `"QR011"`) is the prefix of
  `sample_ID`, `sub("_.*", "", sample_ID)`; it links this table with
  [`volunteer_info`](https://eveliacoss.github.io/GeneExprLupusMx/reference/volunteer_info.md).

- Known inconsistencies of the original records: `SLEDAI_Md` has the
  level `Si` written in two ways (`Si` and `Sí`); in `Sledai_group`, one
  patient with score 0 is classified as `Low` and one with score 10 as
  `High`. Some ages and treatments differ from
  [`volunteer_info`](https://eveliacoss.github.io/GeneExprLupusMx/reference/volunteer_info.md),
  which comes from the recruitment table.

## See also

[`filtered_txi`](https://eveliacoss.github.io/GeneExprLupusMx/reference/filtered_txi.md),
[`volunteer_info`](https://eveliacoss.github.io/GeneExprLupusMx/reference/volunteer_info.md)

## Examples

``` r
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
table(metadata_cleaned$Group, metadata_cleaned$Cell_type)
#>       
#>        moDC monocyte tolDC
#>   Ctrl   10       10    10
#>   SLE    23       23    23

# Age as a number
summary(as.numeric(as.character(metadata_cleaned$Age)))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   19.00   25.00   34.00   32.21   38.00   47.00 
```
