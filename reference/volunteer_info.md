# Clinical information per volunteer

Clinical and recruitment information with one row per volunteer (10
healthy controls and 23 women with systemic lupus erythematosus). Each
volunteer contributed three samples (monocyte, moDC and tolDC) to
`metadata_cleaned`; both tables are linked by `Volunteer_ID`, the prefix
of `metadata_cleaned$sample_ID` (e.g. `"QR011_0"` -\> `"QR011"`).

## Usage

``` r
volunteer_info
```

## Format

A data frame with 33 rows and 16 columns:

- Volunteer_ID:

  Volunteer identifier (e.g. `"QR011"`).

- Volunteer:

  Volunteer number, in order of recruitment.

- Group:

  `Ctrl` (healthy control) or `SLE` (patient).

- Age_years:

  Age in years.

- Years_with_diagnosis:

  Years since SLE diagnosis (`NA` in controls).

- Age_at_diagnosis_years:

  Age at SLE diagnosis (`NA` in controls).

- SLEDAI:

  SLEDAI disease activity: `1` = active (SLEDAI \>= 6 points), `0` = no
  or low activity (0-5 points); `NA` in controls.

- Nephritis:

  Lupus nephritis: `No` or `Yes`.

- Treatment_GCC:

  Receives glucocorticoids (`TRUE`/`FALSE`).

- GCC_drug:

  Glucocorticoid as recorded: `PDN` or `PDS`.

- GCC_dose_mg:

  Glucocorticoid dose in mg (`NA` if not recorded).

- Treatment_AM:

  Receives antimalarials.

- Treatment_MMF:

  Receives mycophenolate mofetil.

- Treatment_AZT:

  Receives azathioprine.

- Treatment_MTT:

  Receives methotrexate.

- Time_of_recruitment:

  Month of recruitment and sample collection (`Date`, first day of the
  month).

## Source

Recruitment table of the study.

## Details

Treatment columns are `NA` for controls (not applicable). The cleaning
script is in `data-raw/volunteer_info.R`.

## Examples

``` r
data(volunteer_info)
head(volunteer_info)
#>   Volunteer_ID Volunteer Group Age_years Years_with_diagnosis
#> 1        QR011         1  Ctrl        25                   NA
#> 2        QR013         2  Ctrl        23                   NA
#> 3        QR014         3   SLE        32                    5
#> 4        QR017         4   SLE        22                    4
#> 5        QR020         5   SLE        42                    5
#> 6        QR021         6   SLE        23                    4
#>   Age_at_diagnosis_years SLEDAI Nephritis Treatment_GCC GCC_drug GCC_dose_mg
#> 1                     NA     NA      <NA>            NA     <NA>          NA
#> 2                     NA     NA      <NA>            NA     <NA>          NA
#> 3                     27      1       Yes          TRUE      PDN           5
#> 4                     18      1        No          TRUE      PDN           8
#> 5                     37      0        No         FALSE     <NA>          NA
#> 6                     19      0        No         FALSE     <NA>          NA
#>   Treatment_AM Treatment_MMF Treatment_AZT Treatment_MTT Time_of_recruitment
#> 1           NA            NA            NA            NA          2022-06-01
#> 2           NA            NA            NA            NA          2022-06-01
#> 3         TRUE          TRUE         FALSE         FALSE          2022-06-01
#> 4        FALSE         FALSE         FALSE         FALSE          2022-06-01
#> 5         TRUE         FALSE          TRUE         FALSE          2022-06-01
#> 6        FALSE         FALSE         FALSE         FALSE          2022-06-01
table(volunteer_info$Group, format(volunteer_info$Time_of_recruitment, "%Y-%m"))
#>       
#>        2022-06 2022-07 2022-08 2022-09 2022-10 2022-12 2023-01 2023-02
#>   Ctrl       2       0       3       0       1       1       2       1
#>   SLE        4       1       1       1       4       1      10       1
```
