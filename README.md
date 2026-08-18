
# GeneExprLupusMx

## Explorando el lupus en México: datos abiertos y análisis reproducibles en R

🎯 El objetivo de `GeneExprLupusMx` es proporcionar un paquete de datos
en R, bilingüe y reproducible, que integre perfiles transcriptómicos y
metadatos clínicos de mujeres mexicanas con **lupus eritematoso
sistémico (SLE) y controles sanos**. Al combinar datos de expresión
génica de *monocitos, células dendríticas derivadas de monocitos (moDCs)
y células dendríticas tolerogénicas (tolDCs)*, el paquete permite:

- **Uso educativo:** servir como recurso didáctico en bioinformática,
  inmunología y ciencias biomédicas.

- **Análisis comparativo:** apoyar la exploración de redes de regulación
  genética y estudios de expresión diferencial.

- **Flujos reproducibles:** ofrecer funciones auxiliares y scripts para
  replicar figuras y análisis del estudio original.

- **Integración global:** hacer visibles y utilizables los conjuntos de
  datos generados localmente en contextos internacionales de
  investigación.

- **Divulgación científica:** promover la ciencia abierta y la difusión
  de datos biomédicos de América Latina.

------------------------------------------------------------------------

## Exploring lupus in Mexico: open data and reproducible analysis in R”

The goal of `GeneExprLupusMx` is to provide a bilingual, reproducible R
data package that integrates transcriptomic profiles and clinical
metadata from Mexican women with **systemic lupus erythematosus (SLE)
and healthy controls**. By combining gene expression data from
*monocytes, monocyte-derived dendritic cells (moDCs), and tolerogenic
dendritic cells (tolDCs)*, the package enables:

- **Educational use:** serving as a teaching resource in bioinformatics,
  immunology, and biomedical sciences.

- **Comparative analysis:** supporting exploration of gene regulation
  networks and differential expression studies.

- **Reproducible workflows:** offering auxiliary functions and scripts
  to replicate figures and analyses from the original study.

- **Global integration:** making locally generated datasets visible and
  usable in international research contexts.

- **Scientific outreach:** promoting open science and the dissemination
  of biomedical data from Latin America.

## Instalación / Installation

Puedes instalar la versión en desarrollo del paquete `GeneExprLupusMx`
desde GitHub con: **/** You can install the development version of
`GeneExprLupusMx` from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("EveliaCoss/GeneExprLupusMx")
```

## Sobre los datos / About the data

Derivado del estudio **“Gene expression profiling of dendritic cell
tolerance dysfunction in women with Systemic lupus erythematosus”**,
este recurso integra datos de expresión génica y metadatos clínicos y
experimentales de mujeres mexicanas con lupus y controles sanas, en los
tres tipos celulares mencionados: *monocitos, células dendríticas
derivadas de monocitos (moDCs) y células dendríticas tolerogénicas
(tolDCs).* **/** Derived from the study **“Gene expression profiling of
dendritic cell tolerance dysfunction in women with Systemic lupus
erythematosus”**, this resource integrates gene expression data and
clinical and experimental metadata from Mexican women with lupus and
healthy controls, across the three mentioned cell types: *monocytes,
monocyte-derived dendritic cells (moDCs), and tolerogenic dendritic
cells (tolDCs).*

El paquete `GeneExprLupusMx` contiene 4 datasets. **/** The
`GeneExprLupusMx` package contains 4 datasets:

- `metadata_cleaned`:
- `metadata`: Archivo suplementario
- `filtered_txi`: Datos sobre expresión diferencial
- `Redes de regulación`:

Cargar datasets **/** Load datasets:

``` r
library(GeneExprLupusMx)
data(package = 'GeneExprLupusMx')
```

One is called `metadata_cleaned`, and is a simplified version of the raw
data; see `?metadata_cleaned` for more info:

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

The second dataset is penguins_raw, and contains all the variables and
original names as downloaded; see ?penguins_raw for more info.

``` r
str(filtered_txi)
#> List of 4
#>  $ abundance          : num [1:21718, 1:99] 7.504 2.357 0 0.228 0.36 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:21718] "A1BG" "A1BG-AS1" "A1CF" "A2M" ...
#>   .. ..$ : chr [1:99] "QR011_0" "QR011_1" "QR011_3" "QR013_0" ...
#>  $ counts             : num [1:21718, 1:99] 179 70 0 16 11 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:21718] "A1BG" "A1BG-AS1" "A1CF" "A2M" ...
#>   .. ..$ : chr [1:99] "QR011_0" "QR011_1" "QR011_3" "QR013_0" ...
#>  $ length             : num [1:21718, 1:99] 1489 1853 8956 4376 1910 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:21718] "A1BG" "A1BG-AS1" "A1CF" "A2M" ...
#>   .. ..$ : chr [1:99] "QR011_0" "QR011_1" "QR011_3" "QR013_0" ...
#>  $ countsFromAbundance: chr "no"
```

Both datasets contain data for 344 penguins. There are 3 different
species of penguins in this dataset, collected from 3 islands in the
Palmer Archipelago, Antarctica.

## Ejemplos / Examples

``` r
library(tidyverse)
#> ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
#> ✔ dplyr     1.2.1     ✔ readr     2.2.0
#> ✔ forcats   1.0.1     ✔ stringr   1.6.0
#> ✔ ggplot2   4.0.3     ✔ tibble    3.3.1
#> ✔ lubridate 1.9.5     ✔ tidyr     1.3.2
#> ✔ purrr     1.2.2     
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
#> ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
metadata_cleaned %>% 
  count(Group, Cell_type)
#>   Group Cell_type  n
#> 1  Ctrl      moDC 10
#> 2  Ctrl  monocyte 10
#> 3  Ctrl     tolDC 10
#> 4   SLE      moDC 23
#> 5   SLE  monocyte 23
#> 6   SLE     tolDC 23
```

Penguins are fun to visualize! For example:

## Cómo citar este paquete / How to cite this package

To cite the palmerpenguins package, please use:

``` r
citation("GeneExprLupusMx")
#> To cite package 'GeneExprLupusMx' in publications use:
#> 
#>   Coss E (2026). _GeneExprLupusMx: GeneExprLupusMx: transcriptómica de
#>   lupus en mujeres mexicanas para investigación y enseñanza_. R package
#>   version 0.0.0.9000, <https://github.com/EveliaCoss/GeneExprLupusMx>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {GeneExprLupusMx: GeneExprLupusMx: transcriptómica de lupus en mujeres mexicanas
#> para investigación y enseñanza},
#>     author = {Evelia Coss},
#>     year = {2026},
#>     note = {R package version 0.0.0.9000},
#>     url = {https://github.com/EveliaCoss/GeneExprLupusMx},
#>   }
```

## Additional data use information

donde vive los datos crudos

## References

Data originally published in:

- Hernández-Ledesma AL, Coss-Navarrete EL, Salazar-Magaña S,
  Ramírez-Espinosa D, Tinajero-Nieto L, Torres-Valdez E, Peña-Ayala AH,
  Félix-Rodriguez G, Frontana-Vázquez G, García Sotelo JS,
  Thomas-Chollier M, Trynka G, Rosetti F, Fernandez-Valverde SL,
  Gutiérrez-Arcelus M, Alpízar-Rodríguez D and Medina-Rivera A (2026)
  Gene expression profiling of dendritic cell tolerance dysfunction in
  women with systemic lupus erythematosus. Front. Immunol. 17:1771959.
  doi: 10.3389/fimmu.2026.1771959
