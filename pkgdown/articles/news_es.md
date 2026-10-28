# Registro de cambios

## GeneExprLupusMx (versión en desarrollo)

### GeneExprLupusMx 0.0.0.9000

Versión inicial de GeneExprLupusMx.

#### Lecciones

- Ruta de aprendizaje con siete lecciones que reproducen las figuras del
  artículo, más una guía de primeros pasos, en español (`*_es`) e inglés
  (`*_en`):
  - `intro`: primeros pasos con el paquete.
  - `clinical` (Lección 1): estadística descriptiva, Tabla 1 y Figura
    S3.
  - `intro_RNAseq` (Lección 2): importación, normalización, *heatmaps* y
    PCA; Figura S2.
  - `l3_DEGs` (Lección 3): expresión diferencial SLE vs Ctrl; Figura 2.
  - `l4_GOterms` (Lección 4): enriquecimiento funcional con g:Profiler;
    Figuras 3 y S4, pruebas hipergeométrica y exacta de Fisher del
    enriquecimiento de ISG, y correlación de Spearman entre los ISG y el
    SLEDAI (tablas suplementarias).
  - `l5_DEGs_monocytes_diff` (Lección 5): diferenciación de monocitos a
    células dendríticas; Figuras 4 y 5.
  - `l6_DEGs_tolDcs` (Lección 6): firma transcripcional de las tolDC;
    Figura 6, y genes concordantes y divergentes con la firma de
    Robertson et al.
  - `l7_DEGs_between_groups` (Lección 7): modelo con interacción
    `Group:Cell_type`; Figuras 7A y 7B, y figura suplementaria de
    trayectorias de expresión de los DEG de interacción.
- Páginas de descarga de datos en español e inglés (`download_es`,
  `download_en`).
- Cada lección incluye glosario, mapa de la lección, retos con solución,
  pies de figura y de tabla con valores calculados en el código y
  botones de navegación por idioma.
- Lecciones homogeneizadas: mismo formato de pasos, glosario unificado,
  prerrequisitos con títulos completos y pies en todas las figuras y
  tablas (incluida la Tabla 1).

#### Datos

- Nuevos datos de referencia: `volunteer_info`, `biotype_summary`,
  `sle_genes`, `isg_genes`, `moDC_genes`, `dc_genes`, `tolDC_genes`.
- Nuevos *checkpoints* entre lecciones: `dds`, `zscore_mat`,
  `de_results`, `go_results`, `de_celltype_results`,
  `go_celltype_results`, `de_interaction_results`,
  `go_interaction_results`.
- Scripts en `data-raw/` para regenerar cada objeto (g:Profiler
  archivado, versión e113).
- Documentación corregida de `filtered_txi` y `metadata_cleaned`.

#### Funciones

- [`get_heatmap_mat()`](https://eveliacoss.github.io/GeneExprLupusMx/reference/get_heatmap_mat.md),
  [`create_heatmap()`](https://eveliacoss.github.io/GeneExprLupusMx/reference/create_heatmap.md)
  (rasterización desactivada por defecto),
  [`processGO()`](https://eveliacoss.github.io/GeneExprLupusMx/reference/processGO.md)
  y
  [`Reduce_GOTerm()`](https://eveliacoss.github.io/GeneExprLupusMx/reference/Reduce_GOTerm.md).

#### Sitio web y documentación

- Nueva página de inicio (`index.Rmd`) con la imagen *Overview* y la
  tabla de lecciones.
- README bilingüe (español e inglés).
- Registro de cambios en español (esta página) e inglés (`NEWS.md`).
- Archivo `inst/CITATION`: `citation("GeneExprLupusMx")` devuelve el
  artículo original.
- Menús de pkgdown para las lecciones en español e inglés, la descarga
  de datos y el registro de cambios.
