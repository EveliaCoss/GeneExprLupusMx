# Changelog

## GeneExprLupusMx (development version)

### GeneExprLupusMx 0.0.0.9000

#### Lecciones / Lessons

- Ruta de aprendizaje en español con siete lecciones que reproducen las
  figuras del artículo, más una guía de primeros pasos / Spanish
  learning path with seven lessons that reproduce the article’s figures,
  plus a getting-started guide:
  - `intro_es`: primeros pasos con el paquete / getting started.
  - `clinical_es` (Lección 1): estadística descriptiva, Tabla 1 y Figura
    S3.
  - `intro_RNAseq_es` (Lección 2): importación, normalización,
    *heatmaps* y PCA; Figura S2.
  - `l3_DEGs_es` (Lección 3): expresión diferencial SLE vs Ctrl; Figura
    2.
  - `l4_GOterms_es` (Lección 4): enriquecimiento funcional con
    g:Profiler; Figuras 3 y S4.
  - `l5_DEGs_monocytes_diff_es` (Lección 5): diferenciación de monocitos
    a células dendríticas; Figuras 4 y 5.
  - `l6_DEGs_tolDcs_es` (Lección 6): firma transcripcional de las tolDC;
    Figura 6.
  - `l7_DEGs_between_groups_es` (Lección 7): modelo con interacción
    `Group:Cell_type`; Figura 7A y 7B.
- Versión en inglés de todas las lecciones (`*_en`) / English version of
  all lessons (`*_en`).
- Páginas de descarga de datos en español e inglés (`download_es`,
  `download_en`) / Data download pages in Spanish and English.
- Cada lección incluye glosario, mapa de la lección, retos con solución,
  pies de figura con valores calculados en el código y botones de
  navegación por idioma / Each lesson includes a glossary, lesson map,
  challenges with solutions, figure captions with values computed in the
  code, and per-language navigation buttons.

#### Datos / Data

- Nuevos datos de referencia / New reference data: `volunteer_info`,
  `biotype_summary`, `sle_genes`, `isg_genes`, `moDC_genes`, `dc_genes`,
  `tolDC_genes`.
- Nuevos *checkpoints* entre lecciones / New checkpoints between
  lessons: `dds`, `zscore_mat`, `de_results`, `go_results`,
  `de_celltype_results`, `go_celltype_results`,
  `de_interaction_results`, `go_interaction_results`.
- Scripts en `data-raw/` para regenerar cada objeto (g:Profiler
  archivado, versión e113) / Scripts in `data-raw/` to regenerate each
  object (archived g:Profiler, version e113).
- Documentación corregida de `filtered_txi` y `metadata_cleaned` /
  Corrected documentation of `filtered_txi` and `metadata_cleaned`.

#### Funciones / Functions

- [`get_heatmap_mat()`](https://eveliacoss.github.io/GeneExprLupusMx/reference/get_heatmap_mat.md),
  [`create_heatmap()`](https://eveliacoss.github.io/GeneExprLupusMx/reference/create_heatmap.md)
  (rasterización desactivada por defecto / rasterization off by
  default),
  [`processGO()`](https://eveliacoss.github.io/GeneExprLupusMx/reference/processGO.md)
  y
  [`Reduce_GOTerm()`](https://eveliacoss.github.io/GeneExprLupusMx/reference/Reduce_GOTerm.md).

#### Sitio web y documentación / Website and documentation

- Nueva página de inicio (`index.Rmd`) con la imagen *Overview* y la
  tabla de lecciones; README actualizado / New home page (`index.Rmd`)
  with the Overview image and lesson table; updated README.
- Archivo `inst/CITATION`: `citation("GeneExprLupusMx")` devuelve el
  artículo original / `citation("GeneExprLupusMx")` returns the original
  article.
- Menús de pkgdown para lecciones en español e inglés y descarga de
  datos / pkgdown menus for Spanish and English lessons and data
  download.
- Initial version of GeneExprLupusMx.
