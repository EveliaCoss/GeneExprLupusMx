# Changelog

## GeneExprLupusMx (development version)

### GeneExprLupusMx 0.0.0.9000

Initial version of GeneExprLupusMx.

#### Lessons

- Learning path with seven lessons that reproduce the article’s figures,
  plus a getting-started guide, in Spanish (`*_es`) and English
  (`*_en`):
  - `intro`: getting started with the package.
  - `clinical` (Lesson 1): descriptive statistics, Table 1 and Figure
    S3.
  - `intro_RNAseq` (Lesson 2): import, normalization, *heatmaps* and
    PCA; Figure S2.
  - `l3_DEGs` (Lesson 3): differential expression SLE vs Ctrl; Figure 2.
  - `l4_GOterms` (Lesson 4): functional enrichment with g:Profiler;
    Figures 3 and S4, hypergeometric and Fisher’s exact tests of ISG
    enrichment, and Spearman correlation between ISGs and SLEDAI
    (supplementary tables).
  - `l5_DEGs_monocytes_diff` (Lesson 5): differentiation of monocytes
    into dendritic cells; Figures 4 and 5.
  - `l6_DEGs_tolDcs` (Lesson 6): transcriptional signature of tolDCs;
    Figure 6, and genes concordant and divergent with the Robertson et
    al. signature.
  - `l7_DEGs_between_groups` (Lesson 7): interaction model
    `Group:Cell_type`; Figures 7A and 7B, and supplementary figure of
    expression trajectories of the interaction DEGs.
- Data download pages in Spanish and English (`download_es`,
  `download_en`).
- Each lesson includes a glossary, lesson map, challenges with
  solutions, figure and table captions with values computed in the code,
  and per-language navigation buttons.
- Lessons harmonized: same step format, unified glossary, prerequisites
  with full titles, and captions for every figure and table (including
  Table 1).

#### Data

- New reference data: `volunteer_info`, `biotype_summary`, `sle_genes`,
  `isg_genes`, `moDC_genes`, `dc_genes`, `tolDC_genes`.
- New checkpoints between lessons: `dds`, `zscore_mat`, `de_results`,
  `go_results`, `de_celltype_results`, `go_celltype_results`,
  `de_interaction_results`, `go_interaction_results`.
- Scripts in `data-raw/` to regenerate each object (archived g:Profiler,
  version e113).
- Corrected documentation of `filtered_txi` and `metadata_cleaned`.

#### Functions

- [`get_heatmap_mat()`](https://eveliacoss.github.io/GeneExprLupusMx/reference/get_heatmap_mat.md),
  [`create_heatmap()`](https://eveliacoss.github.io/GeneExprLupusMx/reference/create_heatmap.md)
  (rasterization off by default),
  [`processGO()`](https://eveliacoss.github.io/GeneExprLupusMx/reference/processGO.md)
  and
  [`Reduce_GOTerm()`](https://eveliacoss.github.io/GeneExprLupusMx/reference/Reduce_GOTerm.md).

#### Website and documentation

- New home page (`index.Rmd`) with the Overview image and lesson table.
- Bilingual README (Spanish and English).
- Changelog in Spanish (`articles/news_es.html`) and English (this
  file).
- `inst/CITATION` file: `citation("GeneExprLupusMx")` returns the
  original article.
- pkgdown menus for Spanish and English lessons, data download and
  changelog.
