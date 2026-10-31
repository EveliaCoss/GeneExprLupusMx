# GeneExprLupusMx (development version)

## Smaller package

- The three largest objects (`filtered_txi`, `dds` and `zscore_mat`) are no longer installed with the package. They are files of the v0.1.0 GitHub release, and the new function `download_data()` downloads them the first time they are used and then reads the saved copy. The package now requires R >= 4.0.
- Only the "Get started" guides (`intro_es`, `intro_en`) are installed as vignettes. The lessons, the functions guide and the data download page are now only on the [website](https://eveliacoss.github.io/GeneExprLupusMx/).
- Smaller figures, and the unused logo and `overview.png` were removed.

## Functions

- New `download_data()` to get `filtered_txi`, `dds` and `zscore_mat`.
- `processGO()` and `Reduce_GOTerm()` were renamed `process_go()` and `reduce_go_terms()`, so that all the functions follow the same style (`snake_case`). The old names still work, but are deprecated and show a warning.

## Website

- New "Shinyapp" section in the README and the home page sidebar, with a link to the [GeneExprLupusMx Shinyapp](https://github.com/EveliaCoss/GeneExprLupusMx_shinyapp) (for now, only in Spanish).

# GeneExprLupusMx 0.1.0

First public release of GeneExprLupusMx, archived in Zenodo.

## Lessons

- Learning path with seven lessons that reproduce the article's figures, plus a getting-started guide, in Spanish (`*_es`) and English (`*_en`):
  - `intro`: getting started with the package.
  - `clinical` (Lesson 1): descriptive statistics, Table 1 and Figure S3.
  - `intro_RNAseq` (Lesson 2): import, normalization, *heatmaps* and PCA; Figure S2.
  - `l3_DEGs` (Lesson 3): differential expression SLE vs Ctrl; Figure 2.
  - `l4_GOterms` (Lesson 4): functional enrichment with g:Profiler; Figures 3 and S4, hypergeometric and Fisher's exact tests of ISG enrichment, and Spearman correlation between ISGs and SLEDAI (supplementary tables).
  - `l5_DEGs_monocytes_diff` (Lesson 5): differentiation of monocytes into dendritic cells; Figures 4 and 5.
  - `l6_DEGs_tolDcs` (Lesson 6): transcriptional signature of tolDCs; Figure 6, and genes concordant and divergent with the Robertson et al. signature.
  - `l7_DEGs_between_groups` (Lesson 7): interaction model `Group:Cell_type`; Figures 7A and 7B, and supplementary figure of expression trajectories of the interaction DEGs.
- Guide on how to use the package functions (`functions_es`, `functions_en`), with real examples of `get_heatmap_mat()`, `create_heatmap()`, `processGO()` and `Reduce_GOTerm()`.
- Lesson 7 ends with an illustrative figure summarizing the findings of the learning path.
- Data download pages in Spanish and English (`download_es`, `download_en`).
- Each lesson opens with an illustrative summary image (after the glossary), in Spanish and English: input data, the lesson steps with its real figures, key findings and checkpoints.
- Each lesson includes a glossary, lesson map, challenges with solutions, figure and table captions with values computed in the code, and per-language navigation buttons.
- Lessons harmonized: same step format, unified glossary, prerequisites with full titles, and captions for every figure and table (including Table 1).

## Data

- New reference data: `volunteer_info`, `biotype_summary`, `sle_genes`, `isg_genes`, `moDC_genes`, `dc_genes`, `tolDC_genes`.
- New checkpoints between lessons: `dds`, `zscore_mat`, `de_results`, `go_results`, `de_celltype_results`, `go_celltype_results`, `de_interaction_results`, `go_interaction_results`.
- Scripts in `data-raw/` to regenerate each object (archived g:Profiler, version e113).
- Corrected documentation of `filtered_txi` and `metadata_cleaned`.

## Functions

- `get_heatmap_mat()`, `create_heatmap()` (rasterization off by default), `processGO()` and `Reduce_GOTerm()`.

## Website and documentation

- New home page (`index.Rmd`) with the Overview image and lesson table.
- Bilingual README (Spanish and English).
- Changelog in Spanish (`articles/news_es.html`) and English (this file).
- `inst/CITATION` file: `citation("GeneExprLupusMx")` returns the package citation; the data citations (article and Zenodo dataset) are listed in the citation sections of the website and README. The package Zenodo DOI is pending.
- pkgdown menus: Getting started, Tutorials (Spanish and English lessons), Functions, data download and changelog; navigation buttons with white text.
