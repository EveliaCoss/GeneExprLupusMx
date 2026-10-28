# Lesson 6: The transcriptional signature of tolerogenic dendritic cells

## Welcome

**Tolerogenic dendritic cells (tolDCs)** are a therapeutic promise for
lupus: instead of activating the immune system, they calm it down. But
to use them in patients, one question must be answered: **do tolDCs
obtained from lupus patients acquire the same tolerogenic program as
those from healthy people?** In this lesson we evaluate it with a
published **transcriptional signature**.

**What you need to know beforehand:** you should have completed the
lesson [Differentiation of monocytes into dendritic
cells](https://eveliacoss.github.io/GeneExprLupusMx/articles/l5_DEGs_monocytes_diff_en.md).

**By the end you will be able to:**

1.  Use a **transcriptional signature** from the literature to evaluate
    your own samples.
2.  Combine in one figure the **expression per sample** (z-score) and
    the **change per comparison** (LFC).
3.  Explain why **the choice of reference** changes the interpretation
    of an LFC.
4.  Calculate the change between two cell types from two comparisons
    that share the same reference.
5.  Reproduce **Figure 6** of the article.
6.  Count genes that agree and diverge from a signature and identify
    their common core.

**Estimated time:** 1 hour.

### Lesson map

``` mermaid

flowchart LR
  A["tolDC_genes
signature from Robertson et al."] --> C["71 signature genes
that are DEGs"]
  B["de_celltype_results
(lesson 5)"] --> C
  C --> D["z-score per sample
(zscore_mat)"]
  C --> E["LFC tolDC
vs monocyte"]
  D --> F["Figure 6"]
  E --> F
  C --> G["Does the direction match?
Choosing the reference"]
```

## Glossary

These are the terms we will use. Come back here whenever you need to.

- **Immune tolerance**: the ability of the immune system not to attack
  the body itself; it is lost in lupus.
- **tolDC**: tolerogenic dendritic cells; they promote immune tolerance.
- **Transcriptional signature**: a set of genes whose expression
  characterizes a cell type or state.
- **Activated / repressed**: a gene that goes up / down in the cell
  state of interest relative to a reference.
- **IL-10, dexamethasone (DXM)**: an anti-inflammatory cytokine and a
  glucocorticoid used to generate tolDCs *in vitro*.
- **Reference**: the cell type against which an LFC is calculated;
  changing it changes the sign and size of the LFC.

## 1. Setup

``` r

library(GeneExprLupusMx) # package data
library(tidyverse)       # table manipulation and plots
library(ComplexHeatmap)  # heatmaps
library(circlize)        # color scales for ComplexHeatmap
```

Colors and thresholds (the same as in the previous lessons):

``` r

cell_color    <- c(monocyte = "black", moDC = "palevioletred", tolDC = "gold1")
group_color   <- c(Ctrl = "#96d4ccff", SLE = "#b493b4ff")
biotype_color <- c(PCG = "#7570B3", lncRNA = "#E7298A")
zscore_color  <- colorRamp2(c(-2, -1, 0, 1, 2), c("darkblue", "dodgerblue", "white", "red", "darkred"))
lfc_color     <- colorRamp2(c(-5, -1, 0, 1, 5), c("darkblue", "dodgerblue", "white", "red", "darkred"))

lfc_threshold  <- 1
padj_threshold <- 0.05
```

## 2. **Step 1:** the tolDC signature

[Robertson et al. (2021)](https://doi.org/10.3389/fimmu.2021.733231)
compared tolDCs generated with different protocols and identified the
genes that distinguish tolDCs. The package includes them in
`tolDC_genes`:

``` r

head(tolDC_genes)
#>        gene                        tolDC_signature biotype
#> 1    RNASE1                   Activated - IL10+DXM     PCG
#> 2    S100A8                   Activated - IL10+DXM     PCG
#> 3     CD163 Activated - IL10+DXM, Activated - IL10     PCG
#> 4   SELENOP                   Activated - IL10+DXM     PCG
#> 5 LINC01094                   Activated - IL10+DXM  lncRNA
#> 6   SLC18B1                   Activated - IL10+DXM     PCG
#>                                                                                      link
#> 1 https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.733231/full
#> 2 https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.733231/full
#> 3 https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.733231/full
#> 4 https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.733231/full
#> 5 https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.733231/full
#> 6 https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.733231/full
dplyr::count(tolDC_genes, tolDC_signature)
#>                                                        tolDC_signature  n
#> 1                                                     Activated - IL10  9
#> 2                                                 Activated - IL10+DXM  8
#> 3                               Activated - IL10+DXM, Activated - IL10  1
#> 4                    Activated - IL10+DXM, Activated - tolDC_signature  1
#> 5                                          Activated - tolDC_signature 23
#> 6                                                     Repressed - IL10  7
#> 7                                                 Repressed - IL10+DXM  8
#> 8                                          Repressed - tolDC_signature 25
#> 9                    Repressed - tolDC_signature, Repressed - IL10+DXM  3
#> 10 Repressed - tolDC_signature, Repressed - IL10+DXM, Repressed - IL10  1
```

Each category has two parts:

- **Direction:** `Activated` (the gene goes up in tolDCs) or `Repressed`
  (it goes down).
- **Protocol:** tolDCs generated with IL-10 (`IL10`), with IL-10 and
  dexamethasone (`IL10+DXM`), or the consensus tolDC signature
  (`tolDC_signature`). Some genes appear in several protocols.

## 3. **Step 2:** which signature genes change in our data?

We use the results from the previous lesson (`de_celltype_results`):
each dendritic cell type compared against monocytes, in controls and
patients. We keep the signature genes that are DEGs in **at least one**
of the four comparisons:

``` r

deg_genes <- function(comparison) {
  df <- de_celltype_results[[comparison]]
  df$gene[!is.na(df$padj) & df$padj < padj_threshold & abs(df$log2FoldChange) >= lfc_threshold]
}
all_deg <- unique(unlist(lapply(names(de_celltype_results), deg_genes)))

signature_deg <- tolDC_genes %>%
  filter(gene %in% all_deg, gene %in% rownames(zscore_mat))
nrow(signature_deg)
#> [1] 71
dplyr::count(signature_deg, biotype)
#>   biotype  n
#> 1     PCG 69
#> 2  lncRNA  2
```

**71 of the 86 signature genes change** during the differentiation of
our cells. Almost all of them encode proteins (PCG); two are long
non-coding RNAs (lncRNA).

## 4. **Step 3:** Figure 6

The figure has two *heatmaps* aligned by gene:

- **Left:** the expression (z-score) of each gene in each monocyte and
  tolDC sample, from controls and patients.
- **Right:** the LFC of tolDCs relative to monocytes, in controls and in
  patients.

We prepare the left part with the `zscore_mat` checkpoint:

``` r

sample_info <- metadata_cleaned %>%
  filter(Cell_type %in% c("monocyte", "tolDC")) %>%
  mutate(condition = factor(paste(Group, Cell_type, sep = "_"),
                            levels = c("Ctrl_monocyte", "SLE_monocyte", "Ctrl_tolDC", "SLE_tolDC"))) %>%
  arrange(condition)

zscore_sig <- zscore_mat[signature_deg$gene, sample_info$sample_ID]
dim(zscore_sig)
#> [1] 71 66
```

And the right part, with the tolDC LFC:

``` r

lfc_sig <- sapply(de_celltype_results[c("Ctrl_tolDC", "SLE_tolDC")], function(df)
  df$log2FoldChange[match(signature_deg$gene, df$gene)])
rownames(lfc_sig) <- signature_deg$gene
head(lfc_sig)
#>           Ctrl_tolDC SLE_tolDC
#> RNASE1     11.397672  9.845166
#> S100A8     -7.470788 -7.799939
#> CD163       2.079173  1.746231
#> LINC01094   2.870405  2.776301
#> SLC18B1    -1.023044 -1.479387
#> MERTK       1.427936  1.628772
```

We fix the **order of the categories and of the genes**. Within each
category, the genes are ordered by similarity (*clustering*); here we
write that order explicitly so that the figure always looks the same:

``` r

fig6_order <- list(
  "Repressed - IL10+DXM" = c("MAFF", "SOCS1", "APOL4", "CCL17", "LAMP3", "ALOX15", "CH25H", "ACOT7"),
  "Repressed - tolDC_signature" = c("HCAR2", "FCER1A", "NRGN", "VCL", "ISYNA1", "TRIM36", "ESYT1",
                                    "KCTD6", "PEA15", "NCAPH", "BCAR3", "NDRG2", "CD1A", "CD1E",
                                    "CD1B", "CALCRL", "IRF4", "SLC27A3", "C1orf115", "GALNT12"),
  "Repressed - tolDC_signature, Repressed - IL10+DXM" = c("MMP12", "SCRN1", "CD1C"),
  "Repressed - IL10" = c("HCAR3", "MAP4K1", "TIFAB", "FABP4", "ATP1B1"),
  "Activated - IL10+DXM, Activated - IL10" = c("CD163"),
  "Activated - IL10" = c("HIVEP2", "FPR1", "IL7", "IGF2BP3", "FCGR3A", "FAM20A", "C1S", "TGFA"),
  "Activated - IL10+DXM" = c("SLC18B1", "S100A8", "MERTK", "LINC01094", "C1QB", "RNASE1", "ADAMDEC1"),
  "Repressed - tolDC_signature, Repressed - IL10+DXM, Repressed - IL10" = c("CDH1"),
  "Activated - IL10+DXM, Activated - tolDC_signature" = c("CD14"),
  "Activated - tolDC_signature" = c("NCOA4", "CD300LF", "C5AR1", "CD93", "TREM1", "SMPDL3A", "NRIP1",
                                    "NOD2", "NCF1C", "RGS18", "NINJ1", "CEBPB", "MS4A4A", "PAPSS2",
                                    "C1QA", "SEMA6B", "DRAM1"))

# Are they the same 71 genes we selected?
setequal(unlist(fig6_order), signature_deg$gene)
#> [1] TRUE

# reorder genes, z-scores, and LFC according to the figure
fig6_genes   <- unlist(fig6_order, use.names = FALSE)
fig6_split   <- factor(rep(names(fig6_order), lengths(fig6_order)), levels = names(fig6_order))
signature_deg <- signature_deg[match(fig6_genes, signature_deg$gene), ]
zscore_sig   <- zscore_sig[fig6_genes, ]
lfc_sig      <- lfc_sig[fig6_genes, ]
```

With everything in order, we draw the two panels:

``` r

ht_zscore <- Heatmap(
  zscore_sig, name = "zscore", col = zscore_color,
  cluster_rows = FALSE, cluster_columns = FALSE, # the order comes from fig6_order
  show_row_names = FALSE, show_column_names = FALSE,
  top_annotation = HeatmapAnnotation(
    Group       = as.character(sample_info$Group),
    "Cell type" = as.character(sample_info$Cell_type),
    col = list(Group = group_color, "Cell type" = cell_color),
    show_annotation_name = FALSE),
  column_split = sample_info$condition,
  column_title = c("", "normalized counts", "", ""), column_title_gp = gpar(fontsize = 11),
  row_split = fig6_split, row_title_rot = 0,
  row_title_gp = gpar(fontsize = 9), row_gap = unit(3, "mm"),
  heatmap_legend_param = list(direction = "horizontal"),
  width = unit(12, "cm"))

ht_lfc <- Heatmap(
  lfc_sig, name = "log2FC", col = lfc_color,
  cluster_columns = FALSE, show_column_names = FALSE,
  column_title = "within-group,\ncompared to monocytes", column_title_gp = gpar(fontsize = 11),
  row_names_gp = gpar(fontsize = 8),
  top_annotation = HeatmapAnnotation(
    Group       = c("Ctrl", "SLE"),
    "Cell type" = c("tolDC", "tolDC"),
    col = list(Group = group_color, "Cell type" = cell_color),
    show_legend = FALSE),
  right_annotation = rowAnnotation(Biotype = signature_deg$biotype,
                                   col = list(Biotype = biotype_color)),
  heatmap_legend_param = list(direction = "horizontal"),
  width = unit(2, "cm"))

ht_fig6 <- ht_zscore + ht_lfc

ht_opt$message <- FALSE
draw(ht_fig6, ht_gap = unit(5, "mm"), merge_legends = TRUE,
     heatmap_legend_side = "top", annotation_legend_side = "top")
```

![](l6_DEGs_tolDcs_en_files/figure-html/fig6-1.png)

**Figure 6.** Transcriptional signatures associated with tolDCs obtained
when comparing lupus patients (SLE) and controls. *Heatmaps* showing 71
DEGs associated with the transcriptional signatures of tolDCs (Robertson
et al., 2021): 39 go up and 30 go down in the tolDCs of both groups
relative to monocytes, and 2 change in opposite directions depending on
the group. Rows represent genes, grouped by signature category, and
columns represent monocyte and *in vitro*-differentiated tolDC samples,
colored by cell type and group (SLE or Ctrl). Left panel: *heatmap* with
normalized gene expression values (z-score). Right panel: *heatmap* with
the LFC values of tolDCs relative to monocytes in each group. DEGs were
selected with \|LFC\| ≥ 1 and padj \< 0.05 in at least one comparison
(moDC or tolDC versus monocytes, in patients or in controls). Biotypes
are indicated with colors: purple for protein-coding genes (PCG) and
pink for long non-coding RNAs (lncRNA). SLE, systemic lupus
erythematosus; Ctrl, healthy control; tolDC, tolerogenic dendritic
cells; LFC, log2 fold change.

The original figure can be found in the
[article](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2026.1771959/full).

To save it in high resolution:

``` r

png("Figure_6.png", width = 13, height = 16, units = "in", res = 300)
draw(ht_fig6, ht_gap = unit(5, "mm"), merge_legends = TRUE,
     heatmap_legend_side = "top", annotation_legend_side = "top")
dev.off()
```

**What do you observe?**

- In the left panel, the patterns of the tolDCs from **controls and
  patients are very similar**: the genes that switch on or off do so in
  both groups.
- In the right panel, the LFCs of controls and patients almost always
  have the same color.
- **The patients’ tolDCs acquire the tolerogenic transcriptional
  program**, just like those of the controls. This is good news for
  their potential therapeutic use.

**Challenge 1.** Which signature genes change in **opposite directions**
in the tolDCs of controls and patients?

*Hint:* compare the sign of the two columns of `lfc_sig`.

Show solution

``` r

lfc_sig[sign(lfc_sig[, "Ctrl_tolDC"]) != sign(lfc_sig[, "SLE_tolDC"]), ]
#>          Ctrl_tolDC SLE_tolDC
#> SLC27A3 -0.06446974 0.2778644
#> IL7     -0.57867180 0.2691748
```

Only two genes, and with small LFCs: the tolerogenic program is
practically the same in both groups.

## 5. **Step 4:** do they change in the direction of the signature? The importance of the reference

The signature says which genes are **activated** or **repressed** in
tolDCs. Does it match our LFCs?

``` r

expected <- ifelse(grepl("^Activated", signature_deg$tolDC_signature), 1, -1)
mean(sign(lfc_sig[, "Ctrl_tolDC"]) == expected)
#> [1] 0.4225352
```

Only 42% match! Are the tolDCs poorly differentiated? No: we are using
**a different reference**. Robertson et al. defined the signature by
comparing tolDCs against **non-tolerogenic dendritic cells** (like our
moDCs), whereas our LFCs compare tolDCs against **monocytes**.

Let us look at `CD14`: it is “activated” in the signature, but goes down
in our tolDCs.

``` r

sapply(de_celltype_results, function(df) round(df$log2FoldChange[df$gene == "CD14"], 2))
#>  Ctrl_moDC Ctrl_tolDC   SLE_moDC  SLE_tolDC 
#>      -4.78      -2.23      -5.77      -3.15
```

`CD14` is a monocyte marker: it goes down in tolDCs (−2.2) and goes down
even more in moDCs (−4.8). **Compared with the moDC, the tolDC retains
more `CD14`**, which is exactly what the signature describes.

### 5.1 From two comparisons to one

We do not need to repeat the analysis to compare tolDCs against moDCs.
Since both were compared against the **same monocytes**, in the same
model, we can simply subtract:

``` math
\text{LFC}_{\text{tolDC vs moDC}} = \text{LFC}_{\text{tolDC vs mono}} - \text{LFC}_{\text{moDC vs mono}}
```

``` r

lfc_all <- sapply(de_celltype_results, function(df)
  df$log2FoldChange[match(signature_deg$gene, df$gene)])
rownames(lfc_all) <- signature_deg$gene

tol_vs_modc <- cbind(Ctrl = lfc_all[, "Ctrl_tolDC"] - lfc_all[, "Ctrl_moDC"],
                     SLE  = lfc_all[, "SLE_tolDC"]  - lfc_all[, "SLE_moDC"])

colMeans(sign(tol_vs_modc) == expected) # proportion that matches the signature
#>     Ctrl      SLE 
#> 0.943662 0.915493
```

Now **more than 90% of the genes match** the signature, in controls and
in patients. The tolerogenic differentiation worked; it just had to be
compared against the right reference.

``` r

tibble(gene = rownames(tol_vs_modc), signature = ifelse(expected == 1, "Activated", "Repressed"),
       Ctrl = tol_vs_modc[, "Ctrl"], SLE = tol_vs_modc[, "SLE"]) %>%
  pivot_longer(c(Ctrl, SLE), names_to = "Group", values_to = "LFC") %>%
  ggplot(aes(x = reorder(gene, LFC), y = LFC, fill = signature)) +
  geom_col() +
  facet_wrap(~ Group, ncol = 1) +
  scale_fill_manual(values = c(Activated = "firebrick3", Repressed = "dodgerblue3")) +
  labs(x = NULL, y = "LFC tolDC vs moDC", fill = "Signature") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6), legend.position = "top")
```

![](l6_DEGs_tolDcs_en_files/figure-html/agreement_plot-1.png)

The “activated” genes (red) are almost all on the right, with positive
LFC, and the “repressed” ones (blue) on the left.

**Subtracting LFCs gives magnitude, not significance.** The subtraction
gives the same LFC that `DESeq2` would estimate, but not its p-value. To
statistically test tolDC against moDC, you have to request it from the
`DESeq2` object, for example with
`results(dds, contrast = c("Cell_type", "tolDC", "moDC"))`.

**Challenge 2.** Which genes do **not** match the signature even when
using the moDC as the reference?

Show solution

``` r

tibble(gene = rownames(tol_vs_modc), signature = signature_deg$tolDC_signature,
       Ctrl = round(tol_vs_modc[, "Ctrl"], 2), SLE = round(tol_vs_modc[, "SLE"], 2)) %>%
  filter(sign(Ctrl) != expected | sign(SLE) != expected)
#> # A tibble: 6 × 4
#>   gene    signature                                          Ctrl   SLE
#>   <chr>   <chr>                                             <dbl> <dbl>
#> 1 SOCS1   Repressed - IL10+DXM                               0.5   0.43
#> 2 MMP12   Repressed - tolDC_signature, Repressed - IL10+DXM -0.6   0.91
#> 3 TREM1   Activated - tolDC_signature                        0.27  0   
#> 4 SMPDL3A Activated - tolDC_signature                       -0.71 -0.32
#> 5 RGS18   Activated - tolDC_signature                       -0.75 -0.47
#> 6 CEBPB   Activated - tolDC_signature                       -0.16 -0.35
```

There are few genes, most with LFCs close to zero. Protocol differences
between studies (cytokines, culture times) readily explain this level of
discrepancy.

### 5.2 How many genes agree and how many diverge?

A common way to summarize the comparison with a signature is to
**count**: for each comparison against monocytes, how many signature
DEGs change in the expected direction (**concordant**) and how many in
the opposite one (**divergent**).

``` r

concordance_df <- imap_dfr(de_celltype_results, function(df, comparison) {
  df %>%
    filter(!is.na(padj), padj < padj_threshold, abs(log2FoldChange) >= lfc_threshold) %>%
    inner_join(tolDC_genes, by = "gene") %>%
    mutate(comparison = comparison,
           Group      = sub("_.*", "", comparison),
           status     = ifelse(sign(log2FoldChange) == ifelse(grepl("^Activated", tolDC_signature), 1, -1),
                               "Concordant", "Divergent"))
})

# per comparison
concordance_df %>%
  dplyr::count(comparison, status) %>%
  pivot_wider(names_from = status, values_from = n)
#> # A tibble: 4 × 3
#>   comparison Concordant Divergent
#>   <chr>           <int>     <int>
#> 1 Ctrl_moDC          17        47
#> 2 Ctrl_tolDC         26        31
#> 3 SLE_moDC           17        46
#> 4 SLE_tolDC          26        33

# per group (sum of the moDC and tolDC comparisons)
concordance_counts <- concordance_df %>%
  dplyr::count(Group, status) %>%
  pivot_wider(names_from = status, values_from = n)
concordance_counts
#> # A tibble: 2 × 3
#>   Group Concordant Divergent
#>   <chr>      <int>     <int>
#> 1 Ctrl          43        78
#> 2 SLE           43        79
```

Adding up both comparisons, in controls 43 genes agree with the
signature and 78 diverge; in patients, 43 agree and 79 diverge. Note:
here each gene can be counted twice (once in moDC and once in tolDC).

Which genes agree with the signature in the tolDCs of **both** groups?

``` r

concordant_tolDC <- function(group) {
  concordance_df$gene[concordance_df$comparison == paste0(group, "_tolDC") &
                        concordance_df$status == "Concordant"]
}
core_genes <- sort(intersect(concordant_tolDC("Ctrl"), concordant_tolDC("SLE")))
length(core_genes)
#> [1] 25
core_genes
#>  [1] "ADAMDEC1"  "C1QA"      "C1QB"      "C1S"       "CD163"     "CDH1"     
#>  [7] "DRAM1"     "FAM20A"    "FCER1A"    "FCGR3A"    "HCAR2"     "HCAR3"    
#> [13] "IGF2BP3"   "LINC01094" "MAFF"      "MAP4K1"    "MERTK"     "MS4A4A"   
#> [19] "NRGN"      "PAPSS2"    "RNASE1"    "SCRN1"     "SEMA6B"    "TGFA"     
#> [25] "VCL"
```

These 25 genes form a **tolerogenic core** that the IL-10 protocol
reproduces equally in controls and patients; among them are `CD163`,
`MERTK`, `C1QA`, `C1QB` and `RNASE1`.

But the counts above use **monocytes** as the reference, and section 5.1
showed that with that reference the direction is misleading. Let us
repeat the count with the signature’s reference, tolDC against moDC, for
the 71 signature genes that are DEGs:

``` r

tibble(Group      = colnames(tol_vs_modc),
       Concordant = colSums(sign(tol_vs_modc) == expected),
       Divergent  = colSums(sign(tol_vs_modc) != expected))
#> # A tibble: 2 × 3
#>   Group Concordant Divergent
#>   <chr>      <dbl>     <dbl>
#> 1 Ctrl          67         4
#> 2 SLE           65         6
```

**Results**

- **The IL-10 protocol reproduces a core of the tolerogenic signature:**
  25 genes agree with the signature in the tolDCs of both groups.
- **Agreement is the same in both groups:** 43 concordant genes in
  controls and 43 in patients.
- **Divergence depends on the reference:** against monocytes many genes
  seem to go the opposite way, but against moDCs, the signature’s
  reference, almost all agree.

**Challenge 3.** How many of the 25 tolerogenic core genes are activated
and how many repressed in the signature? Do they still agree with the
signature when the reference is the moDC (`tol_vs_modc`)?

Show solution

``` r

core_expected <- expected[match(core_genes, signature_deg$gene)]
table(ifelse(core_expected == 1, "Activated", "Repressed"))
#> 
#> Activated Repressed 
#>        16         9

colSums(sign(tol_vs_modc[core_genes, ]) == core_expected) # concordant genes
#> Ctrl  SLE 
#>   25   25
```

All 25 genes also agree with the signature against the moDC, in controls
and in patients: they are the most robust core of the tolerogenic
program.

**Challenge 4.** Which gene agrees with the signature in the controls’
tolDCs but not in the patients’? And the other way around?

Show solution

``` r

setdiff(concordant_tolDC("Ctrl"), concordant_tolDC("SLE"))
#> [1] "NDRG2"
setdiff(concordant_tolDC("SLE"), concordant_tolDC("Ctrl"))
#> [1] "C1orf115"
```

Only one gene in each direction: the tolerogenic program of tolDCs is
practically the same in both groups.

## 6. Summary

| Step | Key function | What for? |
|----|----|----|
| Load a signature | `tolDC_genes` | Use prior knowledge to evaluate your data |
| Select genes | `filter(gene %in% all_deg)` | Keep the signature genes that change |
| Expression per sample | `zscore_mat` | See the variability between samples |
| Change per comparison | LFC from `de_celltype_results` | Summarize the change in each group |
| Count agreement | [`imap_dfr()`](https://purrr.tidyverse.org/reference/map_dfr.html), [`intersect()`](https://generics.r-lib.org/reference/setops.html) | Summarize how many genes follow the signature in each group |
| Change the reference | LFC subtraction | Compare two cell types without repeating the analysis |

### Common beginner mistakes

- **Comparing with a signature without checking its reference**: an LFC
  only makes sense relative to something.
- **Concluding that something “failed” because the direction does not
  match**: first check that you are comparing the same thing.
- **Using LFC subtraction as a statistical test**: it gives the
  magnitude, not the p-value.

## 7. Next steps

The patients’ tolDCs acquire the tolerogenic program described in the
literature, very similarly to the controls’. However, comparing the LFCs
of two separate models does not let us statistically test whether
differentiation differs between groups. The lesson [Do patient and
control cells differentiate in the same
way?](https://eveliacoss.github.io/GeneExprLupusMx/articles/l7_DEGs_between_groups_en.md)
answers this question with an interaction model.

To continue, click the **Next lesson: Comparing differentiation between
groups** button at the bottom of this page.

[← Previous lesson: Monocyte
differentiation](https://eveliacoss.github.io/GeneExprLupusMx/articles/l5_DEGs_monocytes_diff_en.md)
[Next lesson: Comparing differentiation between groups
→](https://eveliacoss.github.io/GeneExprLupusMx/articles/l7_DEGs_between_groups_en.md)

[↑ Back to top](#top "Back to the top of the page")

## Session information

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: aarch64-apple-darwin23
#> Running under: macOS Tahoe 26.6.2
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: America/Mexico_City
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] grid      stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] circlize_0.4.18            ComplexHeatmap_2.28.0     
#>  [3] lubridate_1.9.5            forcats_1.0.1             
#>  [5] stringr_1.6.0              dplyr_1.2.1               
#>  [7] purrr_1.2.2                readr_2.2.0               
#>  [9] tidyr_1.3.2                tibble_3.3.1              
#> [11] ggplot2_4.0.3              tidyverse_2.0.0           
#> [13] GeneExprLupusMx_0.0.0.9000
#> 
#> loaded via a namespace (and not attached):
#>  [1] shape_1.4.6.1       gtable_0.3.6        rjson_0.2.23       
#>  [4] xfun_0.60           bslib_0.11.0        htmlwidgets_1.6.4  
#>  [7] GlobalOptions_0.1.4 tzdb_0.5.0          vctrs_0.7.3        
#> [10] tools_4.6.0         generics_0.1.4      stats4_4.6.0       
#> [13] parallel_4.6.0      cluster_2.1.8.2     pkgconfig_2.0.3    
#> [16] RColorBrewer_1.1-3  S7_0.2.2            desc_1.4.3         
#> [19] S4Vectors_0.50.1    lifecycle_1.0.5     compiler_4.6.0     
#> [22] farver_2.1.2        textshaping_1.0.5   codetools_0.2-20   
#> [25] clue_0.3-68         htmltools_0.5.9     sass_0.4.10        
#> [28] yaml_2.3.12         pillar_1.11.1       pkgdown_2.2.1      
#> [31] crayon_1.5.3        jquerylib_0.1.4     cachem_1.1.0       
#> [34] magick_2.9.1        iterators_1.0.14    foreach_1.5.2      
#> [37] tidyselect_1.2.1    digest_0.6.39       stringi_1.8.9      
#> [40] labeling_0.4.3      fastmap_1.2.0       colorspace_2.1-3   
#> [43] cli_3.6.6           magrittr_2.0.5      utf8_1.2.6         
#> [46] withr_3.0.3         scales_1.4.0        timechange_0.4.0   
#> [49] rmarkdown_2.31      matrixStats_1.5.0   otel_0.2.0         
#> [52] ragg_1.5.2          png_0.1-9           GetoptLong_1.1.1   
#> [55] hms_1.1.4           evaluate_1.0.5      knitr_1.51         
#> [58] IRanges_2.46.0      doParallel_1.0.17   rlang_1.3.0        
#> [61] Rcpp_1.1.2          glue_1.8.1          BiocGenerics_0.58.1
#> [64] rstudioapi_0.19.0   jsonlite_2.0.0      R6_2.6.1           
#> [67] systemfonts_1.3.2   fs_2.1.0
```
