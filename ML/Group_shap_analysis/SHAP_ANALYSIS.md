# Group SHAP, the final set of figures


Group SHAP concept: feature group enters the game **as one player**, and its value is computed by going through all
`2^k` coalitions exactly, with no sampling. This is not the same thing as adding up the importance
of the single features in the group.

Group SHAP values were calculated only for 200 barcodes with 50 background barcodes for each one because of high computational complexity.

**UPD**: `NMD_expected` was removed from the features of both intron datasets, `intron_exp` and
`intron_splicing`. It stays only in `single_combined`.


## Which features are in which group

There are **17 groups and 60 features** in `single_combined`, and **16 groups and 59 features** in
the two intron datasets, which no longer carry `NMD_expected`.

| group | n | features |
|---|---|---|
| `histone_marks` | 13 | H2AFZ, H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K4me2, H3K4me3, H3K79me2, H3K9ac, H3K9me1, H3K9me3, H4K20me1, spKAS-seq |
| `pol2` | 10 | ENCSR000BKR-POLR2AphosphoS5, ENCSR000BMR-POLR2A, ENCSR000DLY-POLR2A, ENCSR000EGF-POLR2AphosphoS2, ENCSR000EHF-POLR2AphosphoS2, ENCSR000EHL-POLR2A, ENCSR000FAJ-POLR2A, ENCSR031TFS-POLR2A, ENCSR388QZF-POLR2A, GSE81553_K562_Smoothed_Pol2_ChIP-Seq |
| `nuclear_org` | 8 | GSE81553_K562_HiC_100kb_Eigenvector, LaminAC_TSA-Seq, LaminB_TSA-Seq, Pol2_TSA-Seq, SON_TSA-Seq, SPIN, Speckle_distance, pSC35_TSA-Seq |
| `gene_annotations` | 7 | DistToTES, DistToTSS, MANE_antisense, MANElength, RelPos, annot_rna_type, gencode_noncoding |
| `transcriptomics` | 5 | GSE81553_K562_GRO-Seq_binned_20kb, GSE81553_K562_RNA-Seq_binned_20kb, RNAseqPolyA, RNAseqPolyAminus, RNAseqPolyAnucleolus |
| `intron_annotations` | 3 | DistToIntron3, DistToIntron5, RelIntronPos |
| `accessibility` | 2 | ATAC, DNAse |
| `gc_content` | 2 | GC, GC_category |
| `global_norm` | 2 | meanMapRPM, meanNormRPM |
| `cage` | 1 | closestCAGEdist |
| `chromhmm` | 1 | ChromHMM |
| `chromosome` | 1 | Chr |
| `conservation` | 1 | phyloP |
| `icshape` | 1 | icSHAPE |
| `mane_expression` | 1 | MANE_all_exp.norm.log10 |
| `nmd` | 1 | NMD_expected. `single_combined` only |
| `repeats` | 1 | RepeatGroup |

### What is dropped for the intergenic table

Outside a MANE gene these 9 features are removed:
`DistToTSS`, `DistToTES`, `RelPos`, `MANElength`, `DistToIntron5`, `DistToIntron3`, `RelIntronPos`,
`MANE_all_exp.norm.log10`, `NMD_expected`.

Three groups disappear completely: `intron_annotations`, `mane_expression` and `nmd`, and in the
intron datasets only the first two, because `nmd` is already absent there. A fourth group,
`gene_annotations`, shrinks from 7 features to 3.

That leaves **14 groups and 51 features**.


## The eight introns

| | | | | | | | |
|---|---|---|---|---|---|---|---|
| PACRGL Control | TTC17 | PACRGL | DDX51 | NBEAL2 | TONSL | USP21 | KPTN |
| `B3-9` | `B3-6` | `H1` | `E1` | `A1` | `B1` | `C1` | `D1` |


## Intragnic and intergenic models

**Intergenic dataset** - integrations  **only outside MANE genes**

**Intergenic dataset** - integrations  **only within MANE genes**

**The splicing classification thresholds** are the tertiles of the pooled sample, **0.288 and
0.692**, which give classes of 33 / 33 / 33 %.

---

# 1. Models quality metrics

**Spearman correlation** for regression and **ROC AUC** for classification. The
`intron_splicing` target has three classes, so its AUC is weighted by class size

An open circle is the training pool (`train`). A filled circle is the 15 % of barcodes hold out
(`test`)

| regression | classification |
|:---:|:---:|
| ![model quality, regression](final_pictures/model_quality_regression.png) | ![model quality, classification](final_pictures/model_quality_classification.png) |

**Regression, Spearman correlation**

| dataset | experiment | intragenic train | intragenic test | intergenic train | intergenic test |
|---|---|---|---|---|---|
| `single_combined` | Regular | 0.917 | 0.833 | 0.767 | 0.710 |
| `single_combined` | PCPA | 0.933 | 0.817 | 0.851 | 0.726 |
| `intron_exp` | PACRGL Control | 0.865 | 0.738 | 0.599 | 0.515 |
| `intron_exp` | TTC17 | 0.807 | 0.711 | 0.523 | 0.421 |
| `intron_exp` | PACRGL | 0.831 | 0.720 | 0.486 | 0.445 |
| `intron_exp` | DDX51 | 0.894 | 0.709 | 0.797 | 0.697 |
| `intron_exp` | NBEAL2 | 0.856 | 0.724 | 0.452 | 0.395 |
| `intron_exp` | TONSL | 0.888 | 0.788 | 0.529 | 0.462 |
| `intron_exp` | USP21 | 0.849 | 0.742 | 0.498 | 0.420 |
| `intron_exp` | KPTN | 0.875 | 0.750 | 0.629 | 0.548 |
| `intron_splicing` | PACRGL Control | 0.614 | 0.363 | 0.590 | 0.382 |
| `intron_splicing` | TTC17 | 0.599 | 0.230 | 0.548 | 0.373 |
| `intron_splicing` | PACRGL | 0.578 | 0.377 | 0.575 | 0.321 |
| `intron_splicing` | DDX51 | 0.624 | 0.322 | 0.565 | 0.253 |
| `intron_splicing` | NBEAL2 | 0.547 | 0.172 | 0.545 | 0.193 |
| `intron_splicing` | TONSL | 0.392 | 0.226 | 0.449 | 0.270 |
| `intron_splicing` | USP21 | 0.682 | 0.270 | 0.510 | 0.310 |
| `intron_splicing` | KPTN | 0.609 | 0.383 | 0.593 | 0.379 |

**Classification, ROC AUC**

| dataset | experiment | intragenic train | intragenic test | intergenic train | intergenic test |
|---|---|---|---|---|---|
| `single_combined` | Regular | 0.973 | 0.938 | 0.885 | 0.873 |
| `single_combined` | PCPA | 0.966 | 0.967 | 0.921 | 0.898 |
| `intron_exp` | PACRGL Control | 0.946 | 0.883 | 0.848 | 0.789 |
| `intron_exp` | TTC17 | 0.929 | 0.898 | 0.794 | 0.736 |
| `intron_exp` | PACRGL | 0.915 | 0.874 | 0.801 | 0.733 |
| `intron_exp` | DDX51 | 0.920 | 0.883 | 0.906 | 0.868 |
| `intron_exp` | NBEAL2 | 0.916 | 0.840 | 0.799 | 0.723 |
| `intron_exp` | TONSL | 0.941 | 0.923 | 0.839 | 0.805 |
| `intron_exp` | USP21 | 0.895 | 0.883 | 0.782 | 0.755 |
| `intron_exp` | KPTN | 0.916 | 0.891 | 0.860 | 0.825 |
| `intron_splicing` | PACRGL Control | 0.885 | 0.703 | 0.865 | 0.700 |
| `intron_splicing` | TTC17 | 0.882 | 0.607 | 0.864 | 0.670 |
| `intron_splicing` | PACRGL | 0.811 | 0.697 | 0.791 | 0.687 |
| `intron_splicing` | DDX51 | 0.775 | 0.656 | 0.808 | 0.656 |
| `intron_splicing` | NBEAL2 | 0.761 | 0.564 | 0.830 | 0.612 |
| `intron_splicing` | TONSL | 0.825 | 0.621 | 0.682 | 0.629 |
| `intron_splicing` | USP21 | 0.758 | 0.618 | 0.733 | 0.639 |
| `intron_splicing` | KPTN | 0.751 | 0.652 | 0.845 | 0.650 |


# 2. Heatmaps of group against experiment

## `single_combined`

| | intragenic | intergenic |
|---|---|---|
| regression | ![](final_pictures/group_shap_heatmap_single_combined_regression_intragenic.png) | ![](final_pictures/group_shap_heatmap_single_combined_regression_intergenic.png) |
| classification | ![](final_pictures/group_shap_heatmap_single_combined_classification_intragenic.png) | ![](final_pictures/group_shap_heatmap_single_combined_classification_intergenic.png) |

## `intron_exp`

| | intragenic | intergenic |
|---|---|---|
| regression | ![](final_pictures/group_shap_heatmap_intron_exp_regression_intragenic.png) | ![](final_pictures/group_shap_heatmap_intron_exp_regression_intergenic.png) |
| classification | ![](final_pictures/group_shap_heatmap_intron_exp_classification_intragenic.png) | ![](final_pictures/group_shap_heatmap_intron_exp_classification_intergenic.png) |

## `intron_splicing`

| | intragenic | intergenic |
|---|---|---|
| regression | ![](final_pictures/group_shap_heatmap_intron_splicing_regression_intragenic.png) | ![](final_pictures/group_shap_heatmap_intron_splicing_regression_intergenic.png) |
| classification | ![](final_pictures/group_shap_heatmap_intron_splicing_classification_intragenic.png) | ![](final_pictures/group_shap_heatmap_intron_splicing_classification_intergenic.png) |

---

# 3. Expression against splicing

The two datasets (intron_exp and intron_splicing) are the same eight introns and the same features. They differ **only in the
target value**. So one point pairs two models that saw the same integrations and the same columns and
were asked two different questions.

Both axes are **logarithmic**. A point is a group,
and the whiskers are the interquartile range of its share across the 8 introns. Regression only.

**With point labels:**

![expression against splicing, regression](final_pictures/group_shap_exp_vs_splicing_regression.png)

**Without point labels** - we can use it to add signatures manually.
![the same, without labels](final_pictures/group_shap_exp_vs_splicing_regression_nolabels.png)

### How the share of a group is computed

SHAP values come out per barcode. Within one experiment we take the mean of the **absolute** value
over barcodes, which gives one number per group. The mean has to be of the absolute value, because
the signed values are centred on the base prediction and would largely cancel each other out. The numbers of all
groups in one model are then normalised to add up to 1, which turns them into shares of that
model's importance budget. The same is done for every experiment, that is for every intron, and the
importance reported for a group is the mean of its share over the experiments.

Two details behind that. The values are not computed on every barcode. Each model explains a random
sample of **200 barcodes** against **50 background donors**. And for a three class
task the mean over barcodes is followed by a mean over the classes, taken after the class margins
are centred.

## The `nmd` feature is a marker of CDS, and must not be read as NMD

This is why the feature was removed from both intron datasets, so the figures above no longer show
an `nmd` group for them. The section stays because it is the reason for the removal.

For the intron datasets `NMD_expected` is **identical to the statement that the integration is
inside a CDS**.
In the earlier run `nmd` led in `intron_splicing` intragenic, which said that the model tells an
insertion into the coding part from an insertion into a UTR, and not anything about NMD.

### How NMD was computed for the intron constructs

The constructs differ between introns, at least in the sequence of the intron itself.

In `single_combined` the construct carried a stop codon only in frames 2 and 3, and not in frame 1.

In the intron datasets the construct always carries a stop codon in all three frames. The insert has
a stop in each of the three frames, and that stop always sits more than 50 nt upstream of the last
exon to exon junction. So NMD is always predicted, and it does not depend on whether the reporter
intron was spliced out.

Two things follow. First, the NMD prediction is always the same for the spliced and the unspliced
form. Second, everything reduces to where the integration landed. In a 5' or 3' UTR, or in the last
exon, there is no NMD. In the CDS there is always NMD, whatever the 50 nt rule would say.

---
# 4. What the models rely on

The share of the importance budget that goes to each group. Every bar is sorted by its own shares,
so a bar can be read on its own, with the strongest group on the left.

![share of the budget per group](final_pictures/group_shap_budget.png)



---

# 5. Spread across experiments

![per experiment, regression intragenic](final_pictures/group_shap_budget_by_experiment_regression_intragenic_fixed_order.png)



