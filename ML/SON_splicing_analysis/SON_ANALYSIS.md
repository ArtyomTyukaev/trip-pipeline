# Splicing change when speckles are disrupted

This proved to be quite a challenging and time-consuming task, as the base model configuration we had previously used did not work, and the new data was significantly noisier. Consequently, I had to determine appropriate noise-filtering thresholds and modify certain model hyperparameters.

# Final pictures

DDX51 and NBEAL2, inside MANE genes.

## Model quality

![model quality](final_pictures/son_final_quality.png)

| construct | classification train | classification test | regression train | regression test |
|---|---|---|---|---|
| DDX51 | 0.841 | 0.721 | 0.783 | 0.382 |
| NBEAL2 | 0.837 | 0.822 | 0.758 | 0.367 |

## Is the AUC real

![AUC against chance](final_pictures/son_final_auc_significance.png)

## Group SHAP

| | inside MANE genes | outside MANE genes |
|---|---|---|
| classification | ![](final_pictures/group_shap_heatmap_final_son_cls_classification_intragenic.png) | ![](final_pictures/group_shap_heatmap_final_son_cls_classification_intergenic.png) |
| regression | ![](final_pictures/group_shap_heatmap_final_son_reg_regression_intragenic.png) | ![](final_pictures/group_shap_heatmap_final_son_reg_regression_intergenic.png) |

## Volcano

![volcano](final_pictures/son_volcano_all.png)

![volcano, adjusted P](final_pictures/son_volcano_all_padj.png)

---

## Thresholds used

**Significance, both tasks:** `minPvalue < 0.005`, that is a loose threshold.

### What lfcSE is, and why the threshold is not the same everywhere

When I picked thresholds on the DESeq2 output more or less at random, the models did not learn at
all. Quality was zero, and there was nothing readable in the interpretation either - everything
just came out white.

So the thresholds had to be picked adaptively.

The value I picked to filter on is `lfcSE`, the standard error of the log2FC estimate. It is an
interval of sorts: how tightly the number is pinned down, the way you would judge agreement between
replicates.

`lfcSE` is the uncertainty of the `meanL2FC` estimate itself, that is how sure DESeq2 is about the
number it produced. The fewer reads a barcode has, the wider that uncertainty.

About the quality-versus-threshold plot below. The models are weak, but some threshold had to be
picked so that they would learn at least something, and it had to be roughly uniform so that the
models stay comparable to each other. So 1.5 was chosen for everything, and a looser 4.2 for two
constructs, PACRGL Control and TONSL, because at stricter thresholds there is simply no data left
for them. Above 4.2 the data looks like junk to me anyway.

**1.5 for the other three** (NBEAL2, DDX51 and the merged `ALL`) - that is where classification
quality peaks: NBEAL2 0.822, DDX51 0.721.

**4.2 for the two others** is the only value at which they have any training examples at all, and it
was chosen on the availability of data, not on the result.

For PACRGL Control and TONSL the 4.2 threshold removes no rows at all. The data for these two
constructs is just noisier.

The volcano plot in the final section is the one place that does not follow this. It applies the same
`lfcSE` cutoff of 1.5 to all four constructs, because panels built at different cutoffs cannot be
compared: the smallest detectable effect scales with the cutoff. At the common cutoff the counts are
226 and 204 selected drops for DDX51 and NBEAL2 against 18 and 23 for PACRGL Control and TONSL, so
the number shown on those two panels is smaller than the class 1 count their models were trained on.

![AUC against the lfcSE filter](other_pictures/son_auc_vs_lfcse.png)

### Labelling rules

A word about the data itself. What do you do about a classification on `meanL2FC < 0`? There are no
exact `log2FC = 0` values in the table, so a zero has to be built artificially. For that I took the
barcodes whose p-value is not significant and whose log2FC is within 0.1 of zero, and called those
zero.

At first I tried a three-class problem - below zero, zero, above zero - and the model did not work
at all. Regression on both positive and negative values did not work either. I talked it over with
Matthias, and since the data is so noisy we decided to predict one thing only: whether splicing goes
down under SON depletion.

**Classification, dropped or not:**

| | rule |
|---|---|
| class 1 | `meanL2FC < 0` **and** `minPvalue < 0.005` |
| class 0 | `minPvalue >= 0.005` **and** \|`meanL2FC`\| `< 0.1` |
| not used for training | everything else |

## Class definition

![class definition](other_pictures/son_class_definition.png)

| construct | `lfcSE` | rows | class 1 | share |
|---|---|---|---|---|
| PACRGL Control | 4.2 | 1 421 | 295 | 21 % |
| DDX51 | 1.5 | 781 | 226 | 29 % |
| NBEAL2 | 1.5 | 897 | 204 | 23 % |
| TONSL | 4.2 | 980 | 195 | 20 % |
| `ALL`, merged | 1.5 | 2 655 | 471 | 18 % |

**Regression, how far it dropped:** rows with `meanL2FC < 0` **and** `minPvalue < 0.005`, so the
same barcodes that classification calls class 1. The target is the measured `meanL2FC` as it is,
negative. No absolute value is taken and there are no artificial zeros.

| construct | `lfcSE` | rows | target |
|---|---|---|---|
| PACRGL Control | 4.2 | 295 | -24.25 to -0.01 |
| DDX51 | 1.5 | 226 | -6.75 to -0.25 |
| NBEAL2 | 1.5 | 204 | -8.14 to -0.29 |
| TONSL | 4.2 | 195 | -26.39 to -0.32 |
| `ALL`, merged | 1.5 | 471 | -8.45 to -0.25 |

An earlier version of this table gave 10 323 / 3 139 / 3 436 / 7 608 / 11 360 rows. Those came from
a build that kept every negative barcode with no significance filter. It was replaced, and the
numbers above are what every current figure is made from.

### What you need to know when reading these numbers

**`ALL` was built at 1.5, so it is almost entirely NBEAL2 plus DDX51.** PACRGL Control and TONSL
contribute nine positives each.

**Classification quality goes hand in hand with how much the prediction tracks coverage. Coverage,
though, is not a pure artefact.** The size of the effect is tied to coverage very strongly, but the
direction is not:

- Spearman of \|`meanL2FC`\| with coverage (`meanBaseMean`), with no selection at all: **-0.479**
- the same thing inside the significant set: **-0.822**
- Spearman of `lfcSE` with coverage: **-0.835**
- but the link between the **direction** of the change (dropped or not) and coverage: **+0.010**,
  which is essentially nothing


---

# Data

## meanL2FC distributions

![passed the threshold](other_pictures/son_l2fc_distributions.png)

![did not pass the threshold](other_pictures/son_l2fc_nonsignificant.png)

## meanL2FC by construct

Three views of the same data, one panel per construct, getting narrower from top to bottom.

First, everything that was measured at all, with no filter of any kind. Note the y axis here goes to
2 500, against 100 on the two below.

![every measured barcode](other_pictures/son_l2fc_per_construct_all.png)

| construct | measured | of those, down | share down |
|---|---|---|---|
| PACRGL Control | 20 326 | 10 323 | 50.8 % |
| DDX51 | 12 430 | 6 388 | 51.4 % |
| NBEAL2 | 12 111 | 6 797 | **56.1 %** |
| TONSL | 15 018 | 7 608 | 50.7 % |

Before any selection all four constructs sit close to 50/50, NBEAL2 being the only one a little off
it. Whatever direction shows up later comes from the significance threshold, not from the raw data.

Then the barcodes that pass the significance threshold.

![significant only](other_pictures/son_l2fc_per_construct_loose.png)

![significant and precisely measured](other_pictures/son_l2fc_per_construct_loose_lfcse.png)

| construct | `p < 0.005` | plus its own `lfcSE` | rows dropped by the `lfcSE` filter |
|---|---|---|---|
| **PACRGL Control** | 569 | 569 | **0** |
| DDX51 | 580 | 356 | **224** |
| **NBEAL2** | 428 | 279 | **149** |
| TONSL | 325 | 325 | **0** |

The share of drops among the significant barcodes runs from about 52 % in PACRGL Control to about
75 % in NBEAL2, and the `lfcSE` filter moves it by no more than a couple of points either way.

`PACRGL Control` is the control construct, and it comes out at exactly **50/50**: there is no
directional effect where there should not be one.

The shapes of the distributions differ. `DDX51` and `NBEAL2` have a sharp peak at zero, while
`PACRGL Control` and `TONSL` have a dip there. Most likely this is coverage: at high coverage a
small effect passes the significance threshold, at low coverage it does not.

## Training set sizes

Exactly the rows the model saw, after the intragenic/intergenic split and the training script's own
internal cutoffs.

![classification](other_pictures/son_dataset_sizes_classification.png)

![regression](other_pictures/son_dataset_sizes_regression.png)

Classification, inside MANE genes:

| model | `lfcSE` | rows | class 1 | share | train | test | class 1 in test |
|---|---|---|---|---|---|---|---|
| PACRGL Control | 4.2 | 653 | 110 | 17 % | 551 | 102 | **20** |
| DDX51 | 1.5 | 354 | 119 | 34 % | 299 | 55 | **17** |
| NBEAL2 | 1.5 | 420 | 90 | 21 % | 354 | 66 | **13** |
| TONSL | 4.2 | 460 | 71 | 15 % | 388 | 72 | **10** |
| `ALL`, merged | 1.5 | 1 542 | 259 | 17 % | 1 300 | 242 | **42** |

Classification, outside MANE genes:

| model | `lfcSE` | rows | class 1 | share | train | test | class 1 in test |
|---|---|---|---|---|---|---|---|
| PACRGL Control | 4.2 | 699 | 172 | 25 % | 589 | 110 | **24** |
| DDX51 | 1.5 | 319 | 65 | 20 % | 269 | 50 | **8** |
| NBEAL2 | 1.5 | 299 | 84 | 28 % | 252 | 47 | **20** |
| TONSL | 4.2 | 389 | 111 | 29 % | 329 | 60 | **17** |
| `ALL`, merged | 1.5 | 1 043 | 195 | 19 % | 880 | 163 | **28** |

Regression, rows and how they split:

| model | `lfcSE` | inside, rows / train / test | outside, rows / train / test |
|---|---|---|---|
| PACRGL Control | 4.2 | 88 / 75 / **13** | 133 / 112 / **21** |
| DDX51 | 1.5 | 67 / 57 / **10** | 36 / 31 / **5** |
| NBEAL2 | 1.5 | 60 / 51 / **9** | 42 / 36 / **6** |
| TONSL | 4.2 | 37 / 31 / **6** | 28 / 24 / **4** |
| `ALL`, merged | 1.5 | 201 / 170 / **31** | 154 / 130 / **24** |

The bold column is the one that limits what can be measured at all. ROC AUC is computed on the
test barcodes, and with 8 to 42 positives there its confidence interval is wide no matter how
good the model is.

---

# Quality

| | |
|:---:|:---:|
| ![quality, regression](other_pictures/son_quality_regression.png) | ![quality, classification](other_pictures/son_quality_classification.png) |

| construct | `lfcSE` | table | classification train | classification test | classification Δ | regression train | regression test | regression Δ |
|---|---|---|---|---|---|---|---|---|
| PACRGL Control | 4.2 | inside genes | 0.601 | 0.592 | 0.008 | 0.578 | -0.177 | 0.755 |
| PACRGL Control | 4.2 | outside genes | 0.654 | 0.604 | 0.050 | 0.774 | -0.021 | 0.795 |
| DDX51 | 1.5 | inside genes | 0.841 | 0.721 | 0.120 | 0.783 | 0.382 | 0.401 |
| DDX51 | 1.5 | outside genes | 0.703 | 0.507 | 0.196 | 0.471 | 0.224 | 0.247 |
| NBEAL2 | 1.5 | inside genes | 0.837 | 0.822 | 0.015 | 0.758 | 0.367 | 0.391 |
| NBEAL2 | 1.5 | outside genes | 0.589 | 0.504 | 0.085 | 0.910 | 0.486 | 0.424 |
| TONSL | 4.2 | inside genes | 0.672 | 0.610 | 0.062 | 0.960 | -0.143 | 1.103 |
| TONSL | 4.2 | outside genes | 0.810 | 0.627 | 0.183 | 0.795 | 0.949 | -0.154 |
| ALL | 1.5 | inside genes | 0.699 | 0.626 | 0.073 | 0.811 | 0.456 | 0.354 |
| ALL | 1.5 | outside genes | 0.666 | 0.588 | 0.078 | 0.595 | 0.136 | 0.459 |

## What the classifier actually predicts

An AUC is one number and it hides what the model does. Here every test barcode is a dot, placed
at the probability the model gave it. The upper row in each lane is true class 0, the lower row is
true class 1. The dashed line is the threshold that separates the two classes best, and the columns
on the right say how much of each class it gets right there.

![predicted probabilities](other_pictures/son_predicted_probabilities.png)

Three things show up here that the AUC does not say.

Six models out of ten produce almost no spread. Four models do work. The merged model ranks correctly but with almost no margin.


## Checking against chance

The bar is the 95 % bootstrap interval, that is the precision of the measurement. The grey band is
the 95 % spread under shuffled labels, that is what a dud would look like.

![is the AUC real](other_pictures/son_auc_significance.png)

The same panel for inside MANE genes only, with the number of positives in the
test under each construct.

![AUC against chance, inside MANE genes](other_pictures/son_auc_significance_intragenic.png)

For regression, the link between the prediction and sequencing depth is drawn next to the Spearman
value. Without it, a high Spearman cannot be told apart from simply predicting coverage.

So does coverage sometimes beat the features as a predictor? Partly yes, and it is worth being
precise about what that does and does not mean.

Coverage here is `meanBaseMean`, the mean of the `baseMean` of the two DESeq2 pools, that is how
many reads a barcode gave. It runs from 3 to 54 516, with a median of 69.

**The model does not see it.** The training table has 77 columns and `meanBaseMean` is not one of
them. It sits apart, in `DATA/son_splicing/targets/*.targets.tsv`.


So we trained a control model that predicts coverage from the same features on the same rows.

One caveat about that control. `NMD_expected` was dropped from the SON features on 2026-09-11, so
the real models now run on 14 feature groups inside genes while the control still runs on 15. The
control was deliberately left as it was. A Shapley value splits the payout between the players, so
strictly speaking the two profiles are no longer the same game, and the comparison below is a little
looser than it was. It does not change the direction of the answer: the match was at the null level
to three decimals, and one group out of fifteen cannot move that.

![is the Spearman real](other_pictures/son_spearman_significance.png)

---

# SHAP

## By feature group

| | inside MANE genes | outside MANE genes |
|---|---|---|
| regression | ![](other_pictures/group_shap_heatmap_son_reg_regression_intragenic.png) | ![](other_pictures/group_shap_heatmap_son_reg_regression_intergenic.png) |
| classification | ![](other_pictures/group_shap_heatmap_son_cls_classification_intragenic.png) | ![](other_pictures/group_shap_heatmap_son_cls_classification_intergenic.png) |

### The same, normalised within each column

Four panels where every column is divided by its own sum. Columns add up to 1, and the colour means
the **share of the model's attention**.

| | inside MANE genes | outside MANE genes |
|---|---|---|
| regression | ![](other_pictures/group_shap_heatmap_share_son_reg_regression_intragenic.png) | ![](other_pictures/group_shap_heatmap_share_son_reg_regression_intergenic.png) |
| classification | ![](other_pictures/group_shap_heatmap_share_son_cls_classification_intragenic.png) | ![](other_pictures/group_shap_heatmap_share_son_cls_classification_intergenic.png) |

**The raw version says where the model found any signal at all**, and a pale column there is honest.
**The share version says how the attention is made up** where there is signal.

The two have to be read together. In a column where only one or two groups are non-zero, a large
share means very little on its own, so it is worth checking the raw panel before reading anything
off the share one.

## By individual feature

![regression, inside genes](other_pictures/son_featureshap_regression_intragenic.png)

![regression, outside genes](other_pictures/son_featureshap_regression_intergenic.png)

![classification, inside genes](other_pictures/son_featureshap_classification_intragenic.png)

![classification, outside genes](other_pictures/son_featureshap_classification_intergenic.png)

### The same, normalised within each column

Same idea as for the groups: every column is divided by its own sum, so the columns add up to 1 and
the colour is the share of the model's attention.

![regression, inside genes](other_pictures/son_featureshap_share_regression_intragenic.png)

![regression, outside genes](other_pictures/son_featureshap_share_regression_intergenic.png)

![classification, inside genes](other_pictures/son_featureshap_share_classification_intragenic.png)

![classification, outside genes](other_pictures/son_featureshap_share_classification_intergenic.png)

The raw version shows where the model found signal in the first place, the share version shows how
the attention is composed where signal exists.
