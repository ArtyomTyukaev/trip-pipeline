# Preprocessing steps

Which file joins to which, by which key, in what row order, for all four datasets.

Project root `/home/projects/ulitsky/artemti/2026/TRIP_SCRIPTS/`. Paths are relative to it
unless they start with `/home/`.

---

## 1. The two notebooks

```
DATA/<dataset>/prepare_dataset.ipynb    reads the read-only source tree
                                        -> DATA/<dataset>/   (6 sub-folders)

1) Data Preprocessing.ipynb             reads DATA/<dataset>/
                                        -> tables_for_training/<dataset>/{exp}.alldata.tab
```

Notebook 1 serves all four datasets; the path block is chosen by uncommenting one cell, and
`feature_types` must match — `['histone','speckle']` for the three new datasets,
`['histone','speckle','u2af2','stranded']` for `set1_combined_1`.

What it does per experiment:

```
annotated_tables/{exp}.annotated_tables.txt    index_col=0        -> index is the barcode
combined/*{exp}*combined*{ftype}*.txt          header=None        -> bare numbers
    column names <- columns_names/{ftype}_names.txt       BY POSITION, no check
    row index    <- the barcode order file of that ftype  BY POSITION, length checked only
barcodes_bed/{exp}*.bed                        index_col=3        -> chr,start,end,score,strand

common = set.intersection(annotated index, each ftype barcode list)
final  = concat(annotated, each ftype), restricted to `common`, each sorted by barcode
final  = final[selected_columns]        concat of every file in selected_columns/
final  = merge(bed, final, left_index=True, right_index=True, how='inner')
```

**So there are two different join mechanisms.** Signal values attach to barcodes purely by row
position — the matrices contain no barcodes at all, and the only guard is that the lengths
match. Everything after that joins by barcode, and a barcode survives only if it is in the
annotated table *and* every feature family *and* the hg38 bed.

`get_barcode_file(exp, ftype)` fixes the family → order-file mapping:

| ftype | barcode order file |
|---|---|
| histone | `{exp}.hg38_barcodes.txt` |
| speckle, u2af2 | `{exp}.hg19_barcodes.txt` |
| stranded | `{exp}.stranded_barcodes.txt` |

Output column order = the 5 bed fields, then the `selected_columns/` files in `os.listdir`
order (`histone, u2af2, stranded, speckle, annot_tab` for set1; `histone, speckle, annot_tab`
for the rest).

---

## 2. Datasets at a glance

| dataset | experiments | source tree |
|---|---|---|
| `single_combined` | 2 (`5ss_noP`, `CPA_noP`) | `igoru/data/human/TRIP/single_combined/` |
| `intron_exp` | 48 = 8 introns × 3 sample types × 2 pools | `shared/data/human/TRIP/IntronTRIP_combined/` |
| `intron_splicing` | 16 = 8 introns × 2 pools | same tree |
| `set1_combined_1` | 40 of 48, see §6 | `lidiag/TRIP_data/set1_combined_1/` |

`igoru/data/human/TRIP/` and `shared/data/human/TRIP/` are the same directory through two mounts.

| dataset | output rows | cols | make-up |
|---|---|---|---|
| `single_combined` | 40,286 and 38,520 | 80 | 5 coords + 30 histone + 21 speckle + 22 features + `NMD_expected` + target |
| `intron_exp` | 25,783 … 159,802 | 79 | 5 coords + 30 histone + 21 speckle + 22 features + target |
| `intron_splicing` | 9,405 … 39,203 | 79 | same as `intron_exp` |
| `set1_combined_1` | 811 … 31,930 | 127 | 5 coords + 30 histone + 41 u2af2 + 15 stranded + 21 speckle + 14 features + target |

The barcode is an unnamed index column, so a header line has one field fewer than a data line.

| dataset | `exp.log10` holds | `LOG_TARGET` | classification |
|---|---|---|---|
| `single_combined` | raw `meanRPM`, 0 … 40,607 | **1** | binary, `BINARY_THRESHOLD = 0.1` |
| `intron_exp` | the sample's own column, already `log10(1 + RPM)` | off | 3 classes, median of non-zeros |
| `intron_splicing` | `exp / (exp + exp_unspliced)`, 0 … 1 | off | 3 classes at 0.6 and 0.8 |
| `set1_combined_1` | log expression | off | 3 classes, median of non-zeros |

The target is named `exp.log10` everywhere because both training scripts hardcode that name. In
`single_combined` and `intron_splicing` the name does not describe the contents.

`CATEGORICAL_COLUMNS` for the three new datasets:
`ChromHMM,SPIN,GC_category,RepeatGroup,annot_rna_type,Chr`. For set1 it is the script default
`ChromHMM`, the only factor among its 14 features. The 0.6/0.8 rule is applied by name from the
`introns_splicing` list in `MODELS/XGB_classifier/XGB_DS_Classifier.py`, currently all 16 pool
names; the script raises if `BINARY_THRESHOLD` is set while one of them is in the batch.

Coordinates, from `barcodes_bed/`:

| dataset | bed | width | integration point |
|---|---|---|---|
| `single_combined`, `intron_exp`, `intron_splicing` | `*.plus_50_each_side.bed` | 101 | `start + 50` |
| `set1_combined_1` | `{exp}.bed` | 2 | `start` |

`make_folds.py` builds overlap components from `[start - 50, end + 50)`, so the effective
window is 201 nt for the first three and 102 nt for set1.

---

## 3. `single_combined`

```
{EXP}.txt.processed_table.txt        40,386 x 92 and 38,614 x 96, index = barcode
    -> annotated_tables/{EXP}.annotated_tables.txt   24 cols
       10 factor/is_na + 12 numeric + NMD_expected + exp.log10 (= raw meanRPM)

NMD/output/nmd_table_{EXP}.tsv       joined on barcode, covers every barcode
    NMD_expected = (nmd_expected == 1); positives 7,927/40,386 and 7,809/38,614

meta/{EXP}.{track}.column_histone.txt   30 files, all kept  -> combined/, 30 cols
meta/{EXP}.{track}.column_speckle.txt   29 files, 21 kept   -> combined/, 21 cols
    each file = field 4 (mean0) of meta/{EXP}.{track}.tab.gz, same order, no barcodes

barcodes_order/{EXP}.hg38_barcodes.txt  = field 0 of any histone .tab.gz   58,172 / 57,303 rows
barcodes_order/{EXP}.hg19_barcodes.txt  = field 0 of any speckle .tab.gz   57,850 / 56,980 rows

{EXP}.plus_50_each_side.bed        -> barcodes_bed/{EXP}.bed   (hg38, 101 nt)
hg19.{EXP}.plus_50_each_side.bed   not copied, only its barcode set is checked
```

The 100 and 94 barcodes lost against `processed_table` are exactly those absent from the hg19
bed. Target spread, which is what drives `LOG_TARGET` and the 0.1 threshold: `5ss_noP` median
0.479, max 40,608, 14,988 exact zeros, 58.9 % above 0.1; `CPA_noP` median 1.062, max 24,941,
12,685 zeros, 62.0 % above 0.1.

---

## 4. `intron_exp`

Introns A1, B1, B3-6, B3-9, C1, D1, E1, H1.

```
{intron}.txt.processed_table.txt      index = barcode, 60 named columns
    sample columns = header fields after `closestCAGEdist`,
        minus anything starting with `mean`,
        minus anything ending in .RPM / .normRPM / .mapRPM
        -> exactly 6 per intron (Spliced, Total, Unspliced x 2 pools)
    experiment name = that column name, with `Toal` -> `Total` (H1 only)
    -> DATA/intron_exp/sample_map.tsv   experiment / intron / source_column
    -> annotated_tables/{exp}.annotated_tables.txt   23 cols
       exp.log10 = that sample column; rows where it is NaN are dropped

meta/{intron}.{track}.column_{histone,speckle}.txt   -> combined/{intron}.combined_{kind}.txt
barcodes_order/{intron}.{hg38,hg19}_barcodes.txt     = field 0 of the matching .tab.gz
{intron}.plus_50_each_side.bed                       -> barcodes_bed/{intron}.bed
```

The six sample columns are found by position because their names follow three different
conventions:

```
IntronTRIP_Spliced_A1_Pool_1_1_S25      A1, H1
Spliced_B1_1_1_PolyA_S19                B1, C1, D1
TRIP_RNA_spliced_B3-6_1-1_S53           B3-6, B3-9, E1
```

**Chromatin belongs to the intron, not the experiment.** The six experiments of one intron
share it, so the real file sits under the intron name and every experiment name is a symlink
beside it: `combined/` 112 files of which 96 symlinks (16 real = 8 introns × 2 families),
`barcodes_order/` the same, `barcodes_bed/` 56 of which 48 symlinks. Only `annotated_tables/`
is 48 real files. Copying instead would duplicate up to 95 MB six times per intron.

A missing sample value means the integration is not in that pool of cells — zero means it made
no RNA. The two pools of one intron share 2.9 % (B3-9) to 8.5 % (H1) of their barcodes, so
roughly half of each sample column is NaN. Those rows are dropped at build time, because
`get_gene_tables.py` runs `fillna(0)` on the intragenic table and would turn them into
fabricated zero measurements.

No `NMD_expected` here — that prediction is built for the 5ss/CPA construct. `meanRPM` is
deliberately kept out of `selected_columns/`; it is a summary of the Total samples.

---

## 5. `intron_splicing`

Same 8 introns and same 22 features as `intron_exp`; one experiment per intron and pool.

```
spliced_unspliced.txt      the driver, 16 lines, no header
    field 0 = experiment name, fields 1..3 = spliced / unspliced / total sample names
    splicing.txt beside it is STALE: 10 lines, the five older introns only

{intron}_{pool}.spliced_unspliced.combined_counts.txt      index = barcode
    columns: norm, exp, Chr, Pos, Strand, map, exp_unspliced, exp_total, exp.log10
    exp.log10 = exp / (exp + exp_unspliced), NA where exp + exp_unspliced < 20
    -> the target, after dropna(), joined to the intron's feature frame on barcode

{intron}.txt.processed_table.txt   -> the same 22 features, no readout taken from it
combined/, barcodes_order/, barcodes_bed/   as in intron_exp, built independently
```

**Why the target comes from a different file:** a fraction needs raw counts, and the sample
columns in `processed_table` are logs of normalised values. `exp.log10` in the counts file is
already the fraction, not a log.

Two sibling ratios sit next to it. Against last year's training target on `A1_1`:
`spliced_unspliced` r = 0.9998, `spliced_to_total` r = 0.76,
`spliced_and_unspliced_to_total` r = −0.18. We use `spliced_unspliced`.

The join with the feature frame is an intersection rather than an assertion because one barcode
of `E1_1` is missing from `E1`'s `processed_table` — see §8, first entry.

Checked at build time: target inside 0 … 1, equal to `exp / (exp + exp_unspliced)` recomputed
from the raw counts to within 1e-12, and `exp + exp_unspliced >= 20` on every kept row.

---

## 6. `set1_combined_1`

`prepare_dataset.ipynb` here rebuilds **only** `combined/` and `columns_names/`. The other four
folders and the two name lists are byte-identical copies of the `lidiag` files, `md5sum`
checked.

```
{exp}.annotated_tables.txt    copied as is, 48 cols, exp.log10 is a real log expression
    selected_columns/annot_tab_selected_columns.txt keeps 15:
        norm, MANE_all_exp.norm.log10, MANElength, MANE_intron_length, DistToTSS, DistToTES,
        RelPos, DistToIntron5, DistToIntron3, RelIntronPos, ChromHMM, GC, MANE_all,
        MANE_exonic, exp.log10

meta/{exp}/{exp}.{track}.column_{kind}.txt      per-track files live in a PER-EXPERIMENT
                                                sub-directory, unlike the intron tree
    -> combined/{exp}.combined_{kind}.txt       192 matrices = 48 x 4 families, ~1 GB,
                                                no symlinks, every experiment its own
{exp}.bed                     -> barcodes_bed/{exp}.bed    2 nt intervals
{exp}.{hg38,hg19,stranded}_barcodes.txt  -> barcodes_order/, copies, not regenerated
```

| family | per-track files | kept | order file |
|---|---|---|---|
| histone | 30 | 30 | `hg38_barcodes` |
| speckle | 29 | 21 | `hg19_barcodes` |
| stranded | 15 | 15 | `stranded_barcodes` |
| u2af2 | 89 | 41 | `hg19_barcodes` |

The 8 speckle tracks dropped are 7 `GSE81553_NoPrimaryCtrl_for_*` backgrounds plus
`GSE81553_SON_TSA-Seq_Condition1_Replicate_2nd_pulldown`.

`{exp}.stranded_barcodes.txt` is the plus-strand block followed by the minus-strand block. On
`3intron_5ss_noP_Pool1_S95` it holds 15,158 + 14,928 = 30,086 barcodes — the same set and the
same length as `hg38_barcodes.txt`, agreeing with it at 2 of 30,086 positions. Swapping the two
attaches every stranded value to the wrong barcode and nothing fails. Row counts on that
experiment: histone and stranded 30,086, speckle and u2af2 29,996, output 29,996. Every
experiment's output row count is its hg19 barcode count.

### 8 of the 48 experiments are excluded

`names_all_clean.txt` lists **40**; the full 48 are kept in `names_all_clean_all48.txt`.
Excluded: `3intron_CPA_noP_Pool2_S94` and `..._Pool2_Rep2_S105` (811 rows each),
`3intron_CPA_noP_Pool1_S93` and `..._Pool1_Rep2_S104` (900), `5ss_Expression_rep1_S56` (1,105),
`5ss_Expression_rep2_S57` (1,151), `CPA_Expression_rep2_S54` (2,017),
`CPA_Expression_rep3_S55` (2,046).

**Too few barcodes.** The next smallest of the remaining 40 has 29,545 rows, so this is not a
judgement call. At 811 rows a 15 % holdout plus 3 folds leaves ~230 rows per fold, which makes
metrics and SHAP noise.

**A mapping artifact, concentrated in those four `*_Expression_*` tables.** In every set1
experiment 1.4–1.8 % of barcodes sit inside one 460 bp window at
`chr2:32,916,183-32,916,645` — 241 distinct positions for 477 rows in `Exp_5ss_noP_rep1_S99`.
There `DNAse` is ~52 against a genome-wide mean of 0.13, `ATAC` ~0.09 against 1.06, and
`ChromHMM = EnhBiv` for all of them. High DNase with no ATAC in a 460 bp spot is a blacklist or
repeat mapping artifact: the annotation lookup is correct, the coordinates are wrong, so the
rows are internally consistent and no formal check flags them. The same barcodes map to
scattered real positions in `single_combined`.

Those rows hold **86 % of the total DNase signal mass** and drag `rho(DNAse, ATAC)` from +0.19
to +0.07. In the four `*_Expression_*` tables they are 5–7 % of rows and the correlation flips
sign, to −0.15 … −0.22. The window is **not** filtered out of the remaining 40, where it is
~1.5 % of rows, so any conclusion about `DNAse` or `H3K9me1` from `set1_combined_1` is
unreliable. Every other feature is unaffected, and the pile-up does not exist in the other
three datasets.

The 8 tables stay on disk; nothing reads them, because the experiment list drives notebook 1
and the training scripts.

---

## 7. Rules shared by the three new `annotated_tables/` builders

`single_combined`, `intron_exp` and `intron_splicing` write their own; set1 copies Igor's.

**R row names.** `processed_table.txt` and `combined_counts.txt` were written by R with row
names, so the header has one field fewer than the data lines. Pandas promotes field 0 to the
index — the barcode — so `pd.read_csv(..., sep='\t')` works. `awk`, `header=None` and R with
`header=TRUE` all shift every column by one.

**Four `is_na` features.** The source column holds a gene name or nothing; the feature is the
0/1 presence flag, and the name is discarded because as a category it would be thousands of
levels. `MANE_antisense` and `gencode_noncoding` are binarised in `prepare_dataset.ipynb`;
`MANE_all` and `MANE_exonic` stay raw strings because the training scripts binarise them
themselves. Both are in `settings/columns_to_drop_in_intragenic.txt`, so they only reach a
model on the intergenic table.

**Missing categories are filled with the literal string `absent`** in `annot_rna_type`, `SPIN`,
`GC_category`, `RepeatGroup`. Not `NA`, which `pd.read_csv` reads back as missing and which
silently gives the intragenic and intergenic tables different one-hot columns for the same
feature. Left empty, one-hot encoding emits a row of zeros and "not in a repeat" becomes
indistinguishable from "no data".

**`ChromHMM` keeps its NaN,** and numeric NaN is left alone; `get_gene_tables.py` drops NaN
`ChromHMM` rows, fills the rest with 0 on the intragenic table, and drops ChromHMM states with
fewer than 10 barcodes. **The target must never be NaN,** because of that same `fillna(0)`.

**`selected_columns/` is a hard filter.** Notebook 1 concatenates every file in the folder and
runs `final_df[selected_columns]`, which raises `KeyError` on a name the table lacks. So only
the families the dataset actually has may be listed, and the target must be listed or it is
dropped without a message. Annotation list lengths: `single_combined` 24 (23 features),
`intron_exp` / `intron_splicing` 23 (22), `set1_combined_1` 15 (14).

---

## 8. Traps

**Every library loses the first barcode of its bed at the annotation stage.**
`igoru/repos/trip-unified-pipeline/scripts/11_table_annotation_genuine_barcodes.R:91` reads the
headerless BED6 with `header = T`, so the first integration is consumed as a column-name row.
The next line overwrites the names with the correct hardcoded ones, so there is no error and no
warning — just one row fewer, and the table starts at bed line 2. Confirmed for all 8 introns
and both `single_combined` libraries. It becomes visible in the splicing join for `B1_1`,
`E1_1` and `E1_2`, where the lost barcode also reached the counts file. Harmless in size (1
barcode of 25k–380k) and it causes no misalignment, since everything downstream joins by
barcode. Needs a one-character fix on Igor's side.

**`meta/{prefix}.combined_{kind}.txt` in every source tree — NOT USED.** Glued with
`paste meta/$p.*column_{kind}.txt`, so the columns are in shell-glob alphabetical order, while
the curated name lists (`histone_marks.txt`, `speckle_names.txt`, `stranded_names.txt`,
`u2af2_names.txt`) are not. Notebook 1 applies names by position and checks nothing. Labels
land on the wrong tracks — measured by matching every matrix column back to the per-track file
it equals byte for byte:

| tree | histone | speckle | stranded | u2af2 |
|---|---|---|---|---|
| `lidiag`, `3intron_5ss_noP_Pool1_S95` | 24 of 30 wrong | 0 of 29 | 13 of 15 wrong | 0 of 89 |
| `igoru`, `5ss_noP` | 27 of 30 wrong | 14 of 29 wrong | — | — |
| `IntronTRIP`, `A1` | 27 of 30 wrong | 14 of 29 wrong | — | — |

The two trees were globbed in different collations, so no single Python sort reproduces both.
`settings/custom_functions/build_combined_from_columns.py` takes each column name from its own
file name, and `write_combined()` writes the matrix and the name list from one column list, so
they cannot drift.

**The bed is not the row order of the per-track files.** The bed is sorted `chr1, chr10,
chr11, …` and `bigWigAverageOverBed` writes `chr1, chr2, chr3, …`. Same barcodes, same count,
so every length and set check passes while the values are permuted. Positional agreement
between bed field 3 and `meta/{prefix}.{track}.tab.gz` field 0 is only ~14 % for hg38 and
0.3–3.2 % for hg19, across all ten prefixes. The order therefore comes from the `.tab.gz`.
set1 keeps Igor's `*_barcodes.txt`, which he regenerated in September 2025 after finding this.

**`hg19` is a liftOver of `hg38`** and is shorter; its interval widths run 97–101 instead of a
flat 101. It is used only for the speckle and u2af2 row order and is never copied into
`barcodes_bed/` — notebook 1 looks the bed up as `{exp}*.bed`, which never matches `hg19.*`.

**The `3intron` library type is missing from the intron data.**
`IntronTRIP_combined/combos.txt` — the driver deciding which samples enter `processed_table` —
is byte-for-byte `run2/combos.txt` + `run4/combos.txt`, so run3 is not wired in at all. That
drops its 10 `3intron` samples plus its own 10 `UNS` and 5 `Total_totalRNA`. And `3intron` is
not a renamed `Spliced`: on shared barcodes the two splicing fractions barely correlate (`A1_1`
r = −0.09 on 5,254; `H1_1` r = 0.23 on 3,041). Six of the ten `3intron.combined_counts.txt`
were rebuilt in the combined folder on 2026-07-26; the four for A1 and E1 were not, because the
Norm libraries the run3 driver names for them do not exist there. Needs Igor.

**A fresh timestamp is not a fix.** Igor rebuilt the glued matrices on 2026-08-12 for an
unrelated reason and the column order did not change. He updated the `processed_table` files
around the same date because of a row-filter change.

**Old `*.annotated_tables.txt` from June–July 2025** are still on disk in the IntronTRIP tree
next to the current files, with names that look like plausible experiment names.

---

## 9. Checks the builders run

All are `assert`, so a mismatch stops the build.

```
every column of every written matrix == its own per-track file
                 exact equality with a NaN-equal clause, not a tolerance
matrix columns   == length of columns_names/{kind}_names.txt
matrix rows      == length of the barcode order file of that family, looked up explicitly
barcode order    unique, agrees across the first 3 tracks of the family, equals the bed as a
                 set, and round-trips through the written file positionally
target           equals its source column, holds no NaN
                 splicing only: inside 0..1, equals exp / (exp + exp_unspliced),
                 exp + exp_unspliced >= 20
meanRPM          not among the features
symlinks         every experiment name resolves
```

Two checks on the numbers themselves, which catch a label shift that survives all of the above:

- `DNAse` is a peak track and must be under 20 % non-zero. Measured 7.1 % in
  `set1_combined_1`, 5.4 % in `single_combined`, 6.2 % in `intron_exp`; `H2AFZ` next to it is
  100 % non-zero in all three.
- Inside the **histone** family, `phyloP` must be the only column with negative values.
  Confirmed in all three. Not applied to speckle, where the TSA-seq and HiC eigenvector tracks
  are signed by design.

---

## 10. Current state

`DATA/` is built for all four datasets. `tables_for_training/` holds `single_combined` (2
tables), `intron_exp` (48), `intron_splicing` (16) and `set1_combined_1` (48).
