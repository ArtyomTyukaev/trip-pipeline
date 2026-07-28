# Figures

The seven figures referenced from [../LAB_NOTEBOOK.md](../LAB_NOTEBOOK.md). Both sets were produced
by an `aggregate_and_plot.R` script that parses the per-job logs into a `summary.csv` and plots it;
the underlying `summary.csv` tables are not in this repository yet.

## Stage 04 — TTAA filter, 10 real Map libraries (run3 + run5)

| File | What it shows |
|---|---|
| `ttaa_fig1_pct_dropped.png` | per library, % of barcodes and % of reads dropped by the standard proximity filter vs by the TTAA filter |
| `ttaa_fig2_barcodes_vs_reads.png` | the same libraries as barcodes-dropped against reads-dropped: the gap between the two is the whole point (~29 % of barcodes, ~7.6 % of reads) |
| `ttaa_fig3_totals.png` | pool totals over 568,687 positions / ~27.9 M reads |

## Stage 07 — `oneMismatch` flag, 16 set1 combos

| File | What it shows |
|---|---|
| `mm_fig1_pct_gain.png` | per job, % gain in barcodes, norm reads and exp reads |
| `mm_fig2_totals.png` | pool totals, without the flag (exact) vs with it (allow 1 mismatch) |
| `mm_fig3_scatter_existing.png` | per-barcode exp reads for barcodes present in **both** outputs (`Exp_5ss_noP_rep1_S99`): exact on X, allow-1-mismatch on Y; above the diagonal = reads gained, the band at `x=0` = pre-existing barcodes that had `exp=0` in exact mode |
| `mm_fig4_newbc_exp_hist.png` | exp distribution of the 728 barcodes that appear **only** with the flag (same sample); the bar at `exp=0` (~410) are barcodes pulled in by rescued norm reads alone |

A note on how these were rendered: direct PNG devices in R (`ragg`, `cairo`) crash with SIGILL /
exit 132 on some cluster nodes, so the plots are drawn to PDF and converted with Ghostscript
(`gs -sDEVICE=png16m`).
