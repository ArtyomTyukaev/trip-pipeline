# Figures

The six figures referenced from [../LAB_NOTEBOOK.md](../LAB_NOTEBOOK.md). Both sets were produced
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
| `mm_fig2_totals.png` | pool totals, exact vs `oneMismatch` |
| `mm_fig3_scatter_exp.png` | per-barcode exp reads, exact vs `oneMismatch` (`Exp_5ss_noP_rep1_S99`); the points along `x=0` are barcodes that exist only once 1-nt neighbours are rescued |

A note on how these were rendered: direct PNG devices in R (`ragg`, `cairo`) crash with SIGILL /
exit 132 on some cluster nodes, so the plots are drawn to PDF and converted with Ghostscript
(`gs -sDEVICE=png16m`).
