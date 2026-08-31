"""Expression against splicing — on the TRUE group Shapley value.

    python analysis/group_shap/code/plots/plot_group_shap_exp_vs_splicing.py

Writes `group_shap_exp_vs_splicing.png` / `.pdf` into `--outdir`
(`analysis/group_shap/pictures/` by default).

This is the scatter Igor asked for in `../../chat_with_igor.txt` - "a scatter plot of feature
importance for splicing against expression" - recomputed from `gshap`. It is the counterpart of
`group_ablation/figures/exp_vs_splicing_groups.png` and keeps its panel structure exactly: 2 x 2,
rows are the task, columns are the table type, so all four combinations of the project's two axes
are on one sheet. x is the group's share of the `intron_exp` budget, y the same in
`intron_splicing`, both on log axes with the diagonal.

## No +-20 % band here, unlike the ablation twin

`band=False`. The stripe asserted a distinction finer than the data's own scatter: measured on this
run, the median interquartile range of a group's share across the 8 introns is **x1.51** while the
band is only **x1.44** wide, and **18 of the 49** points it called "leaning" have the two IQRs
overlapping. The three corner lines that quoted it - the tally and the two lean counts - go with it;
Spearman and the median ratio stay, and neither mentions the band.

Nothing else is lost, because in the GROUP figure the points are coloured by group identity. The
ablation twin keeps its band: its feature version colours points by the band's own verdict.

A claim that survives the removal is one whose whiskers separate. `nmd` is the case worth knowing:
intragenic regression `exp [4.29, 6.32] %` against `spl [11.17, 19.64] %`, classification
`[1.64, 2.46]` against `[11.62, 14.75]` - disjoint in both, so its lean to splicing rests on the
spread rather than on a threshold.

## What changed against the older figure

Only the quantity. The old one plots `sum_i mean|phi_i|` over the group's features - the modulus
taken BEFORE the sum, so members of one group pushing opposite ways on the same barcode add up
instead of cancelling. This one plots the exact group Shapley value with the group as ONE player,
all 2^k coalitions enumerated (`../PLAN.md` §2, quantity C). Same models, same 8 introns, same
round-0 feature sets; a group that moves between the two figures moves because of how its members
interact, not because anything about the model changed.

## Why the two datasets can share a point

The same 8 introns (A1, B1, B3-6, B3-9, C1, D1, E1, H1) and, per `../../CLAUDE.md` §6.2, the
`combined/`, `barcodes_order/` and `barcodes_bed/` files of an intron are byte-identical between
the two `DATA/` trees. The datasets differ in their READOUT and in nothing else, and
`data/group_map.json` confirms the two ran on the same membership - 59 features intragenic and 54
intergenic, identically grouped. So a point pairs two models that saw the same integrations and
the same columns and were asked two different questions.

## The one thing this figure cannot say

**A share is a slice of a budget, and both budgets sum to 1.** A point above the diagonal means the
splicing model spends a larger FRACTION of its attention there, not that it knows more about it.
The two models are not equally good - on holdout `intron_exp` regression reaches R2 ~ 0.42 against
`intron_splicing`'s ~ 0.09, and §6.2 shows the splicing readout's own two pools agree at only
r ~ 0.42, so there is a ceiling on it. Read this as composition, never as amount. The second,
smaller trap: `intergenic` is not a set of intergenic integrations, it is all the rows with the
eight positional features hidden, and five of `annotations`' ten members are among them.

## Why it is drawn by the other analysis's code

`draw()` is imported from `group_ablation/plotting/plot_exp_vs_splicing.py` and called unchanged -
it already takes the output name and the title as arguments, so nothing has to be overridden. The
figure will sit beside its ablation twin in the same document and the group colours have to be the
same colours; running the same function is the only way to guarantee that rather than assert it.
"""

import argparse
import contextlib
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from group_shap_shares import (ANALYSIS, ABLATION_PLOTTING, FIGURES, load_group_shap,   # noqa: E402
                               mean_group_share, spread)
sys.path.insert(0, ABLATION_PLOTTING)
sys.path.insert(0, ANALYSIS)
from group_labels import group_label, table_label, dataset_label, task_label            # noqa: E402

from common import MM                                                         # noqa: E402
import plot_exp_vs_splicing as exp_vs_splicing                                # noqa: E402
from plot_exp_vs_splicing import X_DATASET, Y_DATASET, draw                   # noqa: E402

TITLE = ("Expression against splicing · true group Shapley, share of the importance budget, "
         "one point = one feature group")
TITLE_REGRESSION = ("Expression against splicing · regression · true group Shapley, "
                    "share of the importance budget, one point = one feature group")
NAME = "group_shap_exp_vs_splicing"

# Two extra rows of seats, further from the marker than the borrowed placer's own +-22 pt. The
# group-Shapley budget is more concentrated than the summed-|SHAP| one the ablation figure plots -
# the top three groups hold about half of it - so more markers land in the same corner and the
# default seats run out. Run with the stock list, `global norm` in classification/intragenic and
# `mane expression` in classification/intergenic found no seat and were dropped, which on a
# a scatter of a dozen-odd colours leaves a point identifiable only by its hue.
#
# Appending is safe by construction and that is why it is an append and not a replacement: the
# placer takes the FIRST free seat and stops, so a label that already fitted somewhere in the stock
# list cannot be moved by anything added after it. The new seats are only ever reached by a label
# that was otherwise going to be dropped.
EXTRA_SEATS = [30, -30, 38, -38, 46, -46, 54, -54]

# The placer counts two labels as clashing when their boxes come within `gap` PIXELS, and the layout
# pass runs at the screen dpi, so the stock 2.0 is about one typographic point - enough to stop an
# overlap, not enough to stop two names reading as one word. Widened to 4.0, which is the smallest
# value that separates `gc content` from `nuclear org` in classification/intergenic.
LABEL_GAP = 4.0


@contextlib.contextmanager
def roomier_labels():
    """More seats and a wider clash margin, for the length of one `draw()` call.

    Both are properties of the borrowed placer rather than of the figure, so they are patched on the
    module and put back afterwards. Nothing else in this project imports that module at the same
    time; `run_ablation.py` and the ablation figure scripts each run in their own process.
    """
    seats, hits = exp_vs_splicing.SEATS, exp_vs_splicing._hits
    exp_vs_splicing.SEATS = list(seats) + EXTRA_SEATS
    exp_vs_splicing._hits = lambda box, others, gap=LABEL_GAP: hits(box, others, gap)
    try:
        yield
    finally:
        exp_vs_splicing.SEATS, exp_vs_splicing._hits = seats, hits


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=FIGURES)
    ap.add_argument("--regression-only", action="store_true",
                    help="one row instead of two, and both label variants — the paper pair")
    args = ap.parse_args(argv)

    frame = load_group_shap()
    common = dict(
        x_wide=mean_group_share(frame, X_DATASET), y_wide=mean_group_share(frame, Y_DATASET),
        x_spread=spread(frame, X_DATASET), y_spread=spread(frame, Y_DATASET),
        outdir=args.outdir, label_rule=lambda f: set(f.index),
        whiskers=True, band=False, median=False,
        namer=group_label, table_namer=table_label, task_namer=task_label,
        x_label=dataset_label(X_DATASET), y_label=dataset_label(Y_DATASET))

    if args.regression_only:
        # Two panels, so the sheet is half as tall; the width stays 183 mm because the panels are
        # square and `frame_panel` fixes the aspect - narrowing would shrink the panels, not the gap.
        for labels in (True, False):
            suffix = "" if labels else "_nolabels"
            with roomier_labels():
                draw("group", size=(183 * MM, 105 * MM),
                     name="%s_regression%s" % (NAME, suffix),
                     title=TITLE_REGRESSION, tasks=["regression"], labels=labels, **common)
        return 0

    with roomier_labels():
        frames = draw("group", size=(183 * MM, 183 * MM), name=NAME, title=TITLE, **common)

    for (task, table_type), panel in frames.items():
        both = panel[(panel["exp"] > 0) & (panel["spl"] > 0)]
        counts = panel["lean"].value_counts()
        print("%-14s %-11s n=%2d  rho=%.2f  median x%.2f  splicing %2d / even %2d / expression %2d"
              % (task, table_type, len(panel),
                 both["exp"].corr(both["spl"], method="spearman"), both["ratio"].median(),
                 counts.get("splicing", 0), counts.get("even", 0), counts.get("expression", 0)))

    # Distance from the diagonal, measured the way the axes measure it: the panels are log-log and
    # the band is a ratio, so |log2(ratio)| is the number of doublings a point sits off the line.
    # Ranking on the raw ratio instead would call x2.0 twice as far out as x0.5, which on this
    # figure is the same distance in opposite directions.
    print("\nfurthest from the diagonal, by |log2(splicing share / expression share)|")
    for (task, table_type), panel in frames.items():
        ordered = panel.dropna(subset=["ratio"])
        ordered = ordered.assign(off=np.abs(np.log2(ordered["ratio"]))).sort_values(
            "off", ascending=False)
        print("  %s / %s" % (task, table_type))
        for key, row in ordered.head(4).iterrows():
            print("    %-18s exp %6.3f %%   spl %6.3f %%   x%.2f  (%.2f doublings, %s)"
                  % (key, row["exp"] * 100, row["spl"] * 100, row["ratio"], row["off"],
                     row["lean"]))

    print("\nevery group, mean over the four panels")
    lr = pd.DataFrame({k: np.log2(v["ratio"]) for k, v in frames.items()})
    summary = pd.DataFrame({"mean_log2": lr.mean(axis=1), "distance": lr.abs().mean(axis=1),
                            "panels_leaning_to_splicing": (lr > 0).sum(axis=1)})
    summary = summary.sort_values("distance", ascending=False)
    print("  group                mean x   distance   panels to splicing")
    for group, row in summary.iterrows():
        print("  %-18s x%5.2f     %.2f       %d/4"
              % (group, 2 ** row["mean_log2"], row["distance"],
                 int(row["panels_leaning_to_splicing"])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
