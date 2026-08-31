"""
What the models rely on when the group is ONE PLAYER — the whole group Shapley budget.

    python analysis/group_shap/plotting/plot_group_shap_budget.py

Writes `group_shap_budget.png` / `.pdf` into `--outdir` (default `analysis/group_shap/figures/`).

One horizontal stacked bar per chain, twelve in all, each summing to 1. It is the deliberate
counterpart of `group_ablation/figures/group_budget.png` - same layout, same palette, same
averaging - so the two can be laid side by side. What differs is the quantity in the bar:

  * **there** a group's height is `sum over its features of mean|SHAP|`, the features being the
    players and the group only a label put on them afterwards.
  * **here** it is the exact group Shapley value, the group entered the game as a single player and
    all 2^k coalitions were enumerated - no sampling. `../PLAN.md` §4 is the method and
    §5 the check journal.

The two coincide only when nothing interacts across a group boundary, which is exactly the question
Igor asked (`chat_with_igor.txt`: "can you compute shap of removing them as a group? Summing up
shap values is not the same I think?"). He is right that they are different quantities; this figure
is what the difference looks like on the whole budget.

Three things are load-bearing, all three inherited from the figure this one mirrors:

  * the share is taken INSIDE one model and the mean over experiments divides by the NUMBER of
    experiments of the dataset, 2 / 8 / 8. `group_shap_shares.mean_group_share` asserts that count
    rather than trusting the row count.
  * the groups and their colours are `analysis/feature_groups.py`, unedited and imported, not
    copied - a group may not mean one thing here and another in the ablation figures.
  * `gshap` is `mean|phi_g|` and so non-negative; the share of a total is well defined. The signed
    `bshap` column of the same table is NOT what is normalised here.

**This is attention, not cost.** A group can hold a quarter of the budget and cost nothing when it
is removed, because a neighbour carries the same signal - that is `group_ablation`'s
`leave_one_out_costs.tsv`. And `global_norm` is technical: barcode abundance is a property of the
library, not of the locus.

The table also carries `gshap_se`, a bootstrap SE over the explained rows. It understates the true
noise - only 68 % of values fell within 2 SE across a seed replicate - so it is not drawn.
"""

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import matplotlib.pyplot as plt                                            # noqa: E402
from matplotlib.patches import Patch                                       # noqa: E402
from group_shap_shares import (DATASETS, FIGURES, FS, FS_SMALL, GRID, GROUP_COLORS, HAIR,  # noqa: E402
                          INK, INK_SOFT, MM, SURFACE, TABLE_TYPES, TASKS, ink_on,
                          load_gshap, mean_group_share, members_label, save)

ROW_H = 0.235       # inches of figure height per bar
BLOCK_GAP = 0.55    # extra rows of blank between two datasets
BAR_H = 0.66        # bar thickness, in row units


def collect(frame, datasets):
    """One record per bar, in drawing order."""
    bars = []
    for dataset in datasets:
        wide = mean_group_share(frame, dataset)
        for task in TASKS:
            for table_type in TABLE_TYPES:
                row = wide.loc[(task, table_type)]
                assert abs(float(row.sum()) - 1.0) < 1e-9, "budget sums to %r" % row.sum()
                bars.append(dict(dataset=dataset, task=task, table_type=table_type,
                                 shares=row[row > 0].to_dict()))
    return bars


def plot(bars, datasets, outdir, fixed_order=False):
    total = {}
    for bar in bars:
        for group, value in bar["shares"].items():
            total[group] = total.get(group, 0.0) + value
    # EACH BAR IS SORTED BY ITS OWN SHARES, biggest leftmost. The alternative - one order for all
    # twelve bars, by the total - lines the segments up so a column can be scanned, but then the
    # leading group of a given model is not at the start of its bar and the bar cannot be read on
    # its own. Reading each bar is what this figure is for; the colour key carries the identity.
    # `order` below stays global, and is used ONLY for the legend, so the key has a stable order.
    order = sorted(total, key=lambda g: -total[g])

    def bar_order(bar):
        # `fixed_order` draws the OTHER version of this figure: one order for all twelve bars, so a
        # colour sits at the same place everywhere and a column can be scanned. Both are kept and
        # both are published - they answer different questions, and neither replaces the other.
        if fixed_order:
            return order
        return sorted(bar["shares"], key=lambda g: -bar["shares"][g])

    y, cuts, cursor = [], [], 0.0
    for i, bar in enumerate(bars):
        if i and bar["dataset"] != bars[i - 1]["dataset"]:
            cuts.append(cursor + (BLOCK_GAP - 1.0) / 2.0)
            cursor += BLOCK_GAP
        y.append(cursor)
        cursor += 1.0
    span = cursor - 1.0 + BAR_H

    fig_h = ROW_H * span + 1.45
    fig = plt.figure(figsize=(183 * MM, fig_h))
    # 0.760 and not 0.775: the last x tick label is centred on the axis end, so half of "100" hangs
    # past it and the figure would no longer be the 183 mm it says it is.
    ax = fig.add_axes([0.215, 0.80 / fig_h, 0.760, 1 - 1.20 / fig_h])

    segments = []
    for bar, yy in zip(bars, y):
        left = 0.0
        for group in bar_order(bar):
            width = bar["shares"].get(group, 0.0)
            if width <= 0:
                continue
            ax.barh(yy, width, left=left, height=BAR_H, color=GROUP_COLORS[group],
                    edgecolor=SURFACE, linewidth=0.4, zorder=3)
            segments.append((group, left, width, yy))
            left += width

    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xticklabels(["0", "20", "40", "60", "80", "100"])
    ax.set_xlabel("share of the model's total group Shapley value, %, averaged over the experiments",
                  fontsize=FS)
    ax.set_ylim(max(y) + 0.75, min(y) - 0.75)
    ax.set_yticks(y)
    ax.set_yticklabels(["%s · %s" % (b["task"], b["table_type"]) for b in bars])
    ax.tick_params(direction="out", length=2.2, width=HAIR, colors=INK, labelsize=FS, pad=2)
    ax.tick_params(axis="y", length=0, colors=INK_SOFT, labelsize=FS_SMALL)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_linewidth(HAIR)
    ax.spines["bottom"].set_color(INK)
    ax.grid(axis="x", color=GRID, linewidth=HAIR, zorder=0)
    ax.grid(axis="y", visible=False)
    ax.set_axisbelow(True)

    for cut in cuts:
        ax.axhline(cut, xmin=-0.28, xmax=1.0, color=GRID, linewidth=HAIR, zorder=1, clip_on=False)
    for dataset in datasets:
        block = [yy for yy, b in zip(y, bars) if b["dataset"] == dataset]
        ax.text(-0.235, float(np.mean(block)), "%s" % dataset,
                transform=ax.get_yaxis_transform(), rotation=90, ha="center", va="center",
                fontsize=FS, color=INK)

    # In-segment labels, placed from MEASURED text widths rather than a guess at character width.
    # A 7 pt name needs about 12 mm and the bars are 142 mm wide, so a group under ~8 % of the
    # budget cannot carry its own name here and the key below carries it instead.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    axes_w_in = ax.get_window_extent(renderer).width / fig.dpi
    for group, left, width, yy in segments:
        for text in ("%s %.0f" % (group.replace("_", " "), width * 100),
                     group.replace("_", " ")):
            probe = fig.text(0, 0, text, fontsize=FS_SMALL - 0.5)
            need = probe.get_window_extent(renderer).width / fig.dpi + 0.055
            probe.remove()
            if width * axes_w_in >= need:
                ax.text(left + width / 2.0, yy, text, ha="center", va="center",
                        fontsize=FS_SMALL - 0.5, color=ink_on(GROUP_COLORS[group]), zorder=4)
                break

    # The one colour key in this figure set. A stacked bar cannot name its own thin segments at
    # this size - the smallest groups are one to three per cent of the budget, about 2 mm - and an
    # unnamed colour is worse than a key.
    handles = [Patch(facecolor=GROUP_COLORS[g], edgecolor=SURFACE, linewidth=0.4) for g in order]
    key = fig.legend(handles, [members_label(g) for g in order], loc="lower center",
                     bbox_to_anchor=(0.5, 0.004), ncol=7, frameon=False, handlelength=1.0,
                     handleheight=1.0, columnspacing=1.2, handletextpad=0.4, fontsize=FS_SMALL)
    for text, group in zip(key.get_texts(), order):
        text.set_color(GROUP_COLORS[group])

    save(fig, "group_shap_budget" + ("_fixed_order" if fixed_order else ""), outdir,
         "What the models rely on · group Shapley budget, the group as one player, "
         "every bar sums to 100 %"
         + (" · one order for every bar" if fixed_order else " · each bar sorted by its own shares"))
    plt.close(fig)
    return order


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=FIGURES)
    ap.add_argument("--datasets", nargs="+", default=DATASETS, choices=DATASETS)
    ap.add_argument("--fixed-order", action="store_true", dest="fixed_order",
                    help="one stacking order for every bar instead of sorting each by its own "
                         "shares; writes group_shap_budget_fixed_order.*")
    args = ap.parse_args(argv)

    frame = load_gshap()
    bars = collect(frame, args.datasets)
    order = plot(bars, args.datasets, args.outdir, fixed_order=args.fixed_order)

    print("\ngroup Shapley budget per group, one line per bar")
    for bar in bars:
        top = sorted(bar["shares"].items(), key=lambda kv: -kv[1])[:4]
        print("  %-16s %-15s %-11s sum=%.9f  %s"
              % (bar["dataset"], bar["task"], bar["table_type"], sum(bar["shares"].values()),
                 "  ".join("%s %.1f%%" % (g, v * 100) for g, v in top)))
    print("\nmean over the four chains of a dataset, biggest three:")
    for dataset in args.datasets:
        part = [b for b in bars if b["dataset"] == dataset]
        mean = {g: sum(b["shares"].get(g, 0.0) for b in part) / len(part) for g in order}
        best = sorted(mean.items(), key=lambda kv: -kv[1])[:3]
        print("  %-16s %s" % (dataset, ", ".join("%s %.1f%%" % (g, v * 100) for g, v in best)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
