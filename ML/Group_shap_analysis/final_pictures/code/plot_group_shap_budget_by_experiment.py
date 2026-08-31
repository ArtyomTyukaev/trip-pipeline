"""
The group Shapley budget of every single experiment — one bar per experiment.

    python analysis/group_shap/plotting/plot_group_shap_budget_by_experiment.py

Writes, into `--outdir` (default `analysis/group_shap/figures/`):

    group_shap_budget_by_experiment_<task>_<table_type>.{png,pdf}   — four files, both tasks x
                                                                     both table types

`plot_group_shap_budget.py` draws the same quantity averaged over the experiments of a dataset,
which answers "what does this dataset rely on". This one does not average: every experiment gets
its own bar, all three datasets on one page, so the question becomes "do the experiments of a
dataset actually agree, or does the mean hide a split". A dataset whose bars all look alike is one
where the average means something; a dataset whose bars differ is one where it does not.

It is the counterpart of `group_ablation/figures/group_budget_by_experiment_<task>_<table>.png`,
drawn in the same layout and the same palette so the two can be read against each other. The
quantity differs: there a group's height is the sum of `mean|SHAP|` over its features, here it is
the exact group Shapley value with the group as a single player, all coalitions enumerated
(`../PLAN.md` §4).

**This is attention, not cost.** A group can hold a quarter of the budget and cost nothing when it
is removed, because a neighbouring group carries the same signal. The cost lives in
`group_ablation/<dataset>/data/leave_one_out_costs.tsv`.
"""

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import matplotlib.pyplot as plt                                            # noqa: E402
from matplotlib.patches import Patch                                       # noqa: E402
from group_shap_shares import (ANALYSIS, DATASETS, FIGURES, FS, FS_SMALL,  # noqa: E402
                          GROUP_COLORS, INK, MM, SURFACE, ink_on, journal, load_gshap,
                          per_experiment_group_share, save)

sys.path.insert(0, ANALYSIS)
from experiment_labels import label, order                                  # noqa: E402

LABEL_FLOOR = 0.055        # a segment narrower than this cannot hold its own number


def draw(frame, task, table_type, outdir, fixed_order=False):
    """One stacked bar per experiment, the three datasets separated by a gap.

    `fixed_order` picks between the two versions of this figure, and both are published:

      * default - each bar sorted by ITS OWN shares, biggest at the bottom, so the leading group
        of an experiment is at the foot of its own bar and the bar reads on its own;
      * `fixed_order` - one order for every column, by the mean across them, so a colour sits at
        the same height everywhere and a row can be scanned across experiments.

    Neither replaces the other: the first answers "what leads this experiment", the second "how
    does this group move across experiments".
    """
    bars, blocks, x = [], [], 0
    for dataset in DATASETS:
        try:
            wide = per_experiment_group_share(frame, dataset)
            part = wide.xs((task, table_type), level=("task", "table_type"))
        except KeyError:
            continue
        first = x
        # `order`, not `sorted`: the introns run B3-9, B3-6, H1, E1, A1, B1, C1, D1, and
        # single_combined's two names fall through it unchanged. See analysis/experiment_labels.py.
        for exp in order(part.index):
            bars.append((dataset, exp, part.loc[exp]))
            x += 1
        blocks.append((dataset, first, x - 1))
        x += 1                                       # the gap between datasets

    if not bars:
        print("  nothing for %s / %s" % (task, table_type))
        return

    # `groups` is the global order, by the mean share across all bars. It always drives the legend,
    # so the key is stable, and it drives the stacking too when `fixed_order` is on.
    groups = sorted({g for _, _, s in bars for g in s.index},
                    key=lambda g: -float(np.mean([s.get(g, 0.0) for _, _, s in bars])))

    def bar_order(share):
        present = [g for g in share.index if float(share.get(g, 0.0)) > 0]
        if fixed_order:
            return [g for g in groups if g in present]
        return sorted(present, key=lambda g: -float(share[g]))

    fig, ax = plt.subplots(figsize=(183 * MM, 4.0))
    # bottom was 0.185, sized for a rotated `B3-6`. The gene names are up to `PACRGL Control`,
    # 14 characters at 6 pt = about 0.75 in of vertical room against a 4.0 in canvas, so the
    # margin has to hold 0.19 of the figure for the text alone before any tick or padding.
    fig.subplots_adjust(left=0.055, right=0.80, top=0.86, bottom=0.30)

    # x of every bar, matching the order `bars` was built in, gaps between dataset blocks
    pos = []
    for dataset, first, last in blocks:
        for i in range(last - first + 1):
            pos.append(first + i)
    for (dataset, exp, share), p in zip(bars, pos):
        bottom = 0.0
        for g in bar_order(share):
            v = float(share.get(g, 0.0))
            if v <= 0:
                continue
            ax.bar(p, v, bottom=bottom, width=0.80, color=GROUP_COLORS.get(g, INK),
                   edgecolor=SURFACE, linewidth=0.6, zorder=2)
            if v >= LABEL_FLOOR:
                ax.text(p, bottom + v / 2, "%.0f" % (v * 100), ha="center", va="center",
                        fontsize=FS_SMALL - 1.5, color=ink_on(GROUP_COLORS.get(g, INK)), zorder=4)
            bottom += v

    ax.set_xticks(pos)
    # Gene names on the axis, intron codes everywhere underneath. `PACRGL Control` is three times
    # the width of `B3-9` and these labels are rotated, so the room they need is vertical and
    # `bottom` above was widened for it.
    ax.set_xticklabels([label(e) for _d, e, _s in bars], fontsize=FS_SMALL, rotation=90)
    ax.set_xlim(-0.7, max(pos) + 0.7)
    ax.set_ylim(0, 1.0)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0", "25", "50", "75", "100 %"], fontsize=FS_SMALL)
    ax.set_ylabel("share of the group Shapley budget", fontsize=FS_SMALL)
    journal(ax, ygrid=True)
    for dataset, first, last in blocks:
        ax.annotate(dataset, xy=((first + last) / 2, 1.0),
                    xycoords=("data", "axes fraction"), xytext=(0, 4),
                    textcoords="offset points", ha="center", va="bottom",
                    fontsize=FS, color=INK)

    # A key, because the thinnest groups are a couple of millimetres of bar and cannot carry their
    # own number. Every segment wide enough is labelled in place; the key names the rest.
    handles = [Patch(facecolor=GROUP_COLORS.get(g, INK), linewidth=0.0,
                     label=g.replace("_", " ")) for g in groups]
    legend = ax.legend(handles=handles, fontsize=FS_SMALL - 1, loc="upper left",
                       bbox_to_anchor=(1.005, 1.02), handlelength=1.0, labelspacing=0.30,
                       frameon=False)
    for text, g in zip(legend.get_texts(), groups):
        text.set_color(GROUP_COLORS.get(g, INK))

    name = "group_shap_budget_by_experiment_%s_%s%s" % (
        task, table_type, "_fixed_order" if fixed_order else "")
    save(fig, name, outdir,
         title="What each experiment relies on · group Shapley · %s · %s   ·   "
               "one bar = one experiment, %% of that model's budget" % (task, table_type))
    plt.close(fig)
    print("    %d bars over %d datasets, %d groups" % (len(bars), len(blocks), len(groups)))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=FIGURES)
    ap.add_argument("--tasks", nargs="+", default=["regression", "classification"],
                    choices=["regression", "classification"])
    ap.add_argument("--fixed-order", action="store_true", dest="fixed_order",
                    help="one stacking order for every column instead of sorting each bar by its "
                         "own shares; writes *_fixed_order.png/.pdf alongside, not instead")
    args = ap.parse_args(argv)
    frame = load_gshap()
    for task in args.tasks:
        for table_type in ("intragenic", "intergenic"):
            print("%s · %s" % (task, table_type))
            draw(frame, task, table_type, args.outdir, fixed_order=args.fixed_order)
    return 0


if __name__ == "__main__":
    sys.exit(main())
