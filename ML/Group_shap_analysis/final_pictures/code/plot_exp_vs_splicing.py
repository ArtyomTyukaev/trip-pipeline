"""
Expression against splicing — the same importance budget, split two ways. Round 0.

    python analysis/group_ablation/plotting/plot_exp_vs_splicing.py

Writes into `--outdir`:

  exp_vs_splicing_groups.png     one point per feature GROUP,   13 points per panel
  exp_vs_splicing_features.png   one point per raw FEATURE,     61 points per panel

Both are 2 x 2: rows are the task, columns are the table type, so all four combinations of the
project's two axes are on one sheet. x is the group's / feature's share of the round-0 importance
budget in `intron_exp`, y the same quantity in `intron_splicing`.

This is the scatter Igor asked for in `chat_with_igor.txt` - "feature importance for splicing
against expression".

## Why the two datasets can be put on the same point

They are the same 8 introns (A1, B1, B3-6, B3-9, C1, D1, E1, H1) and, per `../../CLAUDE.md` §6.2,
the `combined/`, `barcodes_order/` and `barcodes_bed/` files of an intron are byte-identical between
the two `DATA/` trees. The two datasets differ in their READOUT and in nothing else. Verified on the
tables actually trained: the round-0 feature sets are identical in all four panels, 59 features
intragenic and 54 intergenic.

So a point pairs two models that saw the same integrations and the same columns and were asked two
different questions. That is as clean a comparison as this project has.

## Why round 0 and only round 0

Shares are renormalised over the survivors every round, so a group that lives to round 5 has been
multiplied up by five removals and a group dropped in round 1 never leaves its round-0 value.
`../../CLAUDE.md` §8.3 measures the inversion this causes. Round 0 is the only round in which every
group is present in both datasets, and therefore the only round whose two budgets are commensurable.

## The one thing this figure cannot say

**A share is a slice of a budget, and both budgets sum to 1.** A point above the diagonal means the
splicing model spends a LARGER FRACTION of its attention there - not that it knows more about it.
The two models are not equally good: on holdout, `intron_exp` regression reaches R2 ~ 0.42 against
`intron_splicing`'s ~ 0.09, and §6.2 shows the splicing readout's own two pools agree at only
r ~ 0.42, so there is a ceiling on it. Read this figure as **composition**, never as amount.

The second, smaller trap: `intergenic` is not a set of intergenic integrations, it is all the rows
with the eight positional features hidden, and five of `annotations`' ten members are among them.

## The +-20 % band

`BAND = 1.2`. A point is coloured only when its two shares differ by more than a fifth, that is when
`y/x >= 1.2` or `<= 1/1.2`; inside that band it stays grey. The band is a ratio, so it is a constant
width on the log axes both figures use - the shares run over two and a half orders of magnitude and
a linear axis would put ten of the thirteen groups in one corner.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import matplotlib.pyplot as plt                                             # noqa: E402
from matplotlib.lines import Line2D                                         # noqa: E402
from common import (FIGURES, FS, FS_SMALL, GRID, GROUP_COLORS, HAIR, INK, INK_SOFT,
                    MM, MUTED, SURFACE, TABLE_TYPES, TASKS)                 # noqa: E402
from group_shares import N_EXPERIMENTS, load_shap, mean_group_share, save   # noqa: E402

X_DATASET, Y_DATASET = "intron_exp", "intron_splicing"
X_LABEL, Y_LABEL = "intron_exp  ·  expression", "intron_splicing  ·  splicing"

# A fifth. Everything inside is "the two readouts weight it the same".
BAND = 1.2

# Which way a point leans. Neither colour may be the reserved green/violet of the table types
# (`group_shares.TABLE_COLOUR`) - the columns of this figure ARE the table types, and a point
# coloured green inside an intragenic panel would read as a third statement about the table.
LEAN = {"expression": "#b8500f", "splicing": "#0f6f8c", "even": "#9b9b96"}


# ------------------------------------------------------------------------------- data
def mean_feature_share(dataset):
    """(task, table_type) x feature, round 0, levels folded up to the raw feature.

    The same three rules `group_shares.mean_group_share` follows, one fold shallower: `share`
    already sums to 1 over the LEVELS of one model, so summing it over the levels of a feature is
    the whole of the fold; and the average over experiments divides by `N_EXPERIMENTS`, so a feature
    a model never used contributes a zero rather than dropping out of its own denominator.

    The fold is done **per panel** and not once over the whole table, which is not a detail. The
    intragenic and the intergenic models do not hold the same columns - `MANE_all` and `MANE_exonic`
    are on the intragenic drop list and five more annotation features are on the intergenic one - so
    one unstack over both would invent a zero-importance row for a feature that was never in the
    model. Zero and absent are different things here: a zero would be plotted on the axis and
    counted in the band tally. Absent comes back NaN and `paired()` drops it.
    """
    shap = load_shap(dataset)
    shap = shap[shap["round"] == 0]
    panels = {}
    for (task, table_type), part in shap.groupby(["task", "table_type"]):
        per = part.groupby(["experiment", "feature"])["share"].sum().unstack("feature").fillna(0.0)
        total = per.sum(axis=1).to_numpy()
        assert np.allclose(total, 1.0, atol=1e-9), "a model's feature shares sum to %r" % total[
            np.argmax(np.abs(total - 1.0))]
        panels[(task, table_type)] = per.sum() / N_EXPERIMENTS[dataset]
    return pd.DataFrame(panels).T


def mean_group_share_round0(dataset):
    """(task, table_type) x group, round 0, off the shared helper."""
    wide = mean_group_share(load_shap(dataset), dataset, key="round")
    wide = wide.loc[0]
    return wide


def spread(dataset, unit):
    """(task, table_type, unit) -> the 8 per-intron shares, for the whiskers of the group figure."""
    shap = load_shap(dataset)
    shap = shap[shap["round"] == 0]
    return (shap.groupby(["task", "table_type", "experiment", unit])["share"]
            .sum().unstack(unit).fillna(0.0))


def paired(x_wide, y_wide, task, table_type):
    """One frame of the two shares, the ratio, and which side the point leans to.

    `dropna()` keeps only what BOTH datasets had in this panel. Verified on the tables as built:
    that is all 59 features intragenic and all 54 intergenic, so it drops nothing here - it is the
    assertion that the two datasets really are the same columns, made where it would fail loudly.
    """
    x, y = x_wide.loc[(task, table_type)], y_wide.loc[(task, table_type)]
    frame = pd.DataFrame({"exp": x, "spl": y}).dropna()
    exp, spl = frame["exp"].to_numpy(), frame["spl"].to_numpy()

    with np.errstate(divide="ignore", invalid="ignore"):
        frame["ratio"] = np.where(exp > 0, spl / np.where(exp > 0, exp, np.nan), np.nan)
    # A zero share is a real answer - the model had the feature and gave it nothing - so it gets a
    # lean rather than being dropped; only 0 against 0 has no direction and stays `even`.
    frame["lean"] = np.select(
        [(exp <= 0) & (spl <= 0),
         exp <= 0,
         spl <= 0,
         frame["ratio"] >= BAND,
         frame["ratio"] <= 1 / BAND],
        ["even", "splicing", "expression", "splicing", "expression"], default="even")
    return frame.sort_values("exp", ascending=False)


def axis_floor(frames):
    """The bottom of the log axes: below the smallest share anyone actually has, on a round decade
    step. An axis that starts two decades under the data spends half the panel on nothing."""
    positive = np.concatenate([f[["exp", "spl"]].to_numpy().ravel() for f in frames])
    positive = positive[positive > 0]
    half_decade = 10 ** (np.floor(np.log10(positive.min()) * 2) / 2)
    # ...but never more than three and a half decades of it. On the feature figure the smallest
    # share is 1e-5, and an axis reaching that spends two thirds of its length on features that
    # hold a thousandth of a per cent of the budget. Those are clipped onto the floor and drawn
    # hollow instead.
    return max(half_decade, 3e-4)


def display_name(key):
    """The name as it goes on the panel.

    One rule, and it is reversible: the four tracks of GSE81553 carry the accession and the cell
    line in front of the track name, which is 14 characters of prefix that is the same on all four
    and says nothing a reader of this figure needs. `GSE81553_K562_GRO-Seq_binned_20kb` at full
    length does not fit beside its marker, and it is one of the two strongest
    expression-leaning features in the panel - dropping ITS label to keep a prefix would be the
    wrong trade. Nothing else is renamed; underscores become spaces.
    """
    return str(key).replace("GSE81553_K562_", "").replace("_", " ")


# ------------------------------------------------------------------------------- labels
# The vertical offsets a label may sit at, in typographic points, tried in this order: level with
# the marker first, then progressively further above and below. Each is tried on the right of the
# marker and then on its left, which makes 22 seats.
SEATS = [0, 8, -8, 15, -15, 22, -22]
# Past this a label is far enough from its marker to be read as belonging to a neighbour, so it
# grows a leader. Under it the label sits against the marker and a leader would only add ink.
LEADER_AT = 13


def _box(anchor, size, dx, dy, ha, px):
    """Pixel bounding box of a label placed at (dx, dy) points from `anchor`."""
    x = anchor[0] + dx * px
    x0 = x if ha == "left" else x - size[0]
    y0 = anchor[1] + dy * px - size[1] / 2.0
    return np.array([x0, y0, x0 + size[0], y0 + size[1]])


def _hits(box, others, gap=2.0):
    for other in others:
        if (box[0] < other[2] + gap and other[0] < box[2] + gap
                and box[1] < other[3] + gap and other[1] < box[3] + gap):
            return True
    return False


def place_labels(ax, xs, ys, texts, colours, fontsize, order, occupied=(), pad=4.5):
    """Greedy label placement: strongest point first, each takes the nearest free seat.

    A force-directed repulsion was tried first and is the wrong tool for this scatter. Sixty
    features sit on one diagonal, so pushing overlapping labels apart resolves the cloud into a
    tall column of names joined to their points by long crossing leaders - legible as a list,
    useless as a scatter.

    Greedy placement inverts the failure mode: a label either sits AGAINST its own marker or is not
    drawn at all. Nothing is ever attached to the wrong point, no leader is needed, and the labels
    that get dropped are by construction the least important ones - `order` runs from the largest
    share down, so the names a reader looks for are placed while the panel is still empty. The
    count of what was dropped is returned, never swallowed.

    `xs`/`ys` are EVERY point, because every marker is an obstacle; `order` lists only the ones to
    be labelled, in the sequence they get first refusal. `occupied` seeds the placer with boxes it
    must avoid - the panel's own corner block of text.
    """
    fig = ax.figure
    if len(order) == 0:
        return []
    px = fig.dpi / 72.0
    # `arrowprops` only ever takes effect if it is passed to the constructor - the arrow patch is
    # built in `__init__` and never afterwards - so every label is born with a leader and the ones
    # that end up sitting against their marker have theirs hidden.
    anns = {i: ax.annotate(texts[i], (xs[i], ys[i]), xytext=(0, 0), textcoords="offset points",
                           fontsize=fontsize, color=colours[i], va="center",
                           annotation_clip=False, zorder=6,
                           arrowprops=dict(arrowstyle="-", color=MUTED, linewidth=HAIR * 0.9,
                                           shrinkA=0.0, shrinkB=3.0))
            for i in order}
    fig.canvas.draw()
    sizes = {i: (a.get_window_extent().width, a.get_window_extent().height)
             for i, a in anns.items()}
    anchors = ax.transData.transform(np.column_stack([xs, ys]))
    frame = np.array(ax.get_window_extent().extents)          # x0, y0, x1, y1

    placed = [np.array([a[0] - 5.0, a[1] - 5.0, a[0] + 5.0, a[1] + 5.0]) for a in anchors]
    placed += [np.asarray(b) for b in occupied]
    dropped = []

    for i in order:
        best = None
        for dy in SEATS:
            for ha in ("left", "right"):
                dx = pad if ha == "left" else -pad
                box = _box(anchors[i], sizes[i], dx, dy, ha, px)
                if box[0] < frame[0] or box[2] > frame[2] or box[1] < frame[1] or box[3] > frame[3]:
                    continue
                if _hits(box, placed):
                    continue
                best = (dx, dy, ha, box)
                break
            if best is not None:
                break
        if best is None:
            dropped.append(texts[i])
            anns[i].remove()
            continue
        dx, dy, ha, box = best
        anns[i].set_ha(ha)
        anns[i].xyann = (dx, dy)
        anns[i].arrow_patch.set_visible(abs(dy) >= LEADER_AT)
        placed.append(box)
    return dropped


# ------------------------------------------------------------------------------- drawing
def frame_panel(ax, lo, hi, ticks, band=True):
    """Log-log square, the diagonal, and - when `band` - the +-20 % stripe around it.

    `band=False` keeps the diagonal and drops the stripe. See `draw()` for why a caller would.
    """
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")

    edge = np.geomspace(lo, hi, 60)
    if band:
        ax.fill_between(edge, edge / BAND, edge * BAND, color="#eceae4", linewidth=0, zorder=0)
    ax.plot(edge, edge, color=MUTED, linewidth=HAIR, zorder=1)
    if band:
        for factor in (BAND, 1 / BAND):
            ax.plot(edge, edge * factor, color=MUTED, linewidth=HAIR * 0.8,
                    linestyle=(0, (2.5, 2.0)), zorder=1)

    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(["%g%%" % (t * 100) for t in ticks], fontsize=FS_SMALL)
    ax.set_yticklabels(["%g%%" % (t * 100) for t in ticks], fontsize=FS_SMALL)
    ax.minorticks_off()
    ax.grid(True, color=GRID, linewidth=HAIR * 0.7, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(HAIR)
        ax.spines[side].set_color(INK)
    ax.tick_params(direction="out", length=2.2, width=HAIR, colors=INK, pad=2)


def panel_corner(ax, frame, unit, band=True, median=True):
    """Spearman of the two budgets, the tally of the band, and the tilt. Bottom-right, inside.

    With `band=False` the three band-dependent lines go too - the tally and the two lean counts are
    statements ABOUT the +-20 % threshold and would be quoting a line the panel no longer draws.
    Spearman and the median ratio stay: neither mentions the band, and the median is the one that
    keeps the figure honest, for the reason in the next paragraph.

    The tilt line is not decoration. Both budgets sum to 1, so they cannot both be concentrated:
    the expression models put a third of their attention on their top three features and the
    splicing models a fifth, and the arithmetic of that alone lifts everything else on the splicing
    side. The median ratio measures the lift, and it is the number that stops "40 of 61 lean to
    splicing" being read as forty separate findings.

    Returns the pixel boxes it wrote into, so the label placer treats them as occupied.
    """
    both = frame[(frame["exp"] > 0) & (frame["spl"] > 0)]
    rho = both["exp"].corr(both["spl"], method="spearman")
    counts = frame["lean"].value_counts()
    n_total = len(frame)
    lines = [("Spearman ρ = %.2f" % rho, INK_SOFT)]
    if median:
        lines.append(("median %s at ×%.2f" % (unit, both["ratio"].median()), INK_SOFT))
    if band:
        lines += [
            ("%d of %d outside the band" % (n_total - counts.get("even", 0), n_total), INK_SOFT),
            ("%d lean to splicing" % counts.get("splicing", 0), LEAN["splicing"]),
            ("%d lean to expression" % counts.get("expression", 0), LEAN["expression"]),
        ]
    anns = []
    for i, (text, colour) in enumerate(lines):
        anns.append(ax.annotate(
            text, xy=(0.975, 0.030 + (len(lines) - 1 - i) * 0.048), xycoords="axes fraction",
            ha="right", va="bottom", fontsize=FS_SMALL - 0.5, color=colour))
    ax.figure.canvas.draw()
    return [ann.get_window_extent().extents for ann in anns]


def draw(unit, x_wide, y_wide, x_spread, y_spread, outdir, label_rule, whiskers, size, name,
         title, label_note="", band=True, median=True, tasks=None, labels=True):
    """`band=False` drops the +-20 % stripe and the three corner lines that quote it.

    Three more switches, all defaulting to the original behaviour so the ablation twin is untouched:

    * `median=False` drops the "median group at x…" line from the corner, leaving Spearman. The
      median is the honest reference line for judging ONE group - both budgets sum to 1, so the
      whole cloud sits off the diagonal by construction - but stating it without drawing it invited
      the opposite reading, that a point above the DIAGONAL is leaning. Off for the paper figures.
    * `tasks=[…]` draws only those task rows. One row is a legal grid here.
    * `labels=False` skips the label placer, leaving bare coloured points.

    Why a caller would drop the band: it asserts a distinction finer than the data's own scatter. Measured
    on the group-Shapley run, the median interquartile range of a group's share across the 8 introns
    is x1.51 while the band is x1.44 wide, and 18 of the 49 points it calls "leaning" have their two
    IQRs overlapping - the introns disagree enough that the lean is not supported.

    It is a flag rather than a deletion because `unit == "feature"` COLOURS its points by the band's
    verdict (`LEAN`), so dropping the band there would leave the colours quoting a threshold the
    panel no longer shows. Only the group figure, whose points are coloured by group identity, can
    lose the band without losing anything else.
    """
    tasks = list(TASKS) if tasks is None else list(tasks)
    frames = {(task, table_type): paired(x_wide, y_wide, task, table_type)
              for task in tasks for table_type in TABLE_TYPES}
    lo = axis_floor(frames.values())
    hi = 0.45
    ticks = [t for t in (0.0001, 0.001, 0.01, 0.1) if lo <= t < hi]
    # squeeze=False: with one task row `axes` would come back 1-D and `axes[r][c]` would index
    # the wrong axis instead of failing.
    fig, axes = plt.subplots(len(tasks), len(TABLE_TYPES), figsize=size, squeeze=False)
    dropped_total = []

    for r, task in enumerate(tasks):
        for c, table_type in enumerate(TABLE_TYPES):
            ax = axes[r][c]
            frame = frames[(task, table_type)]
            frame_panel(ax, lo, hi, ticks, band=band)

            xs = frame["exp"].clip(lower=lo).to_numpy()
            ys = frame["spl"].clip(lower=lo).to_numpy()
            if unit == "group":
                colours = [GROUP_COLORS.get(g, INK) for g in frame.index]
            else:
                colours = [LEAN[v] for v in frame["lean"]]

            if whiskers:
                # The two datasets are the same eight introns, so the whisker is the interquartile
                # range OF THE SAME EIGHT MODELS' opinions - a statement about how much the introns
                # disagree, not about a fit.
                for i, key in enumerate(frame.index):
                    for spread_frame, along in ((x_spread, 0), (y_spread, 1)):
                        q1, q3 = np.quantile(spread_frame.loc[(task, table_type)][key], [0.25, 0.75])
                        q1, q3 = max(q1, lo), max(q3, lo)
                        seg = ([q1, q3], [ys[i], ys[i]]) if along == 0 else ([xs[i], xs[i]], [q1, q3])
                        ax.plot(*seg, color=colours[i], linewidth=HAIR * 1.6, alpha=0.42,
                                solid_capstyle="butt", zorder=2)

            # A share at or under the floor - zero included - cannot be placed on the log axis. It
            # is drawn hollow ON the floor, which says "somewhere at or below here", never a value.
            hollow = (frame[["exp", "spl"]] <= lo).any(axis=1).to_numpy()
            ax.scatter(xs[~hollow], ys[~hollow], s=26 if unit == "group" else 15,
                       c=[colours[i] for i in np.where(~hollow)[0]],
                       edgecolor=SURFACE, linewidth=0.6, zorder=4)
            if hollow.any():
                ax.scatter(xs[hollow], ys[hollow], s=22, facecolor="none", zorder=4,
                           edgecolor=[colours[i] for i in np.where(hollow)[0]], linewidth=0.8)

            corner = panel_corner(ax, frame, unit, band=band, median=median)
            if labels:
                keep = label_rule(frame)
                idx = [i for i, key in enumerate(frame.index) if key in keep]
                rank = sorted(idx, key=lambda i: -max(xs[i], ys[i]))
                dropped_total += place_labels(
                    ax, xs, ys, [display_name(k) for k in frame.index], colours,
                    fontsize=FS_SMALL - (0 if unit == "group" else 0.8), order=rank,
                    occupied=corner)

            if r == 0:
                ax.annotate(table_type, xy=(0.5, 1.045), xycoords="axes fraction", ha="center",
                            va="bottom", fontsize=FS + 0.5, color=INK)
            if c == 0:
                ax.annotate(task, xy=(-0.235, 0.5), xycoords="axes fraction", ha="center",
                            va="center", rotation=90, fontsize=FS + 0.5, color=INK)
                ax.set_ylabel("share of the importance budget\n%s" % Y_LABEL, fontsize=FS_SMALL,
                              linespacing=1.35)
            if r == len(tasks) - 1:
                ax.set_xlabel("share of the importance budget\n%s" % X_LABEL, fontsize=FS_SMALL,
                              linespacing=1.35)

    any_hollow = any((f[["exp", "spl"]] <= lo).any(axis=1).any() for f in frames.values())
    handles = []
    if unit == "feature":
        handles += [Line2D([], [], marker="o", linestyle="none", markersize=4, color=LEAN[k],
                           label=v) for k, v in
                    (("expression", "≥ 20 % more of the budget in intron_exp"),
                     ("even", "within ±20 % — the two readouts weight it the same"),
                     ("splicing", "≥ 20 % more of the budget in intron_splicing"))]
    else:
        handles.append(Line2D([], [], linestyle="-", color=MUTED, linewidth=HAIR * 1.6, alpha=0.55,
                              label="whiskers: interquartile range over the 8 introns"))
    if unit == "feature":
        handles.append(Line2D([], [], linestyle="none", label=label_note))
    if any_hollow:
        handles.append(Line2D([], [], marker="o", linestyle="none", markersize=4,
                              markerfacecolor="none", markeredgecolor=MUTED,
                              label="at or below %g %% — drawn on the floor, value not readable"
                              % (lo * 100)))
    fig.legend(handles=handles, loc="lower center", ncol=2 if len(handles) > 2 else 1,
               fontsize=FS_SMALL, bbox_to_anchor=(0.5, -0.005), handletextpad=0.5,
               columnspacing=1.6, labelcolor=INK_SOFT)

    fig.subplots_adjust(hspace=0.28, wspace=0.20)
    save(fig, name, outdir, title)
    plt.close(fig)
    if dropped_total:
        # Never silent. A reader who cannot find a feature on the panel has to be able to learn
        # from the run whether it was absent or only unlabelled.
        print("  [!] %d label(s) had no free seat and were not drawn: %s"
              % (len(dropped_total), ", ".join(sorted(set(dropped_total)))))
    return frames


# ------------------------------------------------------------------------------- entry
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default=FIGURES)
    ap.add_argument("--label-floor", type=float, default=0.02,
                    help="on the FEATURE figure, label a feature only if it holds at least this "
                         "share of the budget in one of the two datasets (default 2 %%)")
    args = ap.parse_args(argv)

    groups = draw(
        "group",
        mean_group_share_round0(X_DATASET), mean_group_share_round0(Y_DATASET),
        spread(X_DATASET, "group"), spread(Y_DATASET, "group"),
        args.outdir,
        label_rule=lambda frame: set(frame.index),
        whiskers=True, size=(183 * MM, 183 * MM), name="exp_vs_splicing_groups",
        title="Expression against splicing · share of the round-0 importance budget, "
              "one point = one feature group")

    features = draw(
        "feature",
        mean_feature_share(X_DATASET), mean_feature_share(Y_DATASET),
        spread(X_DATASET, "feature"), spread(Y_DATASET, "feature"),
        args.outdir,
        # Labelled by SIZE, not by lean: the colour already carries the lean, and labelling every
        # leaning feature means labelling three quarters of the panel.
        label_rule=lambda frame: set(frame.index[frame[["exp", "spl"]].max(axis=1)
                                                 >= args.label_floor]),
        whiskers=False, size=(215 * MM, 215 * MM), name="exp_vs_splicing_features",
        title="Expression against splicing · share of the round-0 importance budget, "
              "one point = one feature",
        label_note="labelled: every feature holding ≥ %g %% of the budget in one of the two, where "
                   "the name fits.  GSE81553_K562_ prefixes dropped" % (args.label_floor * 100))

    for unit, report in (("group", groups), ("feature", features)):
        print("\n=== %s ===" % unit)
        for (task, table_type), frame in report.items():
            counts = frame["lean"].value_counts()
            both = frame[(frame["exp"] > 0) & (frame["spl"] > 0)]
            print("%-14s %-11s n=%2d  rho=%.2f  median ×%.2f  splicing %2d / even %2d / expression %2d"
                  % (task, table_type, len(frame),
                     both["exp"].corr(both["spl"], method="spearman"), both["ratio"].median(),
                     counts.get("splicing", 0), counts.get("even", 0),
                     counts.get("expression", 0)))
        frame = report[("regression", "intragenic")]
        print("  regression/intragenic, the extremes:")
        ordered = frame.dropna(subset=["ratio"]).sort_values("ratio")
        for key in list(ordered.index[:5]) + list(ordered.index[-5:]):
            row = frame.loc[key]
            print("    %-38s exp %.4f  spl %.4f  ×%.2f  %s"
                  % (key, row["exp"], row["spl"], row["ratio"], row["lean"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
