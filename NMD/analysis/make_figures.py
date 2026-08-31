#!/usr/bin/env python3
"""
Figures for analysis/README.ru.md -- how the NMD prediction relates to the measured CHX response.

Descriptive only: distributions, binned medians, composition. No statistical tests.

    module load matplotlib/3.7.0-gfbf-2022b
    python3 analysis/make_figures.py

in   output/nmd_table_{5ss_noP,CPA_noP}.tsv     -- the prediction + the CHX readout
     work/transcripts.tsv                       -- tx_len, for relative position
     logFC/CHX_{5ss,CPA}.txt                    -- ChromHMM and host-gene expression
out  analysis/pics/*.png
"""

import csv
import math
import os
import random
import sys
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PICS = os.path.join(HERE, "pics")

# Palette. The default blue/orange pair was dropped on request, so this uses slots 6 (green) and
# 7 (violet) of the documented dataviz palette -- adjacent in its documented order, therefore
# validated as a pair -- plus a neutral for the deliberately recessive baseline series. Only two
# categorical hues are ever on screen at once, so the all-pairs cap does not bite. The sequential
# ramp is one hue (violet) stepped light -> dark, as an ordered variable requires.
# Light mode only; these are PNGs embedded in a local document.
BLUE, ORANGE, AQUA = "#008300", "#4a3aa7", "#7a6f3f"
GREY, INK, INK2, GRID = "#8a8a85", "#0b0b0b", "#52514e", "#e3e3e0"
SEQ = ["#ded9f5", "#b6abe6", "#8b7bd2", "#6455b9", "#453398"]   # one hue, light -> dark

LIBS = [("5ss_noP", "5ss"), ("CPA_noP", "CPA")]
# Panel titles. This used to map the libraries to the manuscript's wording, "Regular" and "PCPA",
# and that was rejected: the figures are read next to the tables, the logs and the file names, and
# every one of those says 5ss_noP and CPA_noP, so a second set of names for the same two things
# only makes the reader translate. Kept as a mapping rather than deleted because
# export_panel_data.py imports it and fig1G is drawn from what it writes.
DISPLAY = {"5ss_noP": "5ss_noP", "CPA_noP": "CPA_noP"}
DOTS_UP_TO = 500        # boxplots smaller than this also show every barcode as a dot

# Arial, on request, for every figure this script draws.
#
# It has to be registered from a file. This cluster has no Arial and no Helvetica installed -
# fc-list returns 38 entries, all of them DejaVu Sans, Source Code Pro or Cantarell - so asking
# matplotlib for "Arial" without this would fall through to DejaVu Sans silently: the figure would
# look unchanged and be set in the wrong face, which is the failure mode worth guarding against.
#
# One glyph the file does NOT carry is U+2082, the subscript two of "log2". Every label that needs
# it is written as mathtext instead and composed by matplotlib out of the regular digits, with the
# mathtext faces pointed at Arial below so the formula matches the prose around it.
ARIAL = os.path.normpath(os.path.join(ROOT, os.pardir, "context", "other", "Arial.ttf"))
if os.path.exists(ARIAL):
    fm.fontManager.addfont(ARIAL)
    FONT = fm.FontProperties(fname=ARIAL).get_name()
    resolved = fm.findfont(fm.FontProperties(family=FONT), fallback_to_default=False)
    if os.path.realpath(resolved) != os.path.realpath(ARIAL):
        raise RuntimeError("asked for %s and matplotlib resolved %s" % (FONT, resolved))
else:
    FONT = "DejaVu Sans"
    print("  [!] no Arial at %s - figures will be set in %s" % (ARIAL, FONT))

plt.rcParams.update({
    "font.family": FONT, "font.sans-serif": [FONT],
    # mathtext has its own font set, and left alone it would set every $...$ in DejaVu while the
    # labels around it are Arial.
    "mathtext.fontset": "custom", "mathtext.rm": FONT, "mathtext.it": "%s:italic" % FONT,
    "mathtext.bf": "%s:bold" % FONT, "mathtext.default": "regular",
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.edgecolor": INK2, "axes.linewidth": 0.8,
    "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "grid.color": GRID, "grid.linewidth": 0.6,
    "legend.frameon": False, "figure.dpi": 150, "savefig.bbox": "tight",
    # Vector output (save(..., vector=True)): keep the labels as real text rather than outlines --
    # TrueType embedded in the PDF, plain <text> in the SVG -- so a figure that goes to a journal
    # can still have its type edited and searched. No effect on the PNGs.
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})


# ----------------------------------------------------------------------------- journal style
# The four panels that go into the manuscript (fig01c, fig08, fig11_relpos, and fig1G, which is
# drawn by fig04_nature.R) follow the conventions Nature-family journals ask for: one sans face at
# one size, sentence case, nothing bold, no grid, no frame, no background tint, hairlines a printer
# can still hold, and the sample sizes on the figure rather than only in the caption.
#
# One thing cannot be matched here. fig1G is drawn in R because this cluster has neither Helvetica
# nor Arial, and matplotlib can only embed a font it can find, so the vector output of these three
# is set in DejaVu Sans. The geometry and the line weights match fig1G; the typeface does not. The
# route to real Helvetica is the one fig1G already takes - dump the panel to TSV, draw it in R.
MM = 1 / 25.4                       # matplotlib works in inches, journals in millimetres
COL2 = 183 * MM                     # the double-column measure, so 6 pt is 6 pt on the page
HAIR = 0.5                          # pt; 0.25 pt is the print floor but vanishes on screen at 150 dpi
FS = 7                              # one size for every label
FS_SMALL = 6                        # sample sizes and footnotes

# Four ordered categories of where a stop codon comes from, so a single hue stepped light -> dark
# rather than four unrelated colours: the order 5' -> 3' along the transcript is the information.
# Same violet as the "PTC created" group of fig1G, so the two panels read as one figure.
SEQ_STROKE = ["#8b7bd2", "#6455b9", "#4a3aa7", "#2e2170"]


def journal(ax, xgrid=False, ygrid=False):
    """Strip an axis down to a hairline x rule with outward ticks, and nothing else."""
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_linewidth(HAIR)
    ax.spines["bottom"].set_color(INK)
    ax.tick_params(direction="out", length=2.2, width=HAIR, colors=INK, labelsize=FS, pad=2)
    ax.tick_params(axis="y", length=0)
    for axis, on in (("x", xgrid), ("y", ygrid)):
        # grid(False, colour=...) is read by matplotlib as "enable it, with these properties",
        # so the properties can only be passed on the branch that actually wants a grid.
        if on:
            ax.grid(True, axis=axis, color=GRID, linewidth=HAIR, zorder=0)
        else:
            ax.grid(False, axis=axis)
    ax.set_axisbelow(True)


def x_rule(ax, ticks):
    """Clip the x rule to the tick range.

    A journal axis is a scale, not a border: it has to stop where the numbers stop. matplotlib
    draws the spine across the whole panel, so a panel with a reserved strip on the right - for
    stars, or for end labels - grows a rule into empty space. ggplot has the same habit; the R
    panel fixes it the same way."""
    ax.spines["bottom"].set_bounds(min(ticks), max(ticks))


def stars(p, n_min=None, min_n=5):
    """The conventional significance ladder, which is NOT the number of decimals of p.

    ns / * / ** / *** / **** at 0.05 / 0.01 / 0.001 / 0.0001, and it stops at four: p = 5.8e-42
    and p = 9.0e-57 are both '****'. The exact numbers therefore have to travel in the caption,
    which is what the console output of each figure is for.

    A borderline p is the weakness of the notation: p = 0.05003 prints 'ns' while 0.0499 would
    print '*', and nothing about the data distinguishes them."""
    if n_min is not None and n_min < min_n:
        return ""
    if p is None or p != p:
        return ""
    for cut, mark in ((1e-4, "****"), (1e-3, "***"), (1e-2, "**"), (5e-2, "*")):
        if p < cut:
            return mark
    return "ns"


def holm(pvals):
    """Holm-Bonferroni adjusted p values, returned in the input order.

    Same guarantee on the family-wise error rate as plain Bonferroni and uniformly more powerful,
    so there is no reason to use Bonferroni instead. Sort ascending, multiply the k-th by (m - k),
    then enforce monotonicity so an adjusted value can never fall below an earlier one.

    The family is the set of comparisons against the shared reference WITHIN ONE LIBRARY. The two
    libraries are different constructs measured on different barcodes - separate experiments, not
    two looks at one question - so correcting across them would make each experiment penalise the
    other."""
    m = len(pvals)
    out = [1.0] * m
    running = 0.0
    for k, i in enumerate(sorted(range(m), key=lambda j: pvals[j])):
        running = max(running, (m - k) * pvals[i])
        out[i] = min(1.0, running)
    return out


def mannwhitney(a, b):
    """Two-sided Mann-Whitney U, the two-sample form of the Wilcoxon rank-sum test.

    Written out rather than imported: this module deliberately depends on nothing but matplotlib,
    and the normal approximation is what scipy would use at these sample sizes anyway. Ties get
    average ranks and the variance is tie-corrected."""
    if len(a) < 3 or len(b) < 3:
        return None
    pool = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    ranks = [0.0] * len(pool)
    ties = []
    i = 0
    while i < len(pool):
        j = i
        while j + 1 < len(pool) and pool[j + 1][0] == pool[i][0]:
            j += 1
        r = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[k] = r
        ties.append(j - i + 1)
        i = j + 1
    n1, n2 = len(a), len(b)
    r1 = sum(rk for rk, (_, g) in zip(ranks, pool) if g == 0)
    u1 = r1 - n1 * (n1 + 1) / 2.0
    n = n1 + n2
    mu = n1 * n2 / 2.0
    tie_term = sum(t ** 3 - t for t in ties)
    var = n1 * n2 / 12.0 * ((n + 1) - tie_term / (n * (n - 1.0)))
    if var <= 0:
        return None
    z = (abs(u1 - mu) - 0.5) / math.sqrt(var)
    return math.erfc(z / math.sqrt(2.0))


# ----------------------------------------------------------------------------- data


def med(v):
    v = sorted(v)
    n = len(v)
    if not n:
        return float("nan")
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2


def quantile(v, q):
    v = sorted(v)
    if not v:
        return float("nan")
    i = min(len(v) - 1, max(0, int(round(q * (len(v) - 1)))))
    return v[i]


def load():
    txlen, cds_start, span = {}, {}, {}
    with open(os.path.join(ROOT, "work", "transcripts.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            txlen[r["tx_id"]] = int(r["tx_len"])
            span[r["tx_id"]] = (int(r["tx_start"]), int(r["tx_end"]), r["strand"])
            if r["cds_tx_start"]:
                cds_start[r["tx_id"]] = int(r["cds_tx_start"])

    data = {}
    for lib, chx in LIBS:
        extra = {}
        with open(os.path.join(ROOT, "logFC", "CHX_%s.txt" % chx)) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                extra[r["Barcode"]] = (r["ChromHMM"], r["MANE_all_exp.norm.log10"])
        rows = []
        with open(os.path.join(ROOT, "output", "nmd_table_%s.tsv" % lib)) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                r["lfc"] = float(r["CHXlog2FoldChange"]) if r["CHXlog2FoldChange"] else None
                r["base"] = float(r["CHXbaseMean"]) if r["CHXbaseMean"] else None
                ch, he = extra.get(r["barcode"], ("", ""))
                r["chromhmm"] = ch if ch not in ("NA", "") else ""
                r["host_exp"] = float(he) if he not in ("NA", "") else None
                r["tx_len"] = txlen.get(r["tx_id"])
                r["relpos"] = (int(r["tx_pos"]) / r["tx_len"]
                               if r["tx_pos"] and r["tx_len"] else None)
                sp = span.get(r["tx_id"])
                if sp and r["in_gene"] == "1":
                    a, b, st = sp
                    f = (int(r["pos"]) - a) / (b - a)
                    r["relgene"] = f if st == "+" else 1 - f
                else:
                    r["relgene"] = None
                cs = cds_start.get(r["tx_id"])
                r["atg_dist"] = (int(r["ptc_tx_pos"]) - cs
                                 if r["ptc_tx_pos"] and cs is not None else None)
                r["i5"] = int(r["feat_idx"]) if r["feat_idx"] and r["loc"] == "intronic" else None
                r["i3"] = (int(r["n_exons"]) - int(r["feat_idx"])
                           if r["feat_idx"] and r["loc"] == "intronic" else None)
                r["has_ptc"] = (
                    (r["verdict_reason"] == "cds_intron" and r["ptc_source"] not in ("", "none"))
                    or r["verdict_reason"] in ("internal_exon", "last_exon"))
                # A PTC the last-exon rule and the 50-nt rule both let through, i.e. one that is
                # actually expected to trigger NMD. has_ptc is deliberately wider: it counts a
                # stop in the last exon too, where there is no downstream junction and NMD is not
                # expected. Panels that ask "does making a PTC change the CHX response" want this
                # narrower flag; panels that ask "where do stops appear" want has_ptc.
                #
                # Never use `not nmd_ptc` to build a no-PTC control group: a barcode with a PTC
                # that fails the 50-nt rule would land in it. `not has_ptc` is the control.
                r["nmd_ptc"] = r["has_ptc"] and r["nmd_expected"] == "1"
                rows.append(r)
        data[lib] = rows
    return data


def expressed(rows):
    return [r for r in rows if r["lfc"] is not None]


def save(fig, name, vector=False):
    """PNG always; `vector` additionally writes SVG and PDF of the same figure.

    Figures that go into the paper are asked for as vector art, so those pass vector=True and get
    a 300 dpi PNG on top for quick viewing. Text is kept as text in the SVG/PDF (no path
    conversion) so the labels stay editable downstream."""
    fig.savefig(os.path.join(PICS, name), facecolor=fig.get_facecolor(),
                dpi=300 if vector else None)
    print("  pics/%s" % name)
    if vector:
        stem = os.path.splitext(name)[0]
        for ext in ("svg", "pdf"):
            fig.savefig(os.path.join(PICS, "%s.%s" % (stem, ext)),
                        facecolor=fig.get_facecolor())
            print("  pics/%s.%s" % (stem, ext))
    plt.close(fig)


def tidy(ax, xgrid=False, ygrid=True):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    ax.set_axisbelow(True)
    if ygrid:
        ax.yaxis.grid(True)
    if xgrid:
        ax.xaxis.grid(True)


def box_groups(ax, groups, colors, xlim=(-3.5, 6.5), ylabels=True, min_n=5,
               xlabel="log2 fold change (CHX / EtOH)", fontsize=None):
    """Horizontal boxes over a jittered strip of the individual barcodes.

    For a group of fewer than DOTS_UP_TO barcodes every barcode is also drawn as one dot behind the
    box, so the reader sees the actual sample rather than five summary numbers; the box is kept
    translucent for the dots to show through. Bigger groups get the box alone. Whiskers
    are the 5th/95th percentiles -- with a handful of barcodes at log2FC 12 the full range squashes
    every box into a sliver.

    `xlim=None` leaves the limits to the caller and only returns the (min, max) actually drawn, so
    a caller with two panels can pick one shared range that holds every whisker end -- a fixed
    range clips them, and a clipped whisker reads as a shorter tail than the data has.
    `fontsize` overrides the size of the tick and axis labels of this axis."""
    rnd = random.Random(0)                       # fixed seed: the jitter must not move between runs
    for i, (g, c) in enumerate(zip(groups, colors), start=1):
        v = g[1]
        # Individual barcodes are drawn only for groups small enough for a dot to mean something.
        # Past a few hundred the strip is a smear that hides the box instead of informing about it.
        if not v or len(v) >= DOTS_UP_TO:
            continue
        ax.scatter(v, [i + rnd.uniform(-0.3, 0.3) for _ in v], s=4, alpha=0.30, color=c,
                   linewidths=0, rasterized=True, zorder=1)
    vals = [g[1] for g in groups]
    # A group too small for a box keeps its row (so the two panels stay aligned) but is shown as
    # a bare median dot -- a box drawn from 2 points would read as a real distribution.
    big = [v if len(v) >= min_n else [] for v in vals]
    bp = ax.boxplot([v if v else [0] for v in big], vert=False, showfliers=False, whis=(5, 95),
                    widths=0.6, patch_artist=True, zorder=3,
                    medianprops=dict(color=INK, linewidth=1.6),
                    whiskerprops=dict(color=INK2, linewidth=0.9),
                    capprops=dict(color=INK2, linewidth=0.9),
                    boxprops=dict(linewidth=0.8, edgecolor=INK2))
    for i, (patch, c) in enumerate(zip(bp["boxes"], colors)):
        if big[i]:
            patch.set_facecolor(c)
            patch.set_alpha(0.38)
        else:
            patch.set_visible(False)
            bp["medians"][i].set_visible(False)
            for k in ("whiskers", "caps"):
                bp[k][2 * i].set_visible(False)
                bp[k][2 * i + 1].set_visible(False)
            ax.plot([med(vals[i])], [i + 1], "o", markersize=8, color=c,
                    markeredgecolor="#fcfcfb", markeredgewidth=1.4, zorder=4)
    ax.axvline(0, color=INK2, linewidth=0.9, linestyle=":")
    if ylabels:
        ax.set_yticklabels(["%s  (n=%d)" % (g[0], len(g[1])) for g in groups],
                           fontsize=fontsize)
    else:
        ax.set_yticklabels(["n=%d" % len(g[1]) for g in groups],
                           fontsize=fontsize or 7.5, color=INK2)
    if xlim is not None:
        ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    if fontsize:
        ax.tick_params(axis="both", labelsize=fontsize)
    tidy(ax, xgrid=True, ygrid=False)

    # Extent of what is drawn as a summary: whisker ends of the real boxes, plus the bare median
    # dot of the groups too small for one. The caller uses it to size the axis.
    drawn = [x for w in bp["whiskers"] if w.get_visible() for x in w.get_xdata()]
    drawn += [med(vals[i]) for i, v in enumerate(big) if not v and vals[i]]
    return (min(drawn), max(drawn)) if drawn else (0.0, 0.0)


def binned(ax, xs, ys, nbins, color, label=None, band=True):
    """Median of ys within equal-count bins of xs, with an IQR band."""
    pts = sorted(zip(xs, ys))
    if len(pts) < nbins * 5:
        nbins = max(2, len(pts) // 20)
    step = len(pts) / nbins
    cx, cm, lo, hi, ns = [], [], [], [], []
    for b in range(nbins):
        chunk = pts[int(b * step):int((b + 1) * step)]
        if len(chunk) < 5:
            continue
        cx.append(med([p[0] for p in chunk]))
        v = [p[1] for p in chunk]
        cm.append(med(v))
        lo.append(quantile(v, 0.25))
        hi.append(quantile(v, 0.75))
        ns.append(len(chunk))
    if band:
        ax.fill_between(cx, lo, hi, color=color, alpha=0.15, linewidth=0)
    ax.plot(cx, cm, color=color, linewidth=1.4, marker="o", markersize=3,
            markeredgecolor="#fcfcfb", markeredgewidth=0.7, label=label)
    # The bin centres and medians come back so a caller can label the end of the line instead of
    # putting a legend box on top of the data. No caller used the old return value.
    return cx, cm, ns


# ----------------------------------------------------------------------------- figures


def fig01_composition(data):
    fig, axes = plt.subplots(2, 1, figsize=(8, 3.4), sharex=True)
    order = [("not expressed", "1", GREY),
             ("outside a gene / antisense", "2", "#c9c9c4"),
             ("no NMD expected", "3", BLUE),
             ("NMD expected", "4", ORANGE)]
    for ax, (lib, _) in zip(axes, LIBS):
        c = Counter(r["category"] for r in data[lib])
        left = 0
        tot = sum(c.values())
        for label, key, col in order:
            w = c[key]
            ax.barh([0], [w], left=left, color=col, alpha=0.85, height=0.55,
                    edgecolor="#fcfcfb", linewidth=2, label=label if ax is axes[0] else None)
            if w / tot > 0.03:
                ax.text(left + w / 2, 0, "%d" % w, ha="center", va="center",
                        fontsize=8, color=INK)
            left += w
        ax.set_yticks([0])
        ax.set_yticklabels([lib])
        ax.set_xlim(0, tot)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
    axes[1].set_xlabel("number of barcodes")
    axes[0].legend(ncol=4, loc="lower left", bbox_to_anchor=(0, 1.05), fontsize=8)
    save(fig, "fig01_composition.png")


def fig01b_loc_counts(data):
    """How many integrations landed in introns and in exons, and how many of those are expressed."""
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 2.4), sharex=True)
    for ax, (lib, _) in zip(axes, LIBS):
        rows = [r for r in data[lib] if r["in_gene"] == "1"]
        for i, loc in enumerate(("exonic", "intronic")):
            s = [r for r in rows if r["loc"] == loc]
            e = sum(1 for r in s if r["lfc"] is not None)
            ax.barh([i], [e], color=ORANGE, alpha=0.85, height=0.5, edgecolor="#fcfcfb",
                    linewidth=1.5, label="expressed" if i == 0 else None)
            ax.barh([i], [len(s) - e], left=[e], color=GREY, alpha=0.6, height=0.5,
                    edgecolor="#fcfcfb", linewidth=1.5,
                    label="not expressed" if i == 0 else None)
            ax.text(len(s) + 250, i, "%d of %d expressed  (%.0f%%)"
                    % (e, len(s), 100 * e / len(s)), va="center", fontsize=7.5, color=INK2)
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["exonic", "intronic"] if ax is axes[0] else ["", ""])
        ax.set_xlim(0, 21000)
        ax.set_title(lib)
        ax.set_xlabel("barcodes inside a gene, same strand")
        tidy(ax, xgrid=True, ygrid=False)
    axes[0].legend(fontsize=8, ncol=2, loc="lower left", bbox_to_anchor=(0, 1.12))
    save(fig, "fig01b_loc_counts.png")


def fig01c_position_landscape(data):
    """Where integrations land along a gene, every gene rescaled to 0..1 so genes of different
    length and exon count pool together. 20 equal bins.

    Top row: along the gene as it sits in the genome. Bottom row: along the mature mRNA -- for an
    intronic hit that is the junction the cassette splices into, since the intron itself is not in
    the mRNA. The dotted line at 1.0 is what a uniform landscape would give.

    No legend box anywhere. Each series is named where it is drawn, in its own colour, next to the
    count for THAT panel -- which is the point: the two libraries do not share those counts, so one
    legend for the figure would be wrong, and four legend boxes were just four copies of the same
    two words."""
    fig, axes = plt.subplots(2, 2, figsize=(COL2, 88 * MM), sharex=True, sharey="row")
    fig.subplots_adjust(left=0.085, right=0.895, top=0.935, bottom=0.165, hspace=0.30, wspace=0.11)
    bins = [i / 20 for i in range(21)]
    rows_spec = [("relgene", "along the gene\n(genomic span)"),
                 ("relpos", "along the\nmature mRNA")]
    series = (("all integrations", lambda r: True, GREY),
              ("expressed only", lambda r: r["lfc"] is not None, ORANGE))
    print("  fig01c: integrations per panel")
    for row, (key, strip) in enumerate(rows_spec):
        for col, (lib, _) in enumerate(LIBS):
            ax = axes[row][col]
            base = [r for r in data[lib] if r["in_gene"] == "1" and r.get(key) is not None]
            for k, (label, sel, colour) in enumerate(series):
                v = [r[key] for r in base if sel(r)]
                ax.hist(v, bins=bins, density=True, histtype="step", linewidth=1.1,
                        color=colour, zorder=2)
                # Named once per column, in the top row only: an in-gene barcode always has both
                # relgene and relpos, so the bottom row is the same barcodes and the same counts,
                # and repeating them would be four labels carrying two facts.
                if row == 0:
                    ax.text(0.975, 0.95 - 0.105 * k,
                            "%s   n = %s" % (label, format(len(v), ",")),
                            transform=ax.transAxes, ha="right", va="top",
                            fontsize=FS_SMALL, color=colour)
                print("    %-8s %-22s %-18s n = %s" % (lib, key, label, format(len(v), ",")))
            ax.axhline(1.0, color=INK2, linewidth=HAIR, linestyle=(0, (1, 2.5)), zorder=1)
            journal(ax)
            ax.set_xlim(0, 1)
            ticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
            ax.set_xticks(ticks)
            ax.set_xticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"])
            x_rule(ax, ticks)
            if row == 0:
                ax.set_title(DISPLAY[lib], fontsize=FS, color=INK, pad=4)
            if col == 0:
                # Names the dotted reference once per row, at the left where both rows are empty.
                ax.text(0.015, 1.0, "uniform", transform=ax.get_yaxis_transform(),
                        ha="left", va="bottom", fontsize=FS_SMALL, color=INK2)
            if col == 1:
                # The row identity as a strip on the right, the way a facet label reads.
                ax.text(1.035, 0.5, strip, transform=ax.transAxes, rotation=270,
                        ha="left", va="center", fontsize=FS_SMALL, color=INK2,
                        linespacing=1.35)
    # Two separate texts, and they have to be far enough apart in FIGURE fraction: savefig runs
    # with bbox="tight", which crops to the outermost ink, so two labels 0.02 apart on a 1,000 px
    # canvas overlap in the output even though the margins look generous in the code.
    fig.text(0.49, 0.058, "relative position along the gene: 0 = 5\u2032 end, 1 = 3\u2032 end",
             ha="center", fontsize=FS, color=INK)
    fig.text(0.012, 0.55, "density of integrations", rotation=90, va="center",
             fontsize=FS, color=INK)
    save(fig, "fig01c_position_landscape.png", vector=True)


def fig02_ecdf(data):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), sharey=True)
    series = [("outside a gene / antisense", "2", GREY), ("no NMD expected", "3", BLUE),
              ("NMD expected", "4", ORANGE)]
    for ax, (lib, _) in zip(axes, LIBS):
        for label, key, col in series:
            v = sorted(r["lfc"] for r in expressed(data[lib]) if r["category"] == key)
            ax.plot(v, [i / len(v) for i in range(1, len(v) + 1)], color=col, linewidth=2,
                    label="%s (n=%d)" % (label, len(v)))
        ax.axvline(0, color=INK2, linewidth=0.9, linestyle=":")
        ax.set_xlim(-3, 4)
        ax.set_title(lib)
        ax.set_xlabel("log2 fold change (CHX / EtOH)")
        tidy(ax)
        ax.legend(fontsize=7.5, loc="upper left")
    axes[0].set_ylabel("fraction of barcodes <= x")
    save(fig, "fig02_ecdf_category.png")


def fig03_verdict(data):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    # The computed branch is split by what the rule actually says about each row, so nothing is
    # averaged together that the rule treats differently. "intronic" marks the computed groups.
    def cdsi(r, lo, hi):
        return (r["verdict_reason"] == "cds_intron" and r["dist_to_last_junction"]
                and lo < int(r["dist_to_last_junction"]) <= hi)
    spec = [("outside a gene / antisense", lambda r: r["category"] == "2", GREY),
            ("UTR - no PTC created", lambda r: r["verdict_reason"] == "utr", BLUE),
            ("exonic, last coding exon\nrule: escape",
             lambda r: r["verdict_reason"] == "last_exon", BLUE),
            ("exonic, earlier coding exon\nrule: NMD",
             lambda r: r["verdict_reason"] == "internal_exon", ORANGE),
            ("intronic, PTC >50 nt before\nthe last junction   rule: NMD",
             lambda r: cdsi(r, 50, 10 ** 9), ORANGE),
            ("intronic, PTC <=50 nt before\nthe last junction   rule: escape",
             lambda r: cdsi(r, 0, 50), BLUE),
            ("intronic, PTC in the last\ncoding exon   rule: escape",
             lambda r: cdsi(r, -10 ** 9, 0), BLUE),
            ("intronic, no PTC at all -- the\nframeshift made none",
             lambda r: r["verdict_reason"] == "no_ptc_frameshift", GREY)]
    for ax, (lib, _) in zip(axes, LIBS):
        rows = expressed(data[lib])
        groups, cols = [], []
        for label, sel, col in spec:
            v = [r["lfc"] for r in rows if sel(r)]
            if v:
                groups.append((label, v))
                cols.append(col)
        box_groups(ax, groups, cols, ylabels=(ax is axes[0]))
        ax.set_title(lib)
        ax.tick_params(axis="y", labelsize=7.5)
    save(fig, "fig03_box_verdict.png")


def fig04_ptc(data):
    """Panel 1G of the paper, so it follows the manuscript conventions rather than the defaults of
    this script: every label is 6 pt, the figure is
    laid out at its intended printed width so that 6 pt stays 6 pt on the page, and it is written
    as SVG and PDF next to the PNG.

    The shared x range is taken from the whiskers the two panels actually drew instead of the fixed
    default -- at the default the 95th-percentile whisker of "PTC created" ran off the right edge in
    both panels, which reads as a shorter tail than the data has."""
    fs = 6
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.0))
    spans = []
    for ax, (lib, _) in zip(axes, LIBS):
        rows = expressed(data[lib])
        groups = [("outside a gene (baseline)", [r["lfc"] for r in rows if r["category"] == "2"]),
                  ("no PTC created", [r["lfc"] for r in rows
                                        if r["in_gene"] == "1" and not r["has_ptc"]]),
                  ("PTC created", [r["lfc"] for r in rows if r["has_ptc"]])]
        spans.append(box_groups(ax, groups, [GREY, BLUE, ORANGE], xlim=None,
                                ylabels=(ax is axes[0]), fontsize=fs))
        ax.set_title(DISPLAY[lib], fontsize=fs)
    lo, hi = min(a for a, _ in spans), max(b for _, b in spans)
    pad = 0.04 * (hi - lo)
    for ax in axes:
        ax.set_xlim(lo - pad, hi + pad)
    save(fig, "fig04_ptc_yesno.png", vector=True)


def fig05_dist_hist(data):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.2), sharey=True)
    for ax, (lib, _) in zip(axes, LIBS):
        v = [int(r["dist_to_last_junction"]) for r in data[lib]
             if r["dist_to_last_junction"] and int(r["dist_to_last_junction"]) > 0]
        ax.hist(v, bins=[10 ** (i / 8) for i in range(0, 45)], color=BLUE, alpha=0.75,
                edgecolor="#fcfcfb", linewidth=0.4)
        ax.axvline(50, color=ORANGE, linewidth=2)
        ax.text(55, ax.get_ylim()[1] * 0.9, "50 nt threshold", color=ORANGE, fontsize=8)
        ax.set_xscale("log")
        ax.set_title("%s (n=%d with a positive distance)" % (lib, len(v)))
        ax.set_xlabel("distance from the PTC to the last junction, nt")
        tidy(ax)
    axes[0].set_ylabel("number of barcodes")
    save(fig, "fig05_dist_hist.png")


def fig07_matched(data):
    """Everything that looked like a violation of the rule, before and after matching on position.

    For each comparison, group A is restricted to its own 10-90 percentile band of relative
    position and group B is restricted to the same band, so the two are read at the same place
    along the transcript. Two of the three differences vanish; the one that matters does not.
    """
    tests = [("50 nt rule:\nescape - NMD",
              lambda r: r["verdict_reason"] == "cds_intron" and r["dist_num"] is not None
              and r["dist_num"] <= 50,
              lambda r: r["verdict_reason"] == "cds_intron" and r["dist_num"] is not None
              and r["dist_num"] > 50),
             ("3'UTR exonic -\nno-PTC control",
              lambda r: r["region"] == "3utr" and r["loc"] == "exonic",
              lambda r: r["verdict_reason"] == "utr"),
             ("has PTC -\nno PTC",
              lambda r: r["verdict_reason"] == "cds_intron" and r["dist_num"] is not None,
              lambda r: r["verdict_reason"] == "utr")]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.4), sharey=True)
    for ax, (lib, _) in zip(axes, LIBS):
        rows = [r for r in expressed(data[lib]) if r["in_gene"] == "1" and r["relpos"] is not None]
        for r in rows:
            r["dist_num"] = (int(r["dist_to_last_junction"])
                             if r["dist_to_last_junction"] else None)
        w = 0.38
        for i, (label, A, B) in enumerate(tests):
            a = [r for r in rows if A(r)]
            b = [r for r in rows if B(r)]
            raw = med([r["lfc"] for r in a]) - med([r["lfc"] for r in b])
            lo, hi = quantile([r["relpos"] for r in a], .1), quantile([r["relpos"] for r in a], .9)
            am = [r["lfc"] for r in a if lo <= r["relpos"] <= hi]
            bm = [r["lfc"] for r in b if lo <= r["relpos"] <= hi]
            matched = med(am) - med(bm)
            for j, (val, col, lab) in enumerate(((raw, GREY, "raw"),
                                                 (matched, ORANGE, "matched on position"))):
                ax.bar([i + (j - 0.5) * w], [val], width=w * 0.9, color=col, alpha=0.85,
                       edgecolor="#fcfcfb", linewidth=1.2, label=lab if i == 0 else None)
                ax.text(i + (j - 0.5) * w, val + (0.04 if val >= 0 else -0.11),
                        "%+.2f" % val, ha="center", fontsize=7.5, color=INK)
        ax.axhline(0, color=INK2, linewidth=0.9)
        ax.set_xticks(range(len(tests)))
        ax.set_xticklabels([t[0] for t in tests], fontsize=7.5)
        ax.set_ylim(-0.35, 1.45)
        ax.set_title(lib)
        tidy(ax)
    axes[0].set_ylabel("difference in median log2FC")
    axes[0].legend(fontsize=8, ncol=2, loc="lower left", bbox_to_anchor=(0, 1.10))
    save(fig, "fig07_matched.png")


def fig08_ptc_source(data):
    """Where the premature stop actually comes from, against the CHX response.

    Restricted to PTCs expected to trigger NMD - the last-exon rule and the 50-nt rule both
    applied. Note the selection also gained a PTC filter it never had: the rows used to be picked
    on ptc_source alone, so a barcode with a source recorded but no PTC by verdict was counted.

    Box: 25th to 75th percentile, black line the median. WHISKERS ARE THE 5th AND 95th
    PERCENTILES, not 1.5 x IQR - a handful of barcodes past log2FC 10 would otherwise squash every
    box into a sliver. That is stated on the figure, because a reader assumes Tukey.

    Both panels share one x range, computed from the whisker ends of all eight boxes, so no tail is
    clipped and the two libraries are directly comparable. The old version had xlim hardcoded to
    (-3.5, 6.5) while the data ran to 12, which cut both ends of the CPA panel.

    The four sources are an ordered series along the transcript, from a stop inside the insert to
    reading through into the host, so they take one hue stepped light to dark rather than four
    unrelated colours.

    Each source is compared with the first one - a stop inside the insert, the canonical case and
    by far the most common - with the same two-sided rank-sum test and the same bracket-and-star
    notation fig1G uses, so the two panels can be read against each other. Expect mostly ns: the
    part of the transcript the stop comes from turns out to matter little, which is itself the
    result. Exact p values go to the console, for the caption."""
    # host_insert_seam is left out on purpose. It is a real category - the stop codon sits astride
    # the joint, its first bases the host's last and its final base the insert's first, most often
    # host TA + insert G = TAG - but it holds 4 barcodes in each library, so its row was a median
    # dot with no distribution behind it. Its absence should be stated in the caption rather than
    # left to look like the category does not exist. insert_3p, host_upstream and the other _seam
    # variants have zero NMD-competent barcodes and never had a row.
    names = [("stop inside the insert", "insert_5p"),
             ("stop from the barcode", "barcode"),
             ("read through, stop in the host", "host_downstream")]
    MIN_BOX = 5
    vals = {}
    for lib, _ in LIBS:
        rows = [r for r in expressed(data[lib]) if r["nmd_ptc"]]
        vals[lib] = {key: [r["lfc"] for r in rows if r["ptc_source"] == key]
                     for _, key in names}

    # One range for both panels, from what is actually drawn: whisker ends where there is a box,
    # the bare median where a group is too small for one.
    ends = []
    for lib, _ in LIBS:
        for _, key in names:
            v = vals[lib][key]
            if len(v) >= MIN_BOX:
                ends += [quantile(v, 0.05), quantile(v, 0.95)]
            elif v:
                ends.append(med(v))
    lo, hi = min(ends), max(ends)
    span = hi - lo
    # The rule has to reach past the longest whisker, or the end of that whisker cannot be read
    # off the scale. So the tick grid is chosen first, from the data, and the panel is then made
    # wide enough to contain it - not the other way round, which is how the rule ended up stopping
    # at 5 with whiskers running to 9.
    lo_tick = 5 * math.floor(lo / 5)
    hi_tick = 5 * math.ceil(hi / 5)
    ticks = list(range(lo_tick, hi_tick + 1, 5))
    n_x = lo - 0.26 * span                           # left strip holds the sample sizes
    # Right strip: one spine per comparison, the taller bracket further out, exactly as in fig1G.
    br_x = [hi + 0.05 * span, hi + 0.11 * span]
    br_tick = 0.02 * span                            # arms reaching back towards the boxes
    star_dx = 0.028 * span
    xlim = (min(lo - 0.28 * span, lo_tick - 0.01 * span),
            max(hi_tick + 0.01 * span, br_x[-1] + star_dx + 0.05 * span))

    fig, axes = plt.subplots(1, 2, figsize=(COL2, 62 * MM), sharex=True, sharey=True)
    # wspace has to clear two tick labels meeting head on: at 0.05 the right end of the left panel
    # and the left end of the right one printed as "10-10".
    fig.subplots_adjust(left=0.235, right=0.99, top=0.88, bottom=0.20, wspace=0.12)
    print("  fig08: source of the PTC, NMD-expected only")
    for ax, (lib, _) in zip(axes, LIBS):
        ref = vals[lib][names[0][1]]
        # Every comparison of this library at once, so the correction can see the whole family
        # before any of it is drawn.
        raw = [mannwhitney(ref, vals[lib][key]) for _, key in names[1:]]
        adj = holm([x if x is not None else 1.0 for x in raw])
        adj = [a if r is not None else None for a, r in zip(adj, raw)]
        ys, boxed = [], []
        for i, (label, key) in enumerate(names):
            y = len(names) - i                        # first name on top
            v = vals[lib][key]
            ys.append((y, label, key, v))
            boxed.append(v if len(v) >= MIN_BOX else [])
        bp = ax.boxplot([v if v else [0] for _, _, _, v in ys],
                        positions=[y for y, _, _, _ in ys], vert=False, showfliers=False,
                        whis=(5, 95), widths=0.52, patch_artist=True, zorder=3,
                        medianprops=dict(color=INK, linewidth=1.1),
                        whiskerprops=dict(color=INK2, linewidth=HAIR),
                        capprops=dict(color=INK2, linewidth=HAIR),
                        boxprops=dict(linewidth=HAIR))
        for i, ((y, label, key, v), patch) in enumerate(zip(ys, bp["boxes"])):
            stroke = SEQ_STROKE[i]
            if len(v) >= MIN_BOX:
                patch.set_facecolor(stroke)
                patch.set_alpha(0.30)
                patch.set_edgecolor(stroke)
            else:
                patch.set_visible(False)
                bp["medians"][i].set_visible(False)
                for k in ("whiskers", "caps"):
                    bp[k][2 * i].set_visible(False)
                    bp[k][2 * i + 1].set_visible(False)
                if v:
                    ax.plot([med(v)], [y], "o", markersize=4.5, color=stroke,
                            markeredgecolor="#fcfcfb", markeredgewidth=0.8, zorder=4)
            p = raw[i - 1] if i > 0 else None
            p_adj = adj[i - 1] if i > 0 else None
            mark = stars(p_adj, n_min=min(len(ref), len(v)), min_n=MIN_BOX) if i > 0 else ""
            if i > 0 and mark:
                # A vertical spine in the right strip, two short arms reaching back to the rows it
                # joins, and the mark rotated to read along it. "ns" turns the other way and one
                # point smaller: it is a word, so it has a reading direction, while a stack of
                # asterisks looks the same either way - the same rule fig04_nature.R follows.
                x = br_x[i - 1]
                y_ref = len(names)
                ax.plot([x, x], [y_ref, y], color=INK2, linewidth=HAIR, zorder=5)
                for yy in (y_ref, y):
                    ax.plot([x - br_tick, x], [yy, yy], color=INK2, linewidth=HAIR, zorder=5)
                word = mark == "ns"
                ax.text(x + star_dx, (y_ref + y) / 2, mark, rotation=-90 if word else 90,
                        ha="center", va="center",
                        fontsize=FS_SMALL - 1 if word else FS_SMALL, color=INK2)
            print("    %-8s %-30s n = %-6s median %+6.3f   %s"
                  % (lib, key, format(len(v), ","), med(v) if v else float("nan"),
                     ("p = %.3g, Holm-adjusted %.3g  %s" % (p, p_adj, mark))
                     if p is not None else "reference"))
        ax.axvline(0, color=INK2, linewidth=HAIR, linestyle=(0, (1, 2.5)), zorder=1)
        ax.set_title(DISPLAY[lib], fontsize=FS, color=INK, pad=4)
        journal(ax)
        ax.set_xlim(*xlim)
        ax.set_ylim(0.45, len(names) + 0.75)
        ax.set_yticks([])
        ax.set_xticks(ticks)
        x_rule(ax, ticks)
        # Each panel carries its own counts, inside itself, in the empty room left of the
        # whiskers. Putting the right panel's column outside its axes ran it into the left
        # panel's stars.
        for y, _, key, _ in ys:
            ax.text(n_x, y, "n = %s" % format(len(vals[lib][key]), ","), ha="left",
                    va="center", fontsize=FS_SMALL, color=INK2)

    # The category names once, in the margin: they are shared by both panels.
    for y, label, key, v in ys:
        axes[0].annotate(label, xy=(0, y), xycoords=("axes fraction", "data"),
                         xytext=(-7, 0), textcoords="offset points", ha="right", va="center",
                         fontsize=FS, color=INK)
    fig.text(0.62, 0.055, "$\\mathregular{log_2}$ fold change (CHX / EtOH)", ha="center",
             fontsize=FS, color=INK)
    save(fig, "fig08_ptc_source.png", vector=True)


def fig09_frame(data):
    fig, axes = plt.subplots(1, 2, figsize=(10, 2.8))
    for ax, (lib, _) in zip(axes, LIBS):
        rows = [r for r in expressed(data[lib]) if r["insert_frame"]]
        groups, cols = [], []
        for f, note in (("0", "stop in the insert at nt 54"),
                        ("1", "no stop in the insert"),
                        ("2", "stop in the insert at nt 23")):
            v = [r["lfc"] for r in rows if r["insert_frame"] == f]
            groups.append(("frame %s\n%s" % (f, note), v))
            cols.append(BLUE if f == "1" else ORANGE)
        box_groups(ax, groups, cols, ylabels=(ax is axes[0]))
        ax.set_title(lib)
        ax.tick_params(axis="y", labelsize=7.5)
    save(fig, "fig09_frame.png")


def fig10_intron_index(data):
    """Median response by which intron the insertion landed in -- counted from the 5' end and,
    separately, from the 3' end. Discrete cross-check of the positional effect."""
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.2), sharey=True)
    for ax, (lib, _) in zip(axes, LIBS):
        rows = [r for r in expressed(data[lib])
                if r["verdict_reason"] == "cds_intron" and r["i5"]]
        lbl5 = [str(i) for i in range(1, 7)] + ["7+"]
        v5 = [[r["lfc"] for r in rows if r["i5"] == i] for i in range(1, 7)]
        v5.append([r["lfc"] for r in rows if r["i5"] > 6])
        # i3 = n_exons - feat_idx, so the LAST intron has i3 == 1; there is no 0
        lbl3 = ["last"] + [str(i) for i in range(2, 6)]
        v3 = [[r["lfc"] for r in rows if r["i3"] == i] for i in range(1, 6)]
        x5 = list(range(len(v5)))
        x3 = [len(v5) + 1 + i for i in range(len(v3))]
        ax.bar(x5, [med(v) for v in v5], color=ORANGE, alpha=0.85,
               edgecolor="#fcfcfb", linewidth=1.2, label="intron number from the 5' end")
        ax.bar(x3, [med(v) for v in v3], color=BLUE, alpha=0.85,
               edgecolor="#fcfcfb", linewidth=1.2, label="intron number from the 3' end")
        for x, v in list(zip(x5, v5)) + list(zip(x3, v3)):
            ax.text(x, 0.04, "%d" % len(v), ha="center", va="bottom", fontsize=6,
                    color=INK2, rotation=90)
        ax.set_xticks(x5 + x3)
        ax.set_xticklabels(lbl5 + lbl3, fontsize=8)
        ax.axhline(0, color=INK2, linewidth=0.9)
        ax.set_title(lib)
        ax.set_xlabel("intron (labels = n)")
        tidy(ax)
    axes[0].set_ylabel("median log2FC")
    axes[0].legend(fontsize=8, ncol=2, loc="lower left", bbox_to_anchor=(0, 1.10))
    save(fig, "fig10_intron_index.png")


def fig11b_ptc_from_atg(data):
    """Median response against how far the PTC sits from the start codon, in nucleotides."""
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.2), sharey=True)
    for ax, (lib, _) in zip(axes, LIBS):
        rows = [r for r in expressed(data[lib])
                if r["verdict_reason"] == "cds_intron" and r["atg_dist"] and r["atg_dist"] > 0]
        binned(ax, [r["atg_dist"] for r in rows], [r["lfc"] for r in rows], 12, ORANGE)
        ax.axvspan(0, 200, color=GREY, alpha=0.13, linewidth=0)
        ax.text(45, ax.get_ylim()[0] + 0.15, "first\n~200 nt", fontsize=7.5, color=INK2)
        ax.axhline(0, color=INK2, linewidth=0.9, linestyle=":")
        ax.set_xscale("log")
        ax.set_title("%s (n=%d)" % (lib, len(rows)))
        ax.set_xlabel("distance from the start codon to the PTC, nt")
        tidy(ax)
    axes[0].set_ylabel("median log2FC (band = interquartile range)")
    save(fig, "fig11b_ptc_from_atg.png")


def fig11c_positional_variables(data):
    """Which positional variable orders the data best. Each is split into 8 equal-count bins, so
    the x axes are directly comparable and only the vertical spread matters."""
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.2), sharey=True)
    spec = [("distance ATG -> PTC", "atg_dist", ORANGE),
            ("insertion point in the transcript", "relpos", BLUE),
            ("distance PTC -> last junction", "dist_num", AQUA)]
    for ax, (lib, _) in zip(axes, LIBS):
        rows = [r for r in expressed(data[lib]) if r["verdict_reason"] == "cds_intron"]
        for r in rows:
            r["dist_num"] = (int(r["dist_to_last_junction"])
                             if r["dist_to_last_junction"] else None)
        for label, key, col in spec:
            s = sorted((r[key], r["lfc"]) for r in rows if r.get(key) is not None)
            nb, step = 8, len(s) / 8
            xs, ys = [], []
            for b in range(nb):
                chunk = s[int(b * step):int((b + 1) * step)]
                if len(chunk) < 25:
                    continue
                xs.append(b + 1)
                ys.append(med([y for _, y in chunk]))
            ax.plot(xs, ys, color=col, linewidth=2, marker="o", markersize=5,
                    markeredgecolor="#fcfcfb", markeredgewidth=0.8,
                    label="%s (spread %.2f)" % (label, max(ys) - min(ys)))
        ax.axhline(0, color=INK2, linewidth=0.9, linestyle=":")
        ax.set_title(lib)
        tidy(ax)
        ax.legend(fontsize=7.5, loc="upper left")
    axes[0].set_ylabel("median log2FC")
    fig.supxlabel("8 equal-sized groups of barcodes: 1 = smallest values, 8 = largest",
                  fontsize=9, y=-0.02)
    save(fig, "fig11c_positional_variables.png")


def fig11d_position_control(data):
    """The control that matters: does the SAME positional trend appear where no PTC is created?

    Fixed position bands, not equal-count bins, so the two groups are compared over the same x --
    the no-PTC group is almost entirely 5'-proximal (its introns are 5'UTR introns), so above
    ~0.15 it simply runs out of barcodes and that has to be visible rather than hidden.
    """
    edges = [(0, .05), (.05, .10), (.10, .15), (.15, .30), (.30, 1.01)]
    labels = ["0–0.05", "0.05–0.10", "0.10–0.15", "0.15–0.30", "0.30–1.0"]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.2), sharey=True)
    for ax, (lib, _) in zip(axes, LIBS):
        rows = [r for r in expressed(data[lib]) if r["relpos"] is not None]
        w = 0.38
        for j, (label, sel, col) in enumerate(
                (("PTC created", lambda r: r["verdict_reason"] == "cds_intron", ORANGE),
                 ("no PTC (UTR) - control", lambda r: r["verdict_reason"] == "utr", BLUE))):
            for i, (lo, hi) in enumerate(edges):
                v = [r["lfc"] for r in rows if sel(r) and lo <= r["relpos"] < hi]
                x = i + (j - 0.5) * w
                if len(v) < 25:
                    ax.text(x, 0.05, "n=%d" % len(v), ha="center", va="bottom", fontsize=6,
                            color=INK2, rotation=90)
                    continue
                ax.bar([x], [med(v)], width=w * 0.9, color=col, alpha=0.85,
                       edgecolor="#fcfcfb", linewidth=1.2, label=label if i == 0 else None)
                ax.text(x, med(v) + (0.05 if med(v) >= 0 else -0.17), "%d" % len(v),
                        ha="center", fontsize=6, color=INK2)
        ax.axhline(0, color=INK2, linewidth=0.9)
        ax.set_xticks(range(len(edges)))
        ax.set_xticklabels(labels, fontsize=7.5)
        ax.set_title(lib)
        ax.set_xlabel("insertion point in the transcript (labels = n)")
        tidy(ax)
        ax.legend(fontsize=7.5, loc="upper left")
    axes[0].set_ylabel("median log2FC")
    save(fig, "fig11d_position_control.png")


def fig11_relpos(data):
    """Effect against POSITION IN THE mRNA -- not along the genomic locus.

    The x variable is relpos = tx_pos / tx_len, both in transcript coordinates. The genomic
    equivalent exists as relgene and is not used here; fig01c shows the two side by side, one per
    row. For an insertion into an intron the transcript position is the splice junction the
    cassette goes into, because the intron itself is not part of the mRNA.

    Restricted to PTCs the last-exon rule and the 50-nt rule both let through, so "PTC" here means
    a stop actually expected to trigger NMD. The control group is still `not has_ptc`, the wider
    flag: a PTC that fails the 50-nt rule has a PTC and does not belong among the barcodes that
    make none.

    Each line is labelled at its own end rather than in a legend, so a label cannot land on the
    data - which is what the legend box used to do."""
    fig, axes = plt.subplots(1, 2, figsize=(COL2, 66 * MM), sharey=True)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.90, bottom=0.16, wspace=0.06)
    series = (("PTC expected to trigger NMD", lambda r: r["nmd_ptc"], ORANGE),
              ("no PTC created", lambda r: r["in_gene"] == "1" and not r["has_ptc"], BLUE))
    print("  fig11_relpos: barcodes per series")
    ends = []
    for ax, (lib, _) in zip(axes, LIBS):
        rows = [r for r in expressed(data[lib]) if r["relpos"] is not None]
        for label, sel, colour in series:
            s = [r for r in rows if sel(r)]
            if len(s) <= 60:
                continue
            cx, cm, _ = binned(ax, [r["relpos"] for r in s], [r["lfc"] for r in s], 10, colour,
                               band=False)
            ends.append((ax, "%s   n = %s" % (label, format(len(s), ",")), colour))
            print("    %-8s %-28s n = %s" % (lib, label, format(len(s), ",")))
        ax.axhline(0, color=INK2, linewidth=HAIR, linestyle=(0, (1, 2.5)), zorder=1)
        ax.set_title(DISPLAY[lib], fontsize=FS, color=INK, pad=4)
        journal(ax)
        ax.set_xlim(-0.02, 1.02)
        ticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
        ax.set_xticks(ticks)
        ax.set_xticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1"])
        x_rule(ax, ticks)
    # Labels go into a band above the data that is empty BY CONSTRUCTION - the y limit is set from
    # the data first and the band added on top. Labelling the end of each line was tried and is
    # wrong here: a right-aligned label runs back along its own rising curve.
    ys = [y for ax in axes for l in ax.get_lines() if len(l.get_ydata()) > 2
          for y in l.get_ydata()]
    ylo, yhi = min(ys), max(ys)
    rng = yhi - ylo
    axes[0].set_ylim(ylo - 0.10 * rng, yhi + 0.46 * rng)
    for ax in axes:
        for k, (_, text, colour) in enumerate([e for e in ends if e[0] is ax]):
            ax.text(0.015, 0.985 - 0.085 * k, text, transform=ax.transAxes,
                    ha="left", va="top", fontsize=FS_SMALL, color=colour)
    axes[0].set_ylabel("median $\\mathregular{log_2}$ fold change (CHX / EtOH)", fontsize=FS, color=INK)
    fig.text(0.53, 0.035, "position in the mRNA: 0 = 5\u2032 end, 1 = 3\u2032 end",
             ha="center", fontsize=FS, color=INK)
    save(fig, "fig11_relpos.png", vector=True)


def fig12_region_loc(data):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.0), sharey=True)
    regions = [("5'UTR", "5utr"), ("CDS", "cds"), ("3'UTR", "3utr")]
    for ax, (lib, _) in zip(axes, LIBS):
        rows = expressed(data[lib])
        w = 0.38
        # A bar is drawn only for n >= MIN_N. Without this the 3'UTR-intronic group (n=5 and 7
        # expressed barcodes) sets the whole y-scale at a median of +5 and flattens everything else.
        MIN_N = 20
        for j, (loc, col) in enumerate((("intronic", BLUE), ("exonic", ORANGE))):
            for i, (label, key) in enumerate(regions):
                v = [r["lfc"] for r in rows if r["region"] == key and r["loc"] == loc]
                x = i + (j - 0.5) * w
                if len(v) < MIN_N:
                    ax.text(x, 0.06, "n=%d\ntoo few" % len(v), ha="center", va="bottom",
                            fontsize=6.5, color=INK2)
                    continue
                h = med(v)
                ax.bar([x], [h], width=w * 0.92, color=col, alpha=0.85, edgecolor="#fcfcfb",
                       linewidth=1.5,
                       label=("into an intron" if loc == "intronic" else "into an exon") if i == 1 else None)
                ax.text(x, h + (0.04 if h >= 0 else -0.14), "n=%d" % len(v), ha="center",
                        fontsize=6.5, color=INK2)
        ax.axhline(0, color=INK2, linewidth=0.9)
        ax.set_xticks(range(len(regions)))
        ax.set_xticklabels([r[0] for r in regions])
        ax.set_title(lib)
        tidy(ax)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("median log2FC")
    save(fig, "fig12_region_loc.png")


def fig13_ma(data):
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.4), sharey=True)
    for ax, (lib, _) in zip(axes, LIBS):
        rows = [r for r in expressed(data[lib]) if r["base"] and r["base"] > 0]
        for label, sel, col in (("no PTC created / outside a gene",
                                 lambda r: not r["has_ptc"], BLUE),
                                ("PTC created", lambda r: r["has_ptc"], ORANGE)):
            s = [r for r in rows if sel(r)]
            ax.scatter([r["base"] for r in s], [r["lfc"] for r in s], s=2.5, alpha=0.22,
                       color=col, linewidths=0, rasterized=True,
                       label="%s (n=%d)" % (label, len(s)))
            binned(ax, [r["base"] for r in s], [r["lfc"] for r in s], 10, col, band=False)
        ax.axhline(0, color=INK2, linewidth=0.9, linestyle=":")
        ax.set_xscale("log")
        ax.set_ylim(-4, 5)
        ax.set_xlim(quantile([r["base"] for r in rows], 0.001) * 0.8,
                    quantile([r["base"] for r in rows], 0.999) * 1.3)
        ax.set_title(lib)
        ax.set_xlabel("mean barcode abundance (CHXbaseMean)")
        tidy(ax)
        leg = ax.legend(fontsize=7.5, loc="upper right", markerscale=4)
        for h in leg.legend_handles:
            h.set_alpha(1)
    axes[0].set_ylabel("log2 fold change (CHX / EtOH)")
    save(fig, "fig13_ma.png")


def fig14_host_expr(data):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.2), sharey=True)
    for ax, (lib, _) in zip(axes, LIBS):
        rows = [r for r in expressed(data[lib]) if r["host_exp"] is not None]
        for label, sel, col in (("PTC created", lambda r: r["has_ptc"], ORANGE),
                                ("no PTC created (UTR)",
                                 lambda r: r["in_gene"] == "1" and not r["has_ptc"], BLUE)):
            s = [r for r in rows if sel(r)]
            if len(s) > 60:
                binned(ax, [r["host_exp"] for r in s], [r["lfc"] for r in s], 10, col,
                       label="%s (n=%d)" % (label, len(s)), band=False)
        ax.axhline(0, color=INK2, linewidth=0.9, linestyle=":")
        ax.set_title(lib)
        ax.set_xlabel("host gene expression, log10")
        tidy(ax)
        ax.legend(fontsize=7.5, loc="upper left")
    axes[0].set_ylabel("median log2FC")
    save(fig, "fig14_host_expr.png")


def fig15_chromatin(data):
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharex=True)
    # one state order for both panels, or the rows would not line up between libraries
    combined = Counter()
    for lib, _ in LIBS:
        combined.update(r["chromhmm"] for r in expressed(data[lib]))
    order = [s for s, n in combined.most_common() if s][:8]
    for ax, (lib, _) in zip(axes, LIBS):
        rows = expressed(data[lib])
        states = order
        ys = list(range(len(states)))[::-1]
        for label, sel, col, off in (("PTC created", lambda r: r["has_ptc"], ORANGE, 0.18),
                                     ("no PTC created / outside a gene",
                                      lambda r: not r["has_ptc"], BLUE, -0.18)):
            for y, st in zip(ys, states):
                v = [r["lfc"] for r in rows if r["chromhmm"] == st and sel(r)]
                if len(v) < 10:
                    continue
                ax.plot([quantile(v, .25), quantile(v, .75)], [y + off] * 2, color=col,
                        linewidth=2.5, alpha=0.4, solid_capstyle="round")
                ax.plot([med(v)], [y + off], "o", markersize=6, color=col,
                        markeredgecolor="#fcfcfb", markeredgewidth=1,
                        label=label if st == states[0] else None)
        ax.axvline(0, color=INK2, linewidth=0.9, linestyle=":")
        ax.set_yticks(ys)
        ax.set_yticklabels(states, fontsize=8)
        ax.set_xlabel("median log2FC, bar = interquartile range")
        ax.set_title(lib)
        tidy(ax, xgrid=True, ygrid=False)
    axes[0].legend(fontsize=8, ncol=2, loc="lower left", bbox_to_anchor=(0, 1.08))
    save(fig, "fig15_chromatin.png")


def fig16_reproducibility(data):
    # the offsets are hand-placed: several groups land almost on top of each other near zero
    spec = [("outside a gene / antisense", lambda r: r["category"] == "2", GREY, (10, 6)),
            ("UTR", lambda r: r["verdict_reason"] == "utr", BLUE, (10, -14)),
            ("earlier coding exon", lambda r: r["verdict_reason"] == "internal_exon",
             ORANGE, (12, 0)),
            ("CDS intron, dist > 50", lambda r: r["dist_to_last_junction"]
             and int(r["dist_to_last_junction"]) > 50, ORANGE, (-24, -30)),
            ("CDS intron, 0 < dist <= 50", lambda r: r["dist_to_last_junction"]
             and 0 < int(r["dist_to_last_junction"]) <= 50, AQUA, (-14, -20)),
            ("stop in the last exon", lambda r: r["dist_to_last_junction"]
             and int(r["dist_to_last_junction"]) <= 0, AQUA, (-130, -6))]
    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    a, b = expressed(data["5ss_noP"]), expressed(data["CPA_noP"])
    lim = [-0.7, 3.6]
    ax.plot(lim, lim, color=GREY, linewidth=1, linestyle="--", zorder=1)
    for label, sel, col, off in spec:
        x = med([r["lfc"] for r in a if sel(r)])
        y = med([r["lfc"] for r in b if sel(r)])
        ax.plot([x], [y], "o", markersize=10, color=col, markeredgecolor="#fcfcfb",
                markeredgewidth=1.4, zorder=3)
        ax.annotate(label, (x, y), textcoords="offset points", xytext=off,
                    fontsize=8, color=INK)
    ax.axhline(0, color=INK2, linewidth=0.8, linestyle=":")
    ax.axvline(0, color=INK2, linewidth=0.8, linestyle=":")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("median log2FC in 5ss_noP")
    ax.set_ylabel("median log2FC in CPA_noP")
    ax.set_title("one dot = one prediction group")
    tidy(ax, xgrid=True)
    save(fig, "fig16_reproducibility.png")


def main(argv=()):
    """No arguments redraws everything; any arguments are name fragments and only the figure
    functions matching one of them run (`make_figures.py fig04`), so a single figure can be
    re-exported without rewriting the other twenty PNGs."""
    os.makedirs(PICS, exist_ok=True)
    figures = (fig01_composition, fig01b_loc_counts, fig01c_position_landscape,
               fig02_ecdf, fig03_verdict, fig04_ptc, fig05_dist_hist,
               fig07_matched, fig08_ptc_source, fig09_frame, fig10_intron_index,
               fig11_relpos, fig11b_ptc_from_atg, fig11c_positional_variables,
               fig11d_position_control, fig12_region_loc, fig13_ma, fig14_host_expr,
               fig15_chromatin, fig16_reproducibility)
    if argv:
        figures = tuple(f for f in figures if any(a in f.__name__ for a in argv))
        if not figures:
            print("no figure matches %s" % " ".join(argv))
            return 1
    data = load()
    print("figures -> analysis/pics/")
    for f in figures:
        f(data)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
