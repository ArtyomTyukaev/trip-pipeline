#!/usr/bin/env python3
"""
Figures for the IntronTRIP NMD-inhibition comparison -- the core of analysis/SUMMARY.md, redone on
the eight IntronTRIP constructs, in the style of analysis/make_figures.py.

What is different from the single_combined analysis, and how each figure handles it:

- **Two isoforms.** Igor's tables carry NMDi/Control separately for molecules that excised the
  reporter intron (spliced isoform) and molecules that kept it (unspliced isoform). The verdict is
  the same for both on this construct, so every figure is drawn once per isoform.
- **Eight constructs, each its own DESeq2 run.** Figures that pool constructs (verdict branches,
  position, intron index, the paired scatter, the "all 8" effect) pool the log2 fold changes as
  they are. The constructs' out-of-gene medians differ by up to ~0.4 (spliced), small next to the
  effects shown, and every comparison drawn is between groups of the same pooled set.
- **No 50-nt-threshold or PTC-source figures.** On this construct the threshold decides no row and
  the stop always comes from the insert (PLAN_INTRON_NMD.md section 10), so there is nothing to split.

Groups, as in stage 4: category 2 = outside a gene or antisense (the baseline), 3 = in a gene, no
NMD expected, 4 = NMD expected. Positional figures use `cds_intron` rows, with the PTC position of
the isoform being measured. Boxes: interquartile range; whiskers 5th-95th percentile, no fliers,
as in analysis/make_figures.py. Tests: two-sided Mann-Whitney, category 4 against category 3, one
comparison per construct and isoform, so no multiple-testing correction; exact p in the console
output and in summary_by_category.tsv. Nothing explanatory is written on the figures.

    module load matplotlib/3.7.0-gfbf-2022b
    python3 intron_nmdi/make_figures.py [figure-name-fragment ...]

in   intron_nmdi/output/nmd_nmdi_<construct>.tsv     -- stage 6
     work/transcripts.tsv                            -- tx_len and the CDS start
out  intron_nmdi/pics/*.png, intron_nmdi/summary_by_category.tsv
"""

import csv
import math
import os
import random
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
NMD = os.path.dirname(HERE)
# analysis/ in the working tree, single_combined/analysis/ in the trip-pipeline repository
for _d in (os.path.join(NMD, "analysis"), os.path.join(NMD, "single_combined", "analysis")):
    if os.path.isfile(os.path.join(_d, "make_figures.py")):
        sys.path.insert(0, _d)
        break
# Palette, fonts, rcParams and the box/bin/test helpers of the single_combined figures, so the two
# analyses read as one. Importing it registers Arial and sets the style; it draws nothing.
import make_figures as mf  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

PICS = os.path.join(HERE, "pics")
LIBS = ["A1", "B1", "B3-6", "B3-9", "C1", "D1", "E1", "H1"]
FORMS = [("spliced", "Spliced"), ("unspliced", "Unspliced")]
GREY, GREEN, VIOLET, OLIVE = mf.GREY, mf.BLUE, mf.ORANGE, mf.AQUA
INK, INK2 = mf.INK, mf.INK2
CATS = [("2", "outside a gene / antisense", GREY),
        ("3", "in a gene, no NMD expected", GREEN),
        ("4", "NMD expected", VIOLET)]
# Igor's libraries call the conditions NMDi and Control; they are CHX and EtOH, as in 5ss/CPA
# (confirmed 2026-10-16), so the labels are those of analysis/make_figures.py
LFC = "log2 fold change (CHX / EtOH)"
MLFC = "median log2 fold change (CHX / EtOH)"
ATG_EDGES = [0, 100, 200, 400, 800, 1600, 3200, 10 ** 9]
BOOT = 400                       # bootstrap resamples for the effect-size intervals
med = mf.med


def num(v):
    return None if v in ("", "NA", "NaN", "nan") else float(v)


def iso(f):
    return "%s isoform" % f


def title(fig, text, y=1.0):
    """One plain title per figure saying what it shows; panels keep the isoform as their title,
    as analysis/make_figures.py keeps the library."""
    fig.suptitle(text, fontsize=11, color=INK, y=y)


# ----------------------------------------------------------------------------- data


def load():
    """Expressed rows only (some fold change present), plus the category counts of every row.

    The tables hold every barcode of the BED, ~1.6 M rows of which ~95% are not expressed; only
    the composition figure needs those, and it needs only their count."""
    tx = {}
    with open(os.path.join(NMD, "work", "transcripts.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tx[r["tx_id"]] = (int(r["tx_len"]),
                              int(r["cds_tx_start"]) if r["cds_tx_start"] else None)
    data, comp = {}, {}
    for lib in LIBS:
        rows, cnt = [], Counter()
        with open(os.path.join(HERE, "output", "nmd_nmdi_%s.tsv" % lib)) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                for f, _ in FORMS:
                    cnt[(f, r["category_%s" % f])] += 1
                lfc = {f: num(r["%s_log2FoldChange" % R]) for f, R in FORMS}
                inter = num(r["SplicingInteraction_log2FoldChange"])
                if lfc["spliced"] is None and lfc["unspliced"] is None and inter is None:
                    continue
                d = {"lib": lib, "lfc": lfc, "inter": inter,
                     "cat": {f: r["category_%s" % f] for f, _ in FORMS},
                     "in_gene": r["in_gene"] == "1", "verdict": r["verdict_reason"],
                     "nmd": r["nmd_expected"], "loc": r["loc"]}
                tl, cs = tx.get(r["tx_id"], (None, None))
                d["relpos"] = int(r["tx_pos"]) / tl if r["tx_pos"] and tl else None
                d["has_ptc"], d["nmd_ptc"], d["atg"], d["dist"] = {}, {}, {}, {}
                for f, _ in FORMS:
                    # same two flags as analysis/make_figures.py: has_ptc is the wide one, and
                    # a no-PTC control is built from `not has_ptc`, never from `not nmd_ptc`
                    hp = ((d["verdict"] == "cds_intron"
                           and r["ptc_source_%s" % f] not in ("", "none"))
                          or d["verdict"] in ("internal_exon", "last_exon"))
                    d["has_ptc"][f] = hp
                    d["nmd_ptc"][f] = hp and d["nmd"] == "1"
                    p = r["ptc_tx_pos_%s" % f]
                    d["atg"][f] = int(p) - cs if p and cs is not None else None
                    q = r["dist_to_last_junction_%s" % f]
                    d["dist"][f] = int(q) if q else None
                fi = int(r["feat_idx"]) if r["feat_idx"] else None
                d["i5"] = fi if fi and r["loc"] == "intronic" else None
                d["i3"] = int(r["n_exons"]) - fi if fi and r["loc"] == "intronic" else None
                rows.append(d)
        data[lib], comp[lib] = rows, cnt
    return data, comp


def by_cat(rows, f):
    """category -> fold changes of one isoform."""
    out = {c: [] for c, _, _ in CATS}
    for r in rows:
        v, c = r["lfc"][f], r["cat"][f]
        if v is not None and c in out:
            out[c].append(v)
    return out


def inter_by_cat(rows):
    """category -> spliced/unspliced-ratio fold changes. The ratio belongs to the barcode, not to
    one isoform; its group is the barcode's category, the same in both isoforms whenever expressed."""
    out = {c: [] for c, _, _ in CATS}
    for r in rows:
        c = r["cat"]["spliced"] if r["cat"]["spliced"] != "1" else r["cat"]["unspliced"]
        if r["inter"] is not None and c in out:
            out[c].append(r["inter"])
    return out


def save(fig, name):
    os.makedirs(PICS, exist_ok=True)
    fig.savefig(os.path.join(PICS, name), facecolor=fig.get_facecolor())
    print("  pics/%s" % name)
    plt.close(fig)


def boot_diff(a, b, rnd):
    """Median(a) - median(b) and its 95% percentile-bootstrap interval."""
    d = []
    for _ in range(BOOT):
        d.append(med(rnd.choices(a, k=len(a))) - med(rnd.choices(b, k=len(b))))
    d.sort()
    return med(a) - med(b), d[int(0.025 * BOOT)], d[int(0.975 * BOOT) - 1]


def grouped_boxes(ax, per_lib, colors, width=0.26):
    """Vertical boxes, one cluster per construct. Returns (box x positions, whisker extent).

    A group under 5 barcodes keeps its slot but is drawn as a bare median dot, as box_groups does."""
    k = len(colors)
    vals, pos, cols = [], [], []
    for i, groups in enumerate(per_lib):
        for j, v in enumerate(groups):
            vals.append(v)
            pos.append(i + (j - (k - 1) / 2) * width)
            cols.append(colors[j])
    big = [v if len(v) >= 5 else [] for v in vals]
    bp = ax.boxplot([v if v else [0] for v in big], positions=pos, widths=width * 0.82,
                    showfliers=False, whis=(5, 95), patch_artist=True, zorder=3,
                    medianprops=dict(color=INK, linewidth=1.4),
                    whiskerprops=dict(color=INK2, linewidth=0.8),
                    capprops=dict(color=INK2, linewidth=0.8),
                    boxprops=dict(linewidth=0.7, edgecolor=INK2))
    for i, (patch, c) in enumerate(zip(bp["boxes"], cols)):
        if big[i]:
            patch.set_facecolor(c)
            patch.set_alpha(0.45)
        else:
            patch.set_visible(False)
            bp["medians"][i].set_visible(False)
            for kk in ("whiskers", "caps"):
                bp[kk][2 * i].set_visible(False)
                bp[kk][2 * i + 1].set_visible(False)
            if vals[i]:
                ax.plot([pos[i]], [med(vals[i])], "o", markersize=5, color=c, zorder=4)
    ys = [y for w in bp["whiskers"] if w.get_visible() for y in w.get_ydata()]
    return pos, (min(ys), max(ys))


def box_panel_labels(ax, per_lib, pos, lo, hi):
    """The group sizes below the lowest whisker of the panel, clear of the x axis. Returns the
    category 4 vs 3 p values in construct order, for the console only: no test is drawn, on
    request."""
    rng = hi - lo
    ps = []
    for i in range(len(per_lib)):
        g2, g3, g4 = per_lib[i]
        ps.append(mf.mannwhitney(g4, g3))
        for j, v in enumerate(per_lib[i]):
            ax.text(pos[3 * i + j], lo - 0.04 * rng, format(len(v), ","), rotation=90,
                    ha="center", va="top", fontsize=5, color=INK2)
    return ps


def legend_boxes(ax):
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=c, alpha=0.45, edgecolor=INK2)
               for _, _, c in CATS]
    ax.legend(handles, [l for _, l, _ in CATS], ncol=3, fontsize=7.5, loc="lower left",
              bbox_to_anchor=(0, 1.10))


def fixed_bins(xs, ys, edges):
    out = []
    for lo, hi in zip(edges, edges[1:]):
        v = [y for x, y in zip(xs, ys) if lo <= x < hi]
        out.append((lo, hi, len(v), med(v) if v else float("nan")))
    return out


def octile_spread(pairs):
    s = sorted(pairs)
    step = len(s) / 8
    meds = [med([y for _, y in s[int(b * step):int((b + 1) * step)]]) for b in range(8)
            if int((b + 1) * step) - int(b * step) >= 25]
    return (max(meds) - min(meds)) if meds else float("nan"), meds


def thin_line(ax, xs, ys, nbins, color):
    """Equal-count-bin medians as a plain thin line - the per-construct traces under a pooled one."""
    pts = sorted(zip(xs, ys))
    if len(pts) < nbins * 20:
        nbins = max(3, len(pts) // 20)
    step = len(pts) / nbins
    cx, cm = [], []
    for b in range(nbins):
        chunk = pts[int(b * step):int((b + 1) * step)]
        if len(chunk) >= 10:
            cx.append(med([p[0] for p in chunk]))
            cm.append(med([p[1] for p in chunk]))
    ax.plot(cx, cm, color=color, linewidth=0.8, alpha=0.55, zorder=2)


def pooled(data, sel):
    return [r for lib in LIBS for r in data[lib] if sel(r)]


# ----------------------------------------------------------------------------- figures


def fig01_composition(data, comp):
    """How many barcodes of each construct reach each group, per isoform. Expressed barcodes only:
    the not-expressed ones are 92-97% of every library and would leave the rest as slivers; their
    number is printed at the end of each bar instead."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4), sharex=True, sharey=True)
    top = max(sum(comp[lib][(f, k)] for k, _, _ in CATS) for lib in LIBS for f, _ in FORMS)
    for ax, (f, R) in zip(axes, FORMS):
        for i, lib in enumerate(LIBS):
            c = comp[lib]
            left = 0
            for key, label, col in CATS:
                w = c[(f, key)]
                ax.barh([i], [w], left=left, color=col, alpha=0.85, height=0.62,
                        edgecolor="#fcfcfb", linewidth=1.2,
                        label=label if (i == 0 and ax is axes[0]) else None)
                left += w
            ne = c[(f, "1")]
            ax.text(left + 250, i, "%s  (%s not expressed)" % (format(left, ","), format(ne, ",")),
                    va="center", fontsize=6.5, color=INK2)
        ax.set_xlim(0, top * 1.55)
        ax.set_yticks(range(len(LIBS)))
        ax.set_yticklabels(LIBS)
        ax.set_title(iso(f))
        ax.set_xlabel("number of barcodes")
        mf.tidy(ax, xgrid=True, ygrid=False)
    axes[0].invert_yaxis()           # once: the y axis is shared, a second call would undo it
    axes[0].legend(ncol=3, loc="lower left", bbox_to_anchor=(0, 1.08), fontsize=7.5)
    title(fig, "1. Barcodes with a measured fold change, by group", y=1.08)
    save(fig, "fig01_composition.png")
    print("  fig01: expressed barcodes per category (2 / 3 / 4), not expressed")
    for lib in LIBS:
        print("    %-5s" % lib + "".join(
            "   %s %5d / %5d / %5d, %7d" % (f, comp[lib][(f, "2")], comp[lib][(f, "3")],
                                            comp[lib][(f, "4")], comp[lib][(f, "1")])
            for f, _ in FORMS))


def fig02_ecdf(data, comp):
    fig, axes = plt.subplots(2, 8, figsize=(13, 4.4), sharex=True, sharey=True)
    for row, (f, R) in enumerate(FORMS):
        for col, lib in enumerate(LIBS):
            ax = axes[row][col]
            groups = by_cat(data[lib], f)
            for k, (key, label, c) in enumerate(CATS):
                v = sorted(groups[key])
                if not v:
                    continue
                ax.plot(v, [i / len(v) for i in range(1, len(v) + 1)], color=c, linewidth=1.5,
                        label=label)
                ax.text(0.97, 0.05 + 0.09 * k, "n=%s" % format(len(v), ","), color=c,
                        transform=ax.transAxes, ha="right", fontsize=5.5)
            ax.axvline(0, color=INK2, linewidth=0.8, linestyle=":")
            ax.set_xlim(-4, 4)
            if row == 0:
                ax.set_title(lib)
            mf.tidy(ax)
            ax.tick_params(labelsize=7)
        axes[row][0].set_ylabel("%s\nfraction of barcodes <= x" % iso(f), fontsize=8)
    fig.supxlabel(LFC, fontsize=9, y=0.0)
    axes[0][0].legend(fontsize=7, loc="lower left", bbox_to_anchor=(0, 1.22), ncol=3)
    save(fig, "fig02_ecdf_category.png")


def fig03_box_category(data, comp):
    """The main result: the three groups per construct, both isoforms on one y scale, category 4
    vs 3 tested."""
    fig, axes = plt.subplots(2, 1, figsize=(9, 6.8), sharex=True, sharey=True)
    panels = []
    for ax, (f, R) in zip(axes, FORMS):
        per_lib = []
        for lib in LIBS:
            g = by_cat(data[lib], f)
            per_lib.append([g["2"], g["3"], g["4"]])
        pos, span = grouped_boxes(ax, per_lib, [c for _, _, c in CATS])
        panels.append((ax, f, per_lib, pos, span))
    lo = min(p[4][0] for p in panels)
    hi = max(p[4][1] for p in panels)
    rng = hi - lo
    axes[0].set_ylim(lo - 0.24 * rng, hi + 0.05 * rng)
    print("  fig03: median log2FC per group, and category 4 vs 3 (two-sided Mann-Whitney)")
    for ax, f, per_lib, pos, _ in panels:
        ps = box_panel_labels(ax, per_lib, pos, lo, hi)
        for lib, (g2, g3, g4), p in zip(LIBS, per_lib, ps):
            print("    %-9s %-5s  baseline %+.2f (n=%d)  no NMD %+.2f (n=%d)  NMD %+.2f (n=%d)"
                  "   diff %+.2f   p = %.2g"
                  % (f, lib, med(g2), len(g2), med(g3), len(g3), med(g4), len(g4),
                     med(g4) - med(g3), p if p is not None else float("nan")))
        ax.axhline(0, color=INK2, linewidth=0.8, linestyle=":")
        ax.set_title(iso(f))
        ax.set_ylabel(LFC)
        mf.tidy(ax)
    axes[-1].set_xticks(range(len(LIBS)))
    axes[-1].set_xticklabels(LIBS)
    legend_boxes(axes[0])
    title(fig, "3. CHX response by NMD prediction, per construct", y=1.03)
    save(fig, "fig03_box_category.png")


def fig04_verdict_pooled(data, comp):
    """Which branch of the rule decided the barcode, all constructs pooled."""
    spec = [("outside a gene / antisense", lambda r, f: r["cat"][f] == "2", GREY),
            ("UTR - no PTC created", lambda r, f: r["verdict"] == "utr", GREEN),
            ("exonic, last coding exon\nrule: escape", lambda r, f: r["verdict"] == "last_exon",
             GREEN),
            ("exonic, earlier coding exon\nrule: NMD", lambda r, f: r["verdict"] == "internal_exon",
             VIOLET),
            ("intronic in the CDS\nrule: NMD", lambda r, f: r["verdict"] == "cds_intron", VIOLET),
            ("intronic, no PTC -- the\nframeshift made none",
             lambda r, f: r["verdict"] == "no_ptc_frameshift", GREY)]
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.9))
    spans = []
    print("  fig04: pooled medians by branch")
    for ax, (f, R) in zip(axes, FORMS):
        rows = [r for lib in LIBS for r in data[lib] if r["lfc"][f] is not None]
        groups, cols = [], []
        for label, sel, c in spec:
            v = [r["lfc"][f] for r in rows if sel(r, f)]
            if v:
                groups.append((label, v))
                cols.append(c)
                print("    %-9s %-40s n = %6d  median %+.2f"
                      % (f, label.replace("\n", " "), len(v), med(v)))
        spans.append(mf.box_groups(ax, groups, cols, xlim=None, ylabels=(ax is axes[0]),
                                   xlabel=LFC))
        ax.set_title(iso(f))
        ax.tick_params(axis="y", labelsize=7.5)
    lo, hi = min(s[0] for s in spans), max(s[1] for s in spans)
    for ax in axes:
        ax.set_xlim(lo - 0.05 * (hi - lo), hi + 0.05 * (hi - lo))
    title(fig, "4. CHX response by branch of the NMD rule, all 8 constructs", y=1.03)
    save(fig, "fig04_verdict_pooled.png")


def fig05_ptc_from_atg(data, comp):
    """Response against the ribosome's distance from the start codon to the PTC: each construct
    (thin) and all eight pooled (thick). In single_combined this was a step at ~200 nt, the
    AUG-proximity escape."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
    print("  fig05: median log2FC of cds_intron rows by distance ATG -> PTC")
    hdr = "".join("%11s" % ("%d-%s" % (lo, hi if hi < 10 ** 9 else "")) for lo, hi
                  in zip(ATG_EDGES, ATG_EDGES[1:]))
    for ax, (f, R) in zip(axes, FORMS):
        print("    %s isoform%s" % (f, hdr))
        for lib in LIBS:
            rows = [r for r in data[lib] if r["verdict"] == "cds_intron"
                    and r["atg"][f] and r["atg"][f] > 0 and r["lfc"][f] is not None]
            thin_line(ax, [r["atg"][f] for r in rows], [r["lfc"][f] for r in rows], 8, GREY)
            fb = fixed_bins([r["atg"][f] for r in rows], [r["lfc"][f] for r in rows], ATG_EDGES)
            print("      %-6s" % lib + "".join("%11s" % ("%+.2f (%d)" % (m, n) if n >= 10 else
                                                    "n=%d" % n) for _, _, n, m in fb))
        rows = pooled(data, lambda r: r["verdict"] == "cds_intron" and r["atg"][f]
                      and r["atg"][f] > 0 and r["lfc"][f] is not None)
        mf.binned(ax, [r["atg"][f] for r in rows], [r["lfc"][f] for r in rows], 14, VIOLET,
                  label="all 8 constructs", band=False)
        fb = fixed_bins([r["atg"][f] for r in rows], [r["lfc"][f] for r in rows], ATG_EDGES)
        print("      %-6s" % "all" + "".join("%11s" % ("%+.2f (%d)" % (m, n)) for _, _, n, m in fb))
        ax.axvspan(1, 200, color=GREY, alpha=0.13, linewidth=0)
        ax.axhline(0, color=INK2, linewidth=0.9, linestyle=":")
        ax.set_xscale("log")
        ax.set_xlim(40, 1e4)
        ax.set_title("%s (n = %s)" % (iso(f), format(len(rows), ",")))
        ax.set_xlabel("distance from the start codon to the PTC, nt")
        mf.tidy(ax)
        ax.legend(handles=[Line2D([0], [0], color=VIOLET, linewidth=1.4, marker="o",
                                  markersize=3, label="all 8 constructs"),
                           Line2D([0], [0], color=GREY, linewidth=0.8, alpha=0.8,
                                  label="single constructs")],
                  fontsize=7, loc="upper left")
    axes[0].set_ylabel(MLFC)
    title(fig, "5. CHX response against the distance from the start codon to the PTC", y=1.03)
    save(fig, "fig05_ptc_from_atg.png")


def fig06_positional_variables(data, comp):
    """Which positional variable orders the response best: octile medians, all constructs pooled."""
    spec = [("distance from the start codon to the PTC", lambda r, f: r["atg"][f], VIOLET),
            ("position of the insertion in the mRNA", lambda r, f: r["relpos"], GREEN),
            ("distance from the PTC to the last exon junction", lambda r, f: r["dist"][f], OLIVE)]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4), sharey=True)
    print("  fig06: spread of octile medians (max - min); 'all' pooled, then per construct")
    for ax, (f, R) in zip(axes, FORMS):
        rows = pooled(data, lambda r: r["verdict"] == "cds_intron" and r["lfc"][f] is not None)
        for label, get, c in spec:
            pairs = [(get(r, f), r["lfc"][f]) for r in rows if get(r, f) is not None]
            spread, meds = octile_spread(pairs)
            ax.plot(range(1, len(meds) + 1), meds, color=c, linewidth=2, marker="o",
                    markersize=5, markeredgecolor="#fcfcfb", markeredgewidth=0.8,
                    label="%s (spread %.2f)" % (label, spread))
            per = []
            for lib in LIBS:
                pl = [(get(r, f), r["lfc"][f]) for r in data[lib] if r["verdict"] == "cds_intron"
                      and r["lfc"][f] is not None and get(r, f) is not None]
                per.append("%s %.2f" % (lib, octile_spread(pl)[0]))
            print("    %-9s %-34s all %.2f |  %s" % (f, label, spread, "  ".join(per)))
        ax.axhline(0, color=INK2, linewidth=0.9, linestyle=":")
        ax.set_title("%s (n = %s)" % (iso(f), format(len(rows), ",")))
        ax.set_xticks(range(1, 9))
        mf.tidy(ax)
        # the spliced lines fill the top of the shared y range, the unspliced ones its bottom
        ax.legend(fontsize=7, loc="lower right" if f == "spliced" else "upper left")
    axes[0].set_ylabel(MLFC)
    fig.supxlabel("barcodes sorted by the variable of each line and split into 8 equal groups: "
                  "1 = smallest values, 8 = largest", fontsize=9, y=-0.02)
    title(fig, "6. Which positional variable best orders the CHX response, all 8 constructs\n"
          "each line is one variable: the distance from the start codon to the PTC, the position "
          "of the insertion in the mRNA,\nor the distance from the PTC to the last exon junction",
          y=1.10)
    save(fig, "fig06_positional_variables.png")


def fig07_relpos(data, comp):
    """Position in the mature mRNA, NMD-expected PTC against no PTC created, all constructs pooled."""
    series = (("PTC expected to trigger NMD", lambda r, f: r["nmd_ptc"][f], VIOLET),
              ("no PTC created", lambda r, f: r["in_gene"] and not r["has_ptc"][f], GREEN))
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4), sharey=True)
    print("  fig07: barcodes per series")
    for ax, (f, R) in zip(axes, FORMS):
        rows = pooled(data, lambda r: r["relpos"] is not None and r["lfc"][f] is not None)
        for label, sel, c in series:
            s = [r for r in rows if sel(r, f)]
            mf.binned(ax, [r["relpos"] for r in s], [r["lfc"][f] for r in s], 10, c,
                      label="%s (n = %s)" % (label, format(len(s), ",")), band=False)
            print("    %-9s %-28s n = %s" % (f, label, format(len(s), ",")))
        ax.axhline(0, color=INK2, linewidth=0.9, linestyle=":")
        ax.set_xlim(-0.02, 1.02)
        ax.set_title(iso(f))
        ax.set_xlabel("position in the mRNA: 0 = 5′ end, 1 = 3′ end")
        mf.tidy(ax)
        ax.legend(fontsize=7, loc="lower right" if f == "spliced" else "upper left")
    axes[0].set_ylabel(MLFC)
    title(fig, "7. CHX response along the mRNA, all 8 constructs", y=1.03)
    save(fig, "fig07_relpos.png")


def fig08_intron_index(data, comp):
    """Which host intron the cassette went into, from the 5' and from the 3' end; pooled."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.3), sharey=True)
    tops, lows = [], []
    for ax, (f, R) in zip(axes, FORMS):
        rows = pooled(data, lambda r: r["verdict"] == "cds_intron" and r["i5"]
                      and r["lfc"][f] is not None)
        v5 = [[r["lfc"][f] for r in rows if r["i5"] == i] for i in range(1, 7)]
        v5.append([r["lfc"][f] for r in rows if r["i5"] > 6])
        v3 = [[r["lfc"][f] for r in rows if r["i3"] == i] for i in range(1, 6)]
        x5 = list(range(len(v5)))
        x3 = [len(v5) + 1 + i for i in range(len(v3))]
        ax.bar(x5, [med(v) for v in v5], color=VIOLET, alpha=0.85, edgecolor="#fcfcfb",
               linewidth=1.2, label="intron number from the 5′ end")
        ax.bar(x3, [med(v) for v in v3], color=GREEN, alpha=0.85, edgecolor="#fcfcfb",
               linewidth=1.2, label="intron number from the 3′ end")
        # n sits on top of its bar, where it can be read against the background
        for x, v in list(zip(x5, v5)) + list(zip(x3, v3)):
            ax.text(x, max(med(v), 0) + 0.03, format(len(v), ","), ha="center", va="bottom",
                    fontsize=6, color=INK2, rotation=90)
        tops.append(max(med(v) for v in v5 + v3))
        lows.append(min(med(v) for v in v5 + v3))
        ax.set_xticks(x5 + x3)
        ax.set_xticklabels([str(i) for i in range(1, 7)] + ["7+", "last"]
                           + [str(i) for i in range(2, 6)], fontsize=8)
        ax.axhline(0, color=INK2, linewidth=0.9)
        ax.set_title(iso(f))
        ax.set_xlabel("intron (labels = n)")
        mf.tidy(ax)
    axes[0].set_ylim(min(0, min(lows) * 1.2), max(tops) * 1.35)   # room for n above the bars
    axes[0].set_ylabel(MLFC)
    axes[0].legend(fontsize=7.5, ncol=2, loc="lower left", bbox_to_anchor=(0, 1.10))
    title(fig, "8. CHX response by the host intron the cassette went into, all 8 constructs",
          y=1.10)
    save(fig, "fig08_intron_index.png")


def fig09_effect_spliced_vs_unspliced(data, comp):
    """The size of the NMD effect - median(NMD expected) minus median(no NMD expected) - per
    construct and isoform, with a 95% bootstrap interval. The prediction is the same for both
    isoforms, so if both were exposed to NMD alike the two series would coincide."""
    rnd = random.Random(0)
    fig, ax = plt.subplots(figsize=(8, 3.4))
    w = 0.16
    print("  fig09: effect = median(cat 4) - median(cat 3), 95%% bootstrap CI (%d resamples)" % BOOT)
    for k, ((f, R), c, mk) in enumerate(zip(FORMS, (VIOLET, OLIVE), ("o", "s"))):
        xs, ys, lo, hi = [], [], [], []
        for i, lib in enumerate(LIBS + ["all"]):
            g = by_cat(pooled(data, lambda r: True) if lib == "all" else data[lib], f)
            e, a, b = boot_diff(g["4"], g["3"], rnd)
            xs.append(i + (k - 0.5) * 2 * w)
            ys.append(e)
            lo.append(e - a)
            hi.append(b - e)
            print("    %-9s %-5s %+.2f  [%+.2f, %+.2f]   n = %d vs %d"
                  % (f, lib, e, a, b, len(g["4"]), len(g["3"])))
        ax.errorbar(xs, ys, yerr=[lo, hi], fmt=mk, color=c, markersize=6, capsize=2.5,
                    elinewidth=1, markeredgecolor="#fcfcfb", label=iso(f), zorder=3)
    ax.axhline(0, color=INK2, linewidth=0.9, linestyle=":")
    ax.axvline(len(LIBS) - 0.5, color="#e3e3e0", linewidth=1)
    ax.set_xticks(range(len(LIBS) + 1))
    ax.set_xticklabels(LIBS + ["all 8"])
    ax.set_ylabel("difference in median log2 fold change\n(NMD expected − no NMD expected)",
                  fontsize=8)
    mf.tidy(ax)
    ax.legend(fontsize=7.5, ncol=2, loc="lower left", bbox_to_anchor=(0, 1.02))
    title(fig, "9. Size of the NMD effect per construct, spliced vs unspliced isoform", y=1.08)
    save(fig, "fig09_effect_spliced_vs_unspliced.png")


def fig10_paired_scatter(data, comp):
    """The same barcode in both isoforms: spliced against unspliced response, all constructs
    pooled. Only barcodes with a fold change in both."""
    groups = (("no NMD expected", ("2", "3"), GREEN),
              ("NMD expected", ("4",), VIOLET))
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.8), sharex=True, sharey=True)
    print("  fig10: barcodes measured in both isoforms")
    for ax, (label, cats, c) in zip(axes, groups):
        rows = pooled(data, lambda r: r["lfc"]["spliced"] is not None
                      and r["lfc"]["unspliced"] is not None and r["cat"]["spliced"] in cats)
        x = [r["lfc"]["spliced"] for r in rows]
        y = [r["lfc"]["unspliced"] for r in rows]
        ax.scatter(x, y, s=2, alpha=0.18, color=c, linewidths=0, rasterized=True)
        mx, my = med(x), med(y)
        rho = spearman(x, y)
        a, b = ols(x, y)
        ax.plot([-6, 8], [-6, 8], color=INK2, linewidth=0.8, linestyle="--", label="y = x")
        # the cloud's own colour at full opacity, thin: it should read as part of the cloud
        ax.plot([-6, 8], [a - 6 * b, a + 8 * b], color=c, linewidth=1.0, zorder=3,
                label="linear regression")
        ax.plot([], [], " ", label="Spearman ρ = %.2f" % rho)
        ax.legend(fontsize=7.5, loc="upper left")
        ax.axhline(0, color=INK2, linewidth=0.7, linestyle=":")
        ax.axvline(0, color=INK2, linewidth=0.7, linestyle=":")
        # the same quantity on both axes, so the same range and the same unit length
        ax.set_xlim(-5, 6)
        ax.set_ylim(-5, 6)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title("%s (n = %s)" % (label, format(len(rows), ",")))
        ax.set_xlabel("%s\nin the spliced isoform" % LFC)
        mf.tidy(ax, xgrid=True)
        print("    %-16s n = %6d  median spliced %+.2f  unspliced %+.2f  Spearman %.3f"
              "  OLS slope %.3f intercept %+.3f" % (label, len(rows), mx, my, rho, b, a))
    axes[0].set_ylabel("%s\nin the unspliced isoform" % LFC)
    title(fig, "10. The same barcode in both isoforms (one dot = one barcode), all 8 constructs",
          y=0.98)
    save(fig, "fig10_paired_scatter.png")


def fig11_interaction(data, comp):
    """The splicing:condition interaction of Igor's DESeq2 - how NMDi changed the barcode's
    spliced/unspliced ratio - by group. Positive = the spliced isoform gained more."""
    fig, ax = plt.subplots(figsize=(9, 3.8))
    per_lib = []
    for lib in LIBS:
        g = inter_by_cat(data[lib])
        per_lib.append([g["2"], g["3"], g["4"]])
    pos, (lo, hi) = grouped_boxes(ax, per_lib, [c for _, _, c in CATS])
    rng = hi - lo
    ax.set_ylim(lo - 0.24 * rng, hi + 0.05 * rng)
    ps = box_panel_labels(ax, per_lib, pos, lo, hi)
    print("  fig11: spliced/unspliced ratio log2FC, median per group; category 4 vs 3")
    for lib, (g2, g3, g4), p in zip(LIBS, per_lib, ps):
        print("    %-5s baseline %+.2f (n=%d)  no NMD %+.2f (n=%d)  NMD %+.2f (n=%d)   p = %.2g"
              % (lib, med(g2), len(g2), med(g3), len(g3), med(g4), len(g4),
                 p if p is not None else float("nan")))
    ax.axhline(0, color=INK2, linewidth=0.8, linestyle=":")
    ax.set_xticks(range(len(LIBS)))
    ax.set_xticklabels(LIBS)
    ax.set_ylabel("%s\nof the spliced / unspliced ratio" % LFC)
    mf.tidy(ax)
    legend_boxes(ax)
    title(fig, "11. CHX effect on the spliced / unspliced ratio, per construct", y=1.12)
    save(fig, "fig11_interaction.png")


def ols(x, y):
    """Least-squares line y = a + b x."""
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    b = (sum((u - mx) * (v - my) for u, v in zip(x, y))
         / sum((u - mx) ** 2 for u in x))
    return my - b * mx, b


def spearman(x, y):
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2.0 + 1
            i = j + 1
        return r
    rx, ry = ranks(x), ranks(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    return sxy / math.sqrt(sxx * syy)


def table_summary(data, comp):
    """summary_by_category.tsv: n, median and share > 0 per construct, measure and category, with
    category 4 vs 3 tested. The numbers a caption or a reply to Igor quotes."""
    path = os.path.join(HERE, "summary_by_category.tsv")
    with open(path, "w") as out:
        out.write("construct\tmeasure\tcategory\tn\tmedian_log2FC\tfrac_positive"
                  "\tdiff_4_minus_3\tp_4_vs_3\n")
        for lib in LIBS:
            for name, g in (("spliced isoform", by_cat(data[lib], "spliced")),
                            ("unspliced isoform", by_cat(data[lib], "unspliced")),
                            ("spliced/unspliced ratio", inter_by_cat(data[lib]))):
                p = mf.mannwhitney(g["4"], g["3"])
                diff = med(g["4"]) - med(g["3"])
                for c, label, _ in CATS:
                    v = g[c]
                    out.write("%s\t%s\t%s\t%d\t%.4f\t%.4f\t%s\t%s\n"
                              % (lib, name, label, len(v), med(v) if v else float("nan"),
                                 sum(1 for x in v if x > 0) / len(v) if v else float("nan"),
                                 "%.4f" % diff if c == "4" else "",
                                 "%.3g" % p if (c == "4" and p is not None) else ""))
    print("  summary_by_category.tsv")


def main(argv=()):
    # fig02 is left out of the default run on request (it is not in SUMMARY.md); name it to draw it
    figures = (fig01_composition, fig03_box_category, fig04_verdict_pooled,
               fig05_ptc_from_atg, fig06_positional_variables, fig07_relpos, fig08_intron_index,
               fig09_effect_spliced_vs_unspliced, fig10_paired_scatter, fig11_interaction,
               table_summary)
    if argv:
        figures = tuple(f for f in figures + (fig02_ecdf,)
                        if any(a in f.__name__ for a in argv))
        if not figures:
            print("no figure matches %s" % " ".join(argv))
            return 1
    data, comp = load()
    print("figures -> intron_nmdi/pics/")
    for f in figures:
        f(data, comp)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
