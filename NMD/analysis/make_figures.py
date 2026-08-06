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
import os
import random
import sys
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
DOTS_UP_TO = 500        # boxplots smaller than this also show every barcode as a dot

plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.edgecolor": INK2, "axes.linewidth": 0.8,
    "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "grid.color": GRID, "grid.linewidth": 0.6,
    "legend.frameon": False, "figure.dpi": 150, "savefig.bbox": "tight",
})


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
                rows.append(r)
        data[lib] = rows
    return data


def expressed(rows):
    return [r for r in rows if r["lfc"] is not None]


def save(fig, name):
    fig.savefig(os.path.join(PICS, name), facecolor=fig.get_facecolor())
    plt.close(fig)
    print("  pics/%s" % name)


def tidy(ax, xgrid=False, ygrid=True):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    ax.set_axisbelow(True)
    if ygrid:
        ax.yaxis.grid(True)
    if xgrid:
        ax.xaxis.grid(True)


def box_groups(ax, groups, colors, xlim=(-3.5, 6.5), ylabels=True, min_n=5,
               xlabel="log2 fold change (CHX / EtOH)"):
    """Horizontal boxes over a jittered strip of the individual barcodes.

    For a group of fewer than DOTS_UP_TO barcodes every barcode is also drawn as one dot behind the
    box, so the reader sees the actual sample rather than five summary numbers; the box is kept
    translucent for the dots to show through. Bigger groups get the box alone. Whiskers
    are the 5th/95th percentiles -- with a handful of barcodes at log2FC 12 the full range squashes
    every box into a sliver."""
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
        ax.set_yticklabels(["%s  (n=%d)" % (g[0], len(g[1])) for g in groups])
    else:
        ax.set_yticklabels(["n=%d" % len(g[1]) for g in groups], fontsize=7.5, color=INK2)
    ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel)
    tidy(ax, xgrid=True, ygrid=False)


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
    ax.plot(cx, cm, color=color, linewidth=2, marker="o", markersize=4,
            markeredgecolor="#fcfcfb", markeredgewidth=0.8, label=label)
    return ns


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
    """Where integrations land along a gene, with every gene rescaled to 0..1 so genes of different
    length and exon count can be pooled. Top row: along the gene as it sits in the genome. Bottom
    row: along the mature mRNA -- for an intronic hit that position is the junction the cassette
    splices into, since the intron itself is not in the mRNA. A flat line at 1.0 is what a uniform
    landscape would give."""
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 5.2), sharex=True, sharey="row")
    bins = [i / 20 for i in range(21)]
    rows_spec = [("relgene", "along the gene (genomic span)"),
                 ("relpos", "along the mature mRNA")]
    for row, (key, title) in enumerate(rows_spec):
        for col, (lib, _) in zip(range(2), LIBS):
            ax = axes[row][col]
            base = [r for r in data[lib] if r["in_gene"] == "1" and r.get(key) is not None]
            for label, sel, col_ in (("all integrations", lambda r: True, GREY),
                                     ("expressed only", lambda r: r["lfc"] is not None, ORANGE)):
                v = [r[key] for r in base if sel(r)]
                ax.hist(v, bins=bins, density=True, histtype="step", linewidth=2, color=col_,
                        label="%s (n=%d)" % (label, len(v)))
            ax.axhline(1.0, color=INK2, linewidth=0.9, linestyle=":")
            ax.set_title("%s — %s" % (lib, title), fontsize=9)
            tidy(ax)
            if row == 1:
                ax.set_xlabel("relative position: 0 = 5' end, 1 = 3' end")
            if col == 0:
                ax.set_ylabel("density of integrations")
            if row == 0 and col == 0:
                ax.legend(fontsize=7.5, loc="upper right")
    save(fig, "fig01c_position_landscape.png")


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
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 2.6))
    for ax, (lib, _) in zip(axes, LIBS):
        rows = expressed(data[lib])
        groups = [("outside a gene (baseline)", [r["lfc"] for r in rows if r["category"] == "2"]),
                  ("no PTC created", [r["lfc"] for r in rows
                                        if r["in_gene"] == "1" and not r["has_ptc"]]),
                  ("PTC created", [r["lfc"] for r in rows if r["has_ptc"]])]
        box_groups(ax, groups, [GREY, BLUE, ORANGE], ylabels=(ax is axes[0]))
        ax.set_title(lib)
    save(fig, "fig04_ptc_yesno.png")


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
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.0))
    names = [("stop inside the insert\n(insert_5p)", "insert_5p"),
             ("stop from the barcode", "barcode"),
             ("stop on the host|insert\nseam", "host_insert_seam"),
             ("read through, stop\nfound in the host", "host_downstream")]
    for ax, (lib, _) in zip(axes, LIBS):
        rows = expressed(data[lib])
        groups, cols = [], []
        for label, key in names:
            v = [r["lfc"] for r in rows if r["ptc_source"] == key]
            if len(v) >= 5:
                groups.append((label, v))
                cols.append(ORANGE)
        box_groups(ax, groups, cols, ylabels=(ax is axes[0]))
        ax.set_title(lib)
        ax.tick_params(axis="y", labelsize=7.5)
    save(fig, "fig08_ptc_source.png")


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
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.2), sharey=True)
    for ax, (lib, _) in zip(axes, LIBS):
        rows = [r for r in expressed(data[lib]) if r["relpos"] is not None]
        for label, sel, col in (("PTC created", lambda r: r["has_ptc"], ORANGE),
                                ("no PTC created (UTR)",
                                 lambda r: r["in_gene"] == "1" and not r["has_ptc"], BLUE)):
            s = [r for r in rows if sel(r)]
            if len(s) > 60:
                binned(ax, [r["relpos"] for r in s], [r["lfc"] for r in s], 10, col,
                       label="%s (n=%d)" % (label, len(s)), band=False)
        ax.axhline(0, color=INK2, linewidth=0.9, linestyle=":")
        ax.set_title(lib)
        ax.set_xlabel("insertion point: 0 = 5' end, 1 = 3' end")
        tidy(ax)
        ax.legend(fontsize=7.5, loc="upper right")
    axes[0].set_ylabel("median log2FC")
    save(fig, "fig11_relpos.png")


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


def main():
    os.makedirs(PICS, exist_ok=True)
    data = load()
    print("figures -> analysis/pics/")
    for f in (fig01_composition, fig01b_loc_counts, fig01c_position_landscape,
              fig02_ecdf, fig03_verdict, fig04_ptc, fig05_dist_hist,
              fig07_matched, fig08_ptc_source, fig09_frame, fig10_intron_index,
              fig11_relpos, fig11b_ptc_from_atg, fig11c_positional_variables,
              fig11d_position_control, fig12_region_loc, fig13_ma, fig14_host_expr, fig15_chromatin,
              fig16_reproducibility):
        f(data)
    return 0


if __name__ == "__main__":
    sys.exit(main())
