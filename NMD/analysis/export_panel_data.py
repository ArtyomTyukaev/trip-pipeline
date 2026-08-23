#!/usr/bin/env python3
"""Dump the numbers behind one figure panel to TSV.

The classification rules (what counts as a PTC, which barcodes are expressed) live in
make_figures.py and are imported from there, so a plotting script written in another language can
draw the panel without re-implementing any of them.

    module load matplotlib/3.7.0-gfbf-2022b
    python3 analysis/export_panel_data.py

out  analysis/data/fig04_ptc_yesno.tsv   one row per expressed barcode: library, group, log2FC
"""

import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from make_figures import DISPLAY, LIBS, expressed, load   # noqa: E402

# The three groups of the fig1G panel, in bottom-to-top order.
#
# The PTC group is `nmd_ptc`, so it holds only stops the last-exon rule and the 50-nt rule both let
# through - the ones actually expected to trigger NMD.
#
# The control group deliberately stays on `not has_ptc`, the WIDER flag, and not on `not nmd_ptc`.
# The difference matters: a barcode whose PTC fails the 50-nt rule has a PTC, so `not nmd_ptc`
# would file it among the barcodes that make none and quietly contaminate the control with the very
# thing being tested. Those barcodes are simply dropped from the panel instead - 45 of 2,979 in
# 5ss_noP and 63 of 3,847 in CPA_noP, which is why the medians barely move.
#
# fig04_ptc() in make_figures.py still draws the wider `has_ptc` version; the two panels answer
# slightly different questions and both are worth keeping.
GROUPS = [("outside_gene", "Outside a gene", lambda r: r["category"] == "2"),
          ("no_ptc", "No PTC created", lambda r: r["in_gene"] == "1" and not r["has_ptc"]),
          ("ptc", "PTC expected to trigger NMD", lambda r: r["nmd_ptc"])]


def main():
    out_dir = os.path.join(HERE, "data")
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "fig04_ptc_yesno.tsv")
    data = load()
    n = 0
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["library", "library_label", "group", "group_label", "lfc"])
        for lib, _ in LIBS:
            rows = expressed(data[lib])
            for key, label, sel in GROUPS:
                for r in rows:
                    if sel(r):
                        w.writerow([lib, DISPLAY[lib], key, label, "%.6g" % r["lfc"]])
                        n += 1
    print("  analysis/data/fig04_ptc_yesno.tsv  (%d rows)" % n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
