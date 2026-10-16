#!/usr/bin/env python3
"""
Stage 5 -- IntronTRIP only: fold the two insert forms into one factor.

The IntronTRIP cassette carries its own intron, so it reaches the host mRNA in TWO forms of
different length. Stage 3 is therefore run twice per reporter -- once with the spliced insert
(377 nt, an exon-exon junction at offset 181) and once with the unspliced one (503-886 nt, no
junction inside) -- and this stage joins the two BY BARCODE into the factor Igor asked for:

    spliced  unspliced   level
       0         0       No NMD expected
       1         0       NMD only if spliced
       0         1       NMD only if unspliced
       1         1       NMD either spliced or unspliced

Two things worth knowing before reading the output.

**The factor is a property of the CONSTRUCT, not of the measurement.** It says what happens to
the molecules of this barcode in which the reporter intron was excised, and what happens to those
in which it was not. It never looks at the measured splicing fraction. That is deliberate: in the
`intron_splicing` dataset that fraction IS the target, so a feature built on it would be leakage.

**On this construct the factor is expected to use only two of its four levels.** The insert carries
an in-frame stop in all three reading phases -- by design, confirmed by the reporters' author -- and
the first one lies at offset 22 or 54, inside the cryptic exon and therefore UPSTREAM of the
reporter intron at 181. Both forms consequently place the premature stop at the same position and
both keep a junction below it, so they agree. Measured on A1: 86674 of 86674 barcodes, 100.0000%.
The middle two levels are not a bug when empty; a NON-empty one would be the surprise, and this
script prints the count either way.

    python3 scripts/05_intron_factor.py [--root DIR] [--libs a,b]

in   <root>/work/nmd_<lib>_spliced.tsv, <root>/work/nmd_<lib>_unspliced.tsv
out  <root>/output/nmd_intron_<lib>.tsv
"""

import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import parse_args, split_lib

COLS = [
    "barcode",
    "nmd_expected",     # 1 / 0, empty when not analysed -- the column the ML datasets consume
    "nmd_factor",       # the four-level string, empty when not analysed
    "nmd_spliced",      # the verdict from the spliced-form run, kept so the factor is auditable
    "nmd_unspliced",    # the verdict from the unspliced-form run
    "verdict_reason",   # which branch of the decision tree answered
]

LEVELS = {
    ("0", "0"): "No NMD expected",
    ("1", "0"): "NMD only if spliced",
    ("0", "1"): "NMD only if unspliced",
    ("1", "1"): "NMD either spliced or unspliced",
}


def read_nmd(path):
    """barcode -> row, from a stage-3 output."""
    rows = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(hdr, line.rstrip("\n").split("\t")))
            rows[f["barcode"]] = f
    return rows


def main():
    root, libs = parse_args()
    work = os.path.join(root, "work")
    out_dir = os.path.join(root, "output")
    os.makedirs(out_dir, exist_ok=True)

    print("stage 5 -- IntronTRIP factor   (root %s)" % root)
    grand = Counter()

    for spec in libs:
        lib = split_lib(spec)[0]
        sp = read_nmd(os.path.join(work, "nmd_%s_spliced.tsv" % lib))
        un = read_nmd(os.path.join(work, "nmd_%s_unspliced.tsv" % lib))

        # The two runs read the same barcodes_<lib>.tsv, so a mismatch means one of them was
        # produced from a different placement file - stale output, not a biological finding.
        if set(sp) != set(un):
            sys.exit("%s: the two forms cover different barcode sets (%d vs %d) - one of the "
                     "stage-3 outputs is stale" % (lib, len(sp), len(un)))

        stat = Counter()
        disagree = 0
        path = os.path.join(out_dir, "nmd_intron_%s.tsv" % lib)
        with open(path, "w") as out:
            out.write("\t".join(COLS) + "\n")
            for bc in sp:
                a, b = sp[bc]["nmd_expected"], un[bc]["nmd_expected"]

                if a == "" and b == "":
                    row = dict(barcode=bc, nmd_expected="", nmd_factor="",
                               nmd_spliced="", nmd_unspliced="",
                               verdict_reason=sp[bc]["verdict_reason"])
                    stat["not analysed"] += 1
                else:
                    level = LEVELS[(a, b)]
                    # When the two forms agree the binary flag is simply that value. They are not
                    # expected to disagree on this construct; if they ever do, a single 0/1 is not
                    # well defined, so take 1 - the barcode's transcript population then contains
                    # NMD substrates - and make the count loud rather than silent.
                    if a != b:
                        disagree += 1
                    row = dict(barcode=bc, nmd_expected=("1" if (a == "1" or b == "1") else "0"),
                               nmd_factor=level, nmd_spliced=a, nmd_unspliced=b,
                               verdict_reason=sp[bc]["verdict_reason"])
                    stat[level] += 1
                out.write("\t".join(row[c] for c in COLS) + "\n")

        grand.update(stat)
        print("\n  %s -> output/nmd_intron_%s.tsv   (%d barcodes)" % (lib, lib, len(sp)))
        for level in ("No NMD expected", "NMD only if spliced", "NMD only if unspliced",
                      "NMD either spliced or unspliced", "not analysed"):
            if stat[level]:
                print("    %-34s %7d" % (level, stat[level]))
        analysed = len(sp) - stat["not analysed"]
        print("    forms agreed on %d of %d analysed barcodes (%.4f%%)"
              % (analysed - disagree, analysed,
                 100.0 * (analysed - disagree) / analysed if analysed else 100.0))

    print("\n  all libraries together:")
    for level in ("No NMD expected", "NMD only if spliced", "NMD only if unspliced",
                  "NMD either spliced or unspliced", "not analysed"):
        print("    %-34s %8d" % (level, grand[level]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
