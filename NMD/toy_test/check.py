#!/usr/bin/env python3
"""
Compare the pipeline output on the toy data against the hand calculation.

One FAIL line per differing field, so a failure names the barcode, the column, what was expected and
what came out. Also refuses to run on a malformed expectation file -- a row with the wrong number of
fields is a bug in the test, not in the pipeline, and it must not masquerade as one.

In the expectation files a dot means "this field must be empty".

    python3 toy_test/check.py [A|B]
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def load_expected(path):
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if not l.startswith("#") and l.strip()]
    hdr = lines[0].split("\t")
    rows = {}
    for ln, line in enumerate(lines[1:], start=2):
        f = line.split("\t")
        if len(f) != len(hdr):
            sys.exit("%s line %d: %d fields, header has %d -- fix the expectation file"
                     % (os.path.basename(path), ln, len(f), len(hdr)))
        rows[f[0]] = {k: ("" if v == "." else v) for k, v in zip(hdr, f)}
    return hdr, rows


def load_output(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        rows = {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows[f[0]] = dict(zip(hdr, f))
    return rows


def run(tag):
    hdr, exp = load_expected(os.path.join(HERE, "expected_%s.tsv" % tag))
    got = load_output(os.path.join(HERE, "root_%s" % tag, "output", "nmd_table_5ss_noP.tsv"))
    checked = [c for c in hdr if c != "barcode"]

    fails = []
    for bc, want_row in exp.items():
        if bc not in got:
            fails.append("%s: missing from the output entirely" % bc)
            continue
        for col in checked:
            want, have = want_row[col], got[bc].get(col, "<no such column>")
            if want != have:
                fails.append("%-6s %-22s expected %-18r got %r" % (bc, col, want, have))
    for bc in got:
        if bc not in exp:
            fails.append("%s: in the output but not in the hand calculation" % bc)

    print("variant %s: %d barcodes x %d fields = %d comparisons"
          % (tag, len(exp), len(checked), len(exp) * len(checked)))
    for f in fails:
        print("  FAIL  %s" % f)
    print("  %s" % ("PASS -- output matches the hand calculation exactly" if not fails
                    else "%d FAILURES" % len(fails)))
    return len(fails)


def main():
    tags = sys.argv[1:] or ["A", "B"]
    return 1 if sum(run(t) for t in tags) else 0


if __name__ == "__main__":
    sys.exit(main())
