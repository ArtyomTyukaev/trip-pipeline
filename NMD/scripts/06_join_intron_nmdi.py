#!/usr/bin/env python3
"""
Stage 6 -- IntronTRIP only: join the prediction with the NMD-inhibition fold changes.

The IntronTRIP counterpart of stage 4. Joins placement (stage 2) + the verdict of BOTH insert forms
(stage 3, run twice) + the factor (stage 5) + Igor's NMDi-vs-Control DESeq2 result, BY BARCODE, into
one table per construct. Every barcode of the genuine BED appears; nothing is filtered.

What differs from stage 4, and why:

**Three fold changes instead of one.** The CHX tables of single_combined carry one contrast on total
RNA. The IntronTRIP ones are split by whether the reporter intron was excised, so they carry
`Spliced_*` (NMDi/Control on spliced molecules), `Unspliced_*` (on molecules that kept it) and
`SplicingInteraction_*` (the splicing:condition interaction - did NMDi shift the spliced/unspliced
ratio). Each gets its own expressed/not-expressed call, so there are two categories per barcode,
`category_spliced` and `category_unspliced`, defined exactly as stage 4 defines `category`.

**The per-form verdict details are kept side by side.** The two forms always agree on the verdict
on this construct (stage 5), but not on where the PTC sits: in one reading phase the insert's first
stop lies downstream of the reporter intron, so `ptc_tx_pos` and the distance to the last junction
differ between forms. Those columns are written twice, suffixed `_spliced` / `_unspliced`, so the
spliced isoform's fold change can be held against the spliced geometry and the unspliced one
against its own.

**Igor's annotation is checked against our placement before it is used.** His tables take Chr/Pos,
ChromHMM and host-gene expression from his annotated tables, which for five constructs are an older
run than the genuine BED his counts were joined to. `annot_match` records, per barcode, whether his
coordinates agree with ours: `same`, `shifted` (same chromosome and strand, within 100 nt - a
few-nt TTAA convention difference), `other_locus`, or `absent` (not in his table). ChromHMM and
host expression are carried only where the locus is the same or shifted; elsewhere they describe a
different integration and are left empty.

As in stage 4, only the DESeq2 result and those two annotation values are reused. Gene identity,
exon/intron placement and the verdict are ours.

    python3 scripts/06_join_intron_nmdi.py [--root DIR] [--libs a,b]

in   <root>/work/barcodes_<lib>.tsv
     <root>/work/nmd_<lib>_spliced.tsv, <root>/work/nmd_<lib>_unspliced.tsv
     <root>/output/nmd_intron_<lib>.tsv
     <root>/logFC/CHX_<lib>.txt
out  <root>/output/nmd_nmdi_<lib>.tsv
"""

import csv
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import parse_args, split_lib

FORMS = ("spliced", "unspliced")
MEASURES = ("Spliced", "Unspliced", "SplicingInteraction")
STATS = ("baseMean", "log2FoldChange", "lfcSE", "padj")
SHIFT_MAX = 100          # nt; Igor's Pos vs ours on the same chromosome and strand

FROM_BC = ["chrom", "pos", "strand", "score",
           "in_gene", "excluded_reason", "tx_id", "n_tx_overlap",
           "loc", "feat_idx", "n_exons", "cds_last_exon", "tx_pos", "region"]
# form-independent by construction: where the insert lands sets its reading phase
FROM_NMD_SHARED = ["insert_frame"]
FROM_NMD_FORM = ["insert_len", "bc_has_stop", "ptc_source", "ptc_tx_pos", "last_junction_pos",
                 "n_junc_downstream", "dist_to_last_junction", "dist_gt_50", "selenocysteine",
                 "note"]
FROM_FACTOR = ["nmd_expected", "nmd_factor", "verdict_reason"]
MEASURE_COLS = ["%s_%s" % (r, s) for r in MEASURES for s in STATS]

COLS = (["barcode"] + FROM_BC[:4] + ["tx_id", "gene_symbol"]
        + ["in_gene", "excluded_reason", "n_tx_overlap", "loc", "feat_idx", "n_exons",
           "cds_last_exon", "tx_pos", "region"]
        + FROM_FACTOR + FROM_NMD_SHARED
        + ["%s_%s" % (c, f) for f in FORMS for c in FROM_NMD_FORM]
        + MEASURE_COLS
        + ["chromhmm", "host_exp_log10", "annot_match", "annot_pos_offset"]
        + ["category_spliced", "category_unspliced"])


def load_tsv(path, key="barcode", keep=None):
    """Keyed by barcode, and it refuses to collapse a duplicate key - same guard as stage 4.

    `keep` restricts the stored columns, and every one of them must be in the header. Igor's tables
    are 84 columns wide and this stage needs a sixth of them."""
    out = {}
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if keep is not None:
            missing = [c for c in keep if c not in rd.fieldnames]
            if missing:
                sys.exit("%s lacks %s" % (os.path.basename(path), ", ".join(missing)))
        for r in rd:
            if r[key] in out:
                sys.exit("%s: %s appears twice -- refusing to join" % (os.path.basename(path),
                                                                       r[key]))
            out[r[key]] = r if keep is None else {k: r[k] for k in keep}
    return out


def na(v):
    return v is None or v in ("", "NA", "NaN", "nan", "NULL")


def category(lfc, brow, nmd_expected):
    """Stage 4's definition, applied to one isoform's fold change."""
    if na(lfc):
        return 1
    if brow["in_gene"] != "1":
        return 2
    if na(nmd_expected):
        return ""                         # in a gene but no verdict -- should not happen
    return 4 if nmd_expected == "1" else 3


def annot_match(brow, crow):
    if crow is None:
        return "absent", ""
    if crow["Chr"] != brow["chrom"] or crow["Strand"] != brow["strand"] or na(crow["Pos"]):
        return "other_locus", ""
    off = int(float(crow["Pos"])) - int(brow["pos"])
    if off == 0:
        return "same", "0"
    if abs(off) <= SHIFT_MAX:
        return "shifted", str(off)
    return "other_locus", str(off)


def main():
    root, libs = parse_args()
    work = os.path.join(root, "work")
    outdir = os.path.join(root, "output")
    os.makedirs(outdir, exist_ok=True)

    print("stage 6 -- IntronTRIP prediction x NMDi fold changes   (root %s)" % root)
    for spec in libs:
        lib = split_lib(spec)[0]
        bc = load_tsv(os.path.join(work, "barcodes_%s.tsv" % lib))
        nmd = {f: load_tsv(os.path.join(work, "nmd_%s_%s.tsv" % (lib, f))) for f in FORMS}
        fac = load_tsv(os.path.join(outdir, "nmd_intron_%s.tsv" % lib))
        for name, d in (("spliced", nmd["spliced"]), ("unspliced", nmd["unspliced"]),
                        ("factor", fac)):
            if set(d) != set(bc):
                sys.exit("%s: the %s table covers a different barcode set than barcodes_%s.tsv "
                         "- stale stage output" % (lib, name, lib))

        chx = load_tsv(os.path.join(root, "logFC", "CHX_%s.txt" % lib), key="Barcode",
                       keep=["Barcode", "Chr", "Pos", "Strand", "ChromHMM",
                             "MANE_all_exp.norm.log10", "MANE_all", "MANE_all_symbol"]
                            + MEASURE_COLS)
        symbol = {}
        for r in chx.values():
            if not na(r.get("MANE_all")) and not na(r.get("MANE_all_symbol")):
                symbol.setdefault(r["MANE_all"], r["MANE_all_symbol"])

        # A barcode Igor measured that is not in the BED would be a measurement we cannot place.
        # It must not happen when the BED is the one his counts were joined to; count it loudly.
        measured_unplaced = sum(1 for b, r in chx.items() if b not in bc
                                and any(not na(r["%s_log2FoldChange" % x]) for x in MEASURES))
        annot_only = sum(1 for b in chx if b not in bc) - measured_unplaced

        stat = Counter()
        dst = os.path.join(outdir, "nmd_nmdi_%s.tsv" % lib)
        with open(dst, "w") as out:
            out.write("\t".join(COLS) + "\n")
            for b, brow in bc.items():
                crow = chx.get(b)
                frow = fac[b]
                row = {"barcode": b, "gene_symbol": symbol.get(brow["tx_id"], "")}
                row.update({k: brow[k] for k in FROM_BC})
                row.update({k: frow[k] for k in FROM_FACTOR})
                # The forms are expected to share the phase; if not, the table would be lying.
                if nmd["spliced"][b]["insert_frame"] != nmd["unspliced"][b]["insert_frame"]:
                    sys.exit("%s %s: insert_frame differs between forms" % (lib, b))
                row["insert_frame"] = nmd["spliced"][b]["insert_frame"]
                for f in FORMS:
                    for c in FROM_NMD_FORM:
                        row["%s_%s" % (c, f)] = nmd[f][b][c]
                for c in MEASURE_COLS:
                    v = crow[c] if crow else ""
                    row[c] = "" if na(v) else v

                match, offset = annot_match(brow, crow)
                row["annot_match"], row["annot_pos_offset"] = match, offset
                stat["annot %s" % match] += 1
                usable = match in ("same", "shifted")
                ch = crow["ChromHMM"] if (crow and usable) else ""
                he = crow["MANE_all_exp.norm.log10"] if (crow and usable) else ""
                row["chromhmm"] = "" if na(ch) else ch
                row["host_exp_log10"] = "" if na(he) else he

                for f, r in (("spliced", "Spliced"), ("unspliced", "Unspliced")):
                    cat = category(row["%s_log2FoldChange" % r], brow, frow["nmd_expected"])
                    row["category_%s" % f] = cat
                    stat["%s category %s" % (f, cat)] += 1
                out.write("\t".join(str(row[c]) for c in COLS) + "\n")

        print("\n  %s -> output/nmd_nmdi_%s.tsv   (%d rows, %d columns)"
              % (lib, lib, len(bc), len(COLS)))
        print("    CHX_%s.txt: %d rows; %d not in the BED and unmeasured (annotation only), "
              "%d not in the BED but MEASURED" % (lib, len(chx), annot_only, measured_unplaced))
        for k in sorted(stat):
            print("    %-28s %7d" % (k, stat[k]))
        if measured_unplaced:
            print("    [!] %d measured barcodes have no placement - wrong BED for this table?"
                  % measured_unplaced)
    return 0


if __name__ == "__main__":
    sys.exit(main())
