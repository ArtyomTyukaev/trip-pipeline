#!/usr/bin/env python3
"""
Build the toy dataset: an invented miniature genome small enough to work out by hand.

Two roots are produced, identical except for the insert:

    root_A/   insert 12 nt  -- a multiple of 3, so a read-through does NOT frameshift the host
    root_B/   insert 13 nt  -- not a multiple of 3, so a read-through DOES frameshift the host

Each root is laid out exactly like NMD/ itself, so the production scripts run over it unchanged
with --root. That is the point of the test: it exercises the real code, not a copy of it.

The insert template mimics the real 400-nt one in the only way that matters -- it carries a stop in
two of the three reading frames and is clean in the third:

    G T A A [N N N N] C T A G (C)
    0 1 2 3  4 5 6 7  8 9 ...

    frame 0 -> codons at 0,3,6,9   -> ...,"TAG" at offset 9      stop inside the insert
    frame 1 -> codons at 1,4,7,10  -> "TAA" at offset 1          stop inside the insert
    frame 2 -> codons at 2,5,8,11  -> no stop unless the barcode supplies one at offset 5

Four toy genes, each present to force one thing to be checked:

    TOY_PLUS   7 exons, + strand  -- 5'UTR introns (one of them exactly at the ATG), CDS introns in
                                    all three frames, a 3'UTR intron, exonic cases in every region
    TOY_MINUS  3 exons, - strand  -- intron numbering in TRANSCRIPT order, which is the reverse of
                                    the genomic order; this is the inversion the annotation gets wrong
    TOY_SINGLE 1 exon,  + strand  -- intronless gene: its only exon is the last exon
    TOY_SEAM   3 exons, + strand  -- host ends ...TA right before the junction, so host "TA" plus the
                                    insert's first base "G" makes a TAG straddling the seam

    python3 toy_test/make_toy.py
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
STOPS = ("TAA", "TAG", "TGA")
INSERTS = {"A": "GTAA" + "NNNN" + "CTAG", "B": "GTAA" + "NNNN" + "CTAGC"}

# ---------------------------------------------------------------------------- the toy genes
# exon_lens are in TRANSCRIPT order. tx_seq is the mature mRNA, 5'->3'.
GENES = [
    dict(
        tx_id="NM_900001.1", symbol="TOY_PLUS", chrom="chr1", strand="+", gstart=1000,
        exon_lens=[6, 3, 13, 14, 62, 8, 7],
        tx_seq="CACACACAC" + "ATG" + "GCT" * 29 + "TAA" + "CACACACACAC",
        cds=(9, 102),
    ),
    dict(
        tx_id="NM_900002.1", symbol="TOY_MINUS", chrom="chr1", strand="-", gstart=2000,
        exon_lens=[9, 60, 12],
        tx_seq="ATG" + "GCT" * 22 + "TAA" + "CACACACAC",
        cds=(0, 72),
    ),
    dict(
        tx_id="NM_900003.1", symbol="TOY_SINGLE", chrom="chr1", strand="+", gstart=3000,
        exon_lens=[21],
        tx_seq="CAC" + "ATGGCTGCTTAA" + "CACACA",
        cds=(3, 15),
    ),
    dict(
        # tx9,tx10 are "TA"; tx11 is the first base of exon 2 and is "C", so the unmodified
        # transcript reads TAC there. Insert anything starting with G and it becomes TAG.
        tx_id="NM_900004.1", symbol="TOY_SEAM", chrom="chr1", strand="+", gstart=4000,
        exon_lens=[11, 9, 12],
        tx_seq="ATGGCTGCTTA" + "CGAGCAGCA" + "CTAACACACACA",
        cds=(0, 24),
    ),
]

# ---------------------------------------------------------------------------- the toy barcodes
# (barcode, chrom, pos, strand, what it is there to test)
BARCODES = [
    ("ACCC", "chr1", 1010, "+", "TOY_PLUS intron 1 -- 5'UTR intron, strictly above the ATG"),
    ("ACCG", "chr1", 1024, "+", "TOY_PLUS intron 2 -- junction lands exactly ON the ATG"),
    ("CCCC", "chr1", 1046, "+", "TOY_PLUS intron 3 -- CDS, frame 2, insert read through"),
    ("CTAA", "chr1", 1048, "+", "TOY_PLUS intron 3 -- CDS, frame 2, the BARCODE supplies the stop"),
    ("CCCG", "chr1", 1070, "+", "TOY_PLUS intron 4 -- CDS, frame 0, stop inside the insert"),
    ("CCTC", "chr1", 1142, "+", "TOY_PLUS intron 5 -- CDS, frame 1, stop inside the insert"),
    ("CCTG", "chr1", 1160, "+", "TOY_PLUS intron 6 -- 3'UTR intron"),
    ("ACGC", "chr1", 1002, "+", "TOY_PLUS exon 1 -- exonic in the 5'UTR"),
    ("ACGG", "chr1", 1056, "+", "TOY_PLUS exon 4 -- exonic, internal coding exon"),
    ("ACTC", "chr1", 1150, "+", "TOY_PLUS exon 6 -- exonic in the CDS, not the last exon"),
    ("ACTG", "chr1", 1168, "+", "TOY_PLUS exon 7 -- exonic in the 3'UTR, and the last exon"),
    ("GCCC", "chr1", 2016, "-", "TOY_MINUS genomic intron 1 = TRANSCRIPT intron 2"),
    ("GCCG", "chr1", 2086, "-", "TOY_MINUS genomic intron 2 = TRANSCRIPT intron 1"),
    ("GCTC", "chr1", 2010, "-", "TOY_MINUS last exon, in the CDS"),
    ("TCCC", "chr1", 3006, "+", "TOY_SINGLE -- intronless, its only exon is the last exon"),
    ("TCCG", "chr1", 4016, "+", "TOY_SEAM intron 1 -- the PTC straddles the host|insert seam"),
    ("TCGC", "chr1", 4041, "+", "TOY_SEAM last exon, in the CDS"),
    ("AGGG", "chr1", 5000, "+", "intergenic"),
    ("AGGC", "chr1", 1046, "-", "antisense inside TOY_PLUS"),
    ("TCGG", "chr1", 1070, "+", "same site as CCCG but NA fold change -- category 1 must win"),
    ("TTTT", "chr1", 5100, "+", "absent from logFC altogether -- must not vanish"),
]

NO_FOLDCHANGE = {"TCGG"}      # present in logFC but with NA
ABSENT = {"TTTT"}             # not in logFC at all


def genomic_exons(g, intron=10):
    """Exon intervals in ascending GENOMIC order."""
    glens = g["exon_lens"] if g["strand"] == "+" else list(reversed(g["exon_lens"]))
    exons, p = [], g["gstart"]
    for L in glens:
        exons.append((p, p + L))
        p += L + intron
    return exons


def tx_to_genomic(exons, strand, total, txpos):
    fwd = txpos if strand == "+" else total - 1 - txpos
    off = 0
    for a, b in exons:
        if off <= fwd < off + (b - a):
            return a + (fwd - off)
        off += b - a
    raise AssertionError("tx position %d outside the transcript" % txpos)


def check(g):
    """Refuse to emit a toy gene that is not a legal transcript."""
    total = sum(g["exon_lens"])
    assert len(g["tx_seq"]) == total, \
        "%s: sequence %d nt but exons sum to %d" % (g["symbol"], len(g["tx_seq"]), total)
    cs, ce = g["cds"]
    cds = g["tx_seq"][cs:ce]
    assert len(cds) % 3 == 0, "%s: CDS %d nt, not divisible by 3" % (g["symbol"], len(cds))
    assert cds[:3] == "ATG", "%s: starts with %s" % (g["symbol"], cds[:3])
    assert cds[-3:] in STOPS, "%s: ends with %s" % (g["symbol"], cds[-3:])
    internal = [i for i in range(0, len(cds) - 3, 3) if cds[i:i + 3] in STOPS]
    assert not internal, "%s: in-frame stop inside the CDS at %s" % (g["symbol"], internal)
    return total


def write_root(tag):
    root = os.path.join(HERE, "root_%s" % tag)
    gb = os.path.join(root, "genuine_barcodes")
    lf = os.path.join(root, "logFC")
    for d in (gb, lf):
        os.makedirs(d, exist_ok=True)

    with open(os.path.join(gb, "MANE.GRCh38.v1.0.refseq_genomic.bed"), "w") as bed, \
         open(os.path.join(gb, "MANE.GRCh38.v1.0.refseq_rna.fna"), "w") as fna:
        for g in GENES:
            total = check(g)
            exons = genomic_exons(g)
            cs, ce = g["cds"]
            g1 = tx_to_genomic(exons, g["strand"], total, cs)
            g2 = tx_to_genomic(exons, g["strand"], total, ce - 1)
            thick_start, thick_end = min(g1, g2), max(g1, g2) + 1
            bed.write("\t".join([
                g["chrom"], str(exons[0][0]), str(exons[-1][1]), g["tx_id"], "0", g["strand"],
                str(thick_start), str(thick_end), "0", str(len(exons)),
                ",".join(str(b - a) for a, b in exons) + ",",
                ",".join(str(a - exons[0][0]) for a, b in exons) + ",",
            ]) + "\n")
            fna.write(">%s %s toy\n%s\n" % (g["tx_id"], g["symbol"], g["tx_seq"]))

    with open(os.path.join(gb, "exon_and_5arm_sequence.txt"), "w") as fh:
        fh.write(INSERTS[tag] + "\n")

    with open(os.path.join(gb, "5ss_noP.txt.retained_barcodes_genuine.bed"), "w") as fh:
        for bc, chrom, pos, strand, _ in BARCODES:
            fh.write("\t".join([chrom, str(pos), str(pos + 1), bc, "1", strand]) + "\n")

    # a toy logFC: only the columns stage 4 actually reads
    cols = ["Barcode", "MANE_all", "MANE_all_symbol",
            "CHXbaseMean", "CHXlog2FoldChange", "CHXlfcSE", "CHXpadj"]
    sym = {g["tx_id"]: g["symbol"] for g in GENES}
    with open(os.path.join(lf, "CHX_5ss.txt"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for i, (bc, chrom, pos, strand, _) in enumerate(BARCODES):
            if bc in ABSENT:
                continue
            host = ""
            for g in GENES:
                e = genomic_exons(g)
                if g["strand"] == strand and e[0][0] <= pos < e[-1][1]:
                    host = g["tx_id"]
            lfc = "NA" if bc in NO_FOLDCHANGE else "%.3f" % (0.5 - 0.1 * i)
            fh.write("\t".join([bc, host or "NA", sym.get(host, "NA"),
                                "100.0", lfc, "0.20", "0.05"]) + "\n")
    return root


def summary():
    """The mechanical facts of the toy genes, printed so the hand calculation can be set up."""
    print("toy genes")
    for g in GENES:
        total = sum(g["exon_lens"])
        exons = genomic_exons(g)
        order = g["exon_lens"]
        junc, off = [], 0
        for L in order[:-1]:
            off += L
            junc.append(off)
        cs, ce = g["cds"]
        print("  %-11s %s %s  tx_len %d  exons %s" % (g["symbol"], g["chrom"], g["strand"],
                                                      total, g["exon_lens"]))
        print("    genomic exons  %s" % ", ".join("%d-%d" % e for e in exons))
        print("    junctions (tx) %s" % junc)
        print("    cds (tx)       %d..%d   3'UTR starts at %d" % (cs, ce, ce))
    print("\ninserts")
    for tag, s in sorted(INSERTS.items()):
        print("  %s: %-22s len %d, len %% 3 = %d, slot at %d-%d"
              % (tag, s, len(s), len(s) % 3, s.find("N"), s.rfind("N") + 1))


def main():
    for tag in sorted(INSERTS):
        root = write_root(tag)
        print("wrote %s" % root)
    print()
    summary()
    return 0


if __name__ == "__main__":
    sys.exit(main())
