# toy_test — verifying the pipeline on invented data

The idea: build a tiny invented genome, small enough that every answer can be worked out by hand,
run **the production scripts themselves** over it (not a copy), and compare against that hand
calculation. If they agree, the coordinate, reading-frame and junction logic is right.

```bash
bash toy_test/run_toy.sh
```

The scripts are called with `--root toy_test/root_A`, so they read the toy files instead of the real
ones. Nothing touches the real data.

Two runs, differing only in the insert:

| root | insert | length | why |
|---|---|---|---|
| `root_A` | `GTAANNNNCTAG` | 12 nt, **a multiple of 3** | reading the insert through does NOT frameshift the host |
| `root_B` | `GTAANNNNCTAGC` | 13 nt, **not a multiple of 3** | reading it through DOES frameshift the host |

**The real insert behaves like B**: it is 400 nt and `400 mod 3 = 1`. `root_A` exercises the branch
that never fires on real data.

Running both also proves nothing about the insert is hard-coded: its length and the barcode slot are
read from the file, the slot being located as the run of `N`s.

---

## The toy insert

The real insert is 400 nt. The toy one mimics it in the single property that matters: **a stop
exists in two of the three reading frames and is absent in the third.**

```
    G  T  A  A [ N  N  N  N ] C  T  A  G (C)
    0  1  2  3   4  5  6  7    8  9  10 11 12
                 └─ barcode slot ─┘

frame 0 -> codons at 0,3,6,9  -> GTA, A+bc, bc+C, TAG  -> stop at offset 9
frame 1 -> codons at 1,4,7,10 -> TAA                   -> stop at offset 1
frame 2 -> codons at 2,5,8,11 -> AAb, bbb, CTA, ...    -> no stop unless the barcode supplies one
```

Exactly as in the real construct: frame 1 is clean (no stop across all 400 nt), frames 0 and 2 stop
inside the insert. And the 4-nt barcode occupies whole codons, so it can contribute a stop itself.

---

## The toy genes

Four genes, each present to force one thing to be checked. Their actual BED12 lines:

```
chr1  1000 1173  NM_900001.1  0  +  1029 1152  0  7  6,3,13,14,62,8,7,  0,16,29,52,76,148,166,
chr1  2000 2101  NM_900002.1  0  -  2009 2101  0  3  12,60,9,           0,22,92,
chr1  3000 3021  NM_900003.1  0  +  3003 3015  0  1  21,                0,
chr1  4000 4052  NM_900004.1  0  +  4000 4044  0  3  11,9,12,           0,21,40,
```

### TOY_PLUS (`NM_900001.1`), plus strand, 7 exons

```
genome: 1000──1006  1016─1019  1029─────1042  1052──────1066  1076────────────1138  1148───1156  1166──1173
        └ E1, 6 ┘   └ E2, 3 ┘  └ E3, 13 ┘     └ E4, 14 ┘      └ E5, 62 ┘           └ E6, 8 ┘    └ E7, 7 ┘
            i1        i2          i3             i4               i5                   i6

transcript (113 nt):
  0    6  9              22            36                            98      106    113
  ├─E1─┼E2┼──────E3──────┼─────E4──────┼──────────E5────────────────┼──E6───┼──E7──┤
  └── 5'UTR ──┤
              ATG (nt 9)                                       stop TAA (nt 99-101)
                                                                            └ 3'UTR from nt 102 ┘

junctions (transcript coordinates): 6, 9, 22, 36, 98, 106
CDS: 9 to 102 (exclusive); stop codon TAA at nt 99-101, i.e. in EXON 6
cds_last_exon = 6 while n_exons = 7 -- exon 7 is a pure 3'UTR exon behind a 3'UTR intron
```

Sequence, which is exactly the FASTA record:
`CACACACAC` + `ATG` + `GCT`×29 + `TAA` + `CACACACACAC`

What it is for: **two 5'UTR introns** (one whose junction lands exactly on the ATG — the boundary
case), **three CDS introns in all three frames**, a **3'UTR intron**, and exonic sites in every
region. Above all, that 3'UTR intron makes `cds_last_exon = 6 ≠ n_exons = 7`, so this gene is what
checks that "last exon" in the rule means the exon holding the stop codon, not the transcript's last.

### TOY_MINUS (`NM_900002.1`), minus strand, 3 exons

Here intron order in the transcript is the **reverse** of the genomic order — the very inversion the
annotation's `MANE_intron_ind` gets wrong.

```
genome (left to right):  2000──2012      2022────────2082      2092──2101
                         genomic E1      genomic E2            genomic E3
                         = transcript E3 = transcript E2       = transcript E1
                               genomic i1         genomic i2
                               = transcript i2    = transcript i1

transcript (81 nt), running right to left along the genome:
  0                     9                                69          81
  ├────── E1, 9 ────────┼───────── E2, 60 ───────────────┼── E3, 12 ─┤
  ATG (nt 0)                                       stop TAA (69-71), 3'UTR from 72

junctions: 9, 69      CDS: 0 to 72
```

A barcode in genomic intron 1 (2012–2022) is **transcript intron 2**, insertion point = junction 69.
A barcode in genomic intron 2 (2082–2092) is **transcript intron 1**, insertion point = junction 9.

### TOY_SINGLE (`NM_900003.1`), intronless

`CAC` + `ATGGCTGCTTAA` + `CACACA`, 21 nt, one exon 3000–3021, CDS 3–15. Its only exon is also the
last exon, so the decision tree always returns "no NMD".

### TOY_SEAM (`NM_900004.1`), plus strand, 3 exons — a stop on the seam

Built so that transcript nt 9 and 10 are `TA` while nt 11 (the first base of exon 2) is `C`:

```
transcript: A T G G C T G C T T A │ C G A G C A G C A │ C T A A C A C A C A C A
position:   0 1 2 3 4 5 6 7 8 9 10│11                 │20
codons:     ATG  GCT  GCT  T A C ← nt 9,10,11 = TAC, an ordinary codon, not a stop
                              └──── exon junction at 11
```

But insert the construct into that intron and something starting with `G` now sits between nt 10 and
11, so the codon becomes `T` `A` `G` = **TAG, a stop straddling the host|insert seam**. This is the
mechanism that produces 20 (5ss) and 15 (CPA) rows in the real data.

---

## Every case, one at a time

21 barcodes. Below: the actual line from `5ss_noP.txt.retained_barcodes_genuine.bed`, then what the
pipeline must produce and why.

### UTR integrations — short-circuited, nothing is computed

| BED line | where it landed | verdict |
|---|---|---|
| `chr1 1010 1011 ACCC 1 +` | TOY_PLUS intron 1, junction 6, ATG at 9 → **above the start** | `utr`, no NMD |
| `chr1 1024 1025 ACCG 1 +` | intron 2, junction **exactly 9 = the ATG** — boundary case | `utr`, no NMD |
| `chr1 1160 1161 CCTG 1 +` | intron 6, junction 106, normal stop ends at 102 → **below the stop** | `utr`, no NMD |
| `chr1 1002 1003 ACGC 1 +` | exon 1, transcript nt 2 → 5'UTR | `utr`, no NMD |
| `chr1 1168 1169 ACTG 1 +` | exon 7, nt 108 → 3'UTR (and it is the last exon) | `utr`, no NMD |

All five have **empty** NMD columns: no PTC arises, so there is nothing to measure. The order of the
tree matters — `ACTG` would also match "last exon", but the UTR check comes first.

### Exonic integrations — by rule, not computed

| BED line | where it landed | verdict |
|---|---|---|
| `chr1 1056 1057 ACGG 1 +` | exon 4 of 7, nt 26, CDS → **internal coding exon** | `internal_exon`, **NMD** |
| `chr1 1150 1151 ACTC 1 +` | exon 6 of 7, nt 100 < 102 → CDS; the stop codon sits at nt 99–101, i.e. **inside exon 6**, so this IS the last coding exon | `last_exon`, no NMD |
| `chr1 4041 4042 TCGC 1 +` | TOY_SEAM exon 3 of 3, nt 21 < 24 → CDS, last exon | `last_exon`, no NMD |
| `chr1 2010 2011 GCTC 1 -` | TOY_MINUS exon 3 of 3 (minus strand!), nt 70 < 72 → CDS | `last_exon`, no NMD |
| `chr1 3006 3007 TCCC 1 +` | TOY_SINGLE, its only exon is the last | `last_exon`, no NMD |

NMD columns are empty here too: how a split host exon actually splices is unknown, so no PTC
position is invented.

### CDS-intronic integrations — the one computed branch

The insert goes in at junction `p`, that junction splits into two (`p` and `p + L`), everything
below shifts by `L`, the first in-frame stop from the CDS start is found, and the distance is the
last junction minus the end of the stop codon.

**`chr1 1070 1071 CCCG 1 +` — frame 0, stop inside the insert**

```
intron 4 -> p = 36.  (36 - 9) % 3 = 0 -> frame 0 -> stop at insert offset 9 -> PTC = 36 + 9 = 45
junctions before: 6, 9, 22, 36, 98, 106
junctions after : 6, 9, 22, 36, 48, 110, 118   <- 36 split into 36 and 48; 98->110, 106->118
last junction 118, end of the stop 45+3 = 48 -> distance 70 > 50 -> NMD
```

**`chr1 1046 1047 CCCC 1 +` — frame 2, the insert is read through**

```
intron 3 -> p = 22.  (22 - 9) % 3 = 1 -> frame 2 -> no stop in the insert in this frame,
and the barcode CCCC supplies none, so the ribosome runs through into the host
variant A: 12 nt, a multiple of 3 -> host frame survives -> the HOST's own stop is found (99 -> 111)
           last junction 118, end of the stop 114 -> distance 4 <= 50 -> no NMD
variant B: 13 nt -> host frame shifts -> no stop anywhere at all -> ptc_source = none, no NMD
```

This is the row the two insert variants exist for.

**`chr1 1048 1049 CTAA 1 +` — same intron, but the barcode supplies the stop**

```
p = 22, frame 2. The codon at insert offset 5 is barcode[1..3] = T A A = a stop.
PTC = 22 + 5 = 27, ptc_source = barcode, bc_has_stop = 1
last junction 118, end of the stop 30 -> distance 88 > 50 -> NMD
```

The same integration site as `CCCC`, opposite verdict — purely because of the barcode's letters.

**`chr1 1142 1143 CCTC 1 +` — frame 1**

```
intron 5 -> p = 98.  (98 - 9) % 3 = 2 -> frame 1 -> stop at offset 1 -> PTC = 99
last junction 118, end of the stop 102 -> distance 16 <= 50 -> no NMD
```

**`chr1 2086 2087 GCCG 1 -` and `chr1 2016 2017 GCCC 1 -` — minus strand, the inversion check**

```
GCCG in genomic intron 2 -> TRANSCRIPT intron 1 -> p = 9
      frame 0 -> PTC = 18; junctions become 9, 21, 81 -> distance 60 > 50 -> NMD
GCCC in genomic intron 1 -> TRANSCRIPT intron 2 -> p = 69
      frame 0 -> PTC = 78; junctions become 9, 69, 81 -> distance 0 -> no NMD
```

Had the intron index been taken in genomic order these two would swap and both verdicts would
invert. This is the main minus-strand check.

**`chr1 4016 4017 TCCG 1 +` — a stop on the seam**

```
TOY_SEAM intron 1 -> p = 11.  (11 - 0) % 3 = 2 -> frame 1
But the codon at position 9 is reached before any insert stop: host nt 9, 10 = T, A
plus the insert's first base G -> TAG. PTC = 9, which is BELOW the insertion point.
ptc_source = host_insert_seam
distance 20 <= 50 -> no NMD
```

### Out of gene and antisense

| BED line | why |
|---|---|
| `chr1 5000 5001 AGGG 1 +` | nowhere — `intergenic`, `in_gene = 0` |
| `chr1 1046 1047 AGGC 1 -` | the same site as `CCCC` but on the minus strand, while TOY_PLUS is plus → `antisense` |

Note that `AGGC` and `CCCC` sit at the same genomic base and differ only in strand.

### Category precedence

| BED line | logFC | category |
|---|---|---|
| `chr1 1070 1071 TCGG 1 +` | `NA` | **1**, even though its placement is identical to `CCCG` (category 4) |
| `chr1 5100 5101 TTTT 1 +` | absent from logFC | **1**, but the row is still in the table |

`TCGG` proves "not expressed" outranks everything: its placement and NMD columns are filled exactly
like `CCCG`'s, yet the category is 1. `TTTT` proves a barcode missing from logFC is not lost — which
matters, because in the real data exactly one barcode per library did drop out of logFC (the BED's
first line was consumed as column names).

---

## What the check prints

```
variant A: 21 barcodes x 21 fields = 441 comparisons
  PASS -- output matches the hand calculation exactly
variant B: 21 barcodes x 21 fields = 441 comparisons
  PASS -- output matches the hand calculation exactly
```

On a mismatch: one line per field, naming the barcode, the column, the expected value and the
observed one.

`expected_A.tsv` and `expected_B.tsv` are written **by hand** from the walkthrough above. They must
never be generated from pipeline output — that would only prove the code agrees with itself. A dot
in a field means "this must be empty". `check.py` separately verifies that every row of an
expectation file has as many fields as its header: a bug in the test must not masquerade as a bug in
the pipeline (that check caught exactly that once).

## One precedent

The first version of this test caught a real bug: in the "no stop found" branch the script failed to
fill `frame_offset`. The fix belonged in the **code**, not in the expectation. When the test fails,
work it out by hand first and decide who is wrong.
