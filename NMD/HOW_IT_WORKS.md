# How the NMD pipeline works

What the scripts read, what comes out, and what every column of the result means. All examples are
real lines from real files.

[the output table](#the-output-table) · [inputs](#inputs) · [the decision tree](#the-decision-tree) ·
[the four stages](#the-four-stages) · [the toy test](#the-toy-test) · [checks](#what-this-was-checked-against) · [how to run](#how-to-run)

---

## The output table

`output/nmd_table_5ss_noP.tsv` (58172 rows) and `nmd_table_CPA_noP.tsv` (57303 rows). One row per
barcode, 34 tab-separated columns, nothing filtered out.

**An empty cell means "not applicable", never zero.** A barcode outside a gene has no exon index; a
barcode we refuse to compute a stop for has no distance. Read a blank as "we are not saying anything
here."

| # | column | meaning and possible values |
|---|---|---|
| 1 | `barcode` | 18 nt of `ACGT`. The key every join runs on. Never repeated. |
| 2 | `chrom` | `chr1`…`chr22`, `chrX`, `chrY`, or a `chrUn_*` scaffold |
| 3 | `pos` | integration coordinate, 0-based, as in the BED |
| 4 | `strand` | `+` or `-` |
| 5 | `score` | read count from BED column 5 |
| 6 | `tx_id` | the MANE transcript we picked, e.g. `NM_030937.6`. Empty out of gene |
| 7 | `gene_symbol` | e.g. `CCNL2`, looked up from the DESeq2 table. Empty out of gene |
| 8 | `in_gene` | `1` only for an overlap **on the same strand**, else `0` |
| 9 | `excluded_reason` | why it was not analysed: `intergenic` / `antisense` / `unplaced_contig`. Empty when `in_gene = 1` |
| 10 | `n_tx_overlap` | how many same-strand MANE genes contained it. `>1` = a tie, broken by transcript length |
| 11 | `loc` | `exonic` / `intronic`. Empty out of gene |
| 12 | `feat_idx` | which exon or intron, 1-based, **in transcript order** |
| 13 | `n_exons` | exons in the transcript |
| 14 | `cds_last_exon` | the exon **holding the host's normal stop codon**. `feat_idx == cds_last_exon` is what "last exon" means in the rule — not `feat_idx == n_exons` |
| 15 | `tx_pos` | insertion point in transcript coordinates; for an intronic hit, the junction position |
| 16 | `region` | `5utr` / `cds` / `3utr`, from comparing `tx_pos` with the CDS bounds |
| 17 | `insert_frame` | `0` / `1` / `2` — which phase the ribosome enters the insert in |
| 18 | `insert_len` | `400`. Written on every computed row so the table cannot be confused with a run on a different insert |
| 19 | `bc_has_stop` | `1` if the barcode's own 18 nt carry an in-frame stop. True for about a fifth of computed rows |
| 20 | `ptc_source` | which part of the mRNA the stop came from: `insert_5p` / `barcode` / `host_insert_seam` / `host_downstream` / `none`. The code can also emit `insert_3p`, `host_upstream` and a `_seam` suffix; none of those occur on this data |
| 21 | `ptc_tx_pos` | stop position in the **modified** transcript |
| 22 | `last_junction_pos` | last junction between **coding** exons of the modified transcript. A 3'UTR junction never counts |
| 23 | `n_junc_downstream` | coding junctions at or below the end of the stop codon |
| 24 | `dist_to_last_junction` | `last_junction_pos − (ptc_tx_pos + 3)`. Negative = the stop is past that junction |
| 25 | `dist_gt_50` | `1` when column 24 exceeds 50. This is the 50-nt rule itself |
| 26 | `selenocysteine` | `1` if a read-through selenoprotein `TGA` was skipped on the way |
| 27 | `note` | plain-text reason; one of five fixed strings, filled on the branches that were answered structurally |
| 28 | `CHXbaseMean` | **how abundant the barcode is** — mean read count over all four samples (2 CHX, 2 EtOH) after correcting for depth. **Not** a fold change. Median 93 (5ss) / 166 (CPA) |
| 29 | `CHXlog2FoldChange` | CHX vs EtOH, log2. Positive = stabilised by CHX. **Empty means not expressed** — the barcode never reached the ≥80-read cut before DESeq2 ran |
| 30 | `CHXlfcSE` | standard error of column 29 |
| 31 | `CHXpadj` | adjusted p-value. Also empty on rows where DESeq2's independent filtering dropped it, though those barcodes are expressed |
| 32 | `nmd_expected` | the prediction, `0` or `1`. Empty on rows we do not analyse |
| 33 | `verdict_reason` | which branch answered: `not_analysed` / `utr` / `last_exon` / `internal_exon` / `cds_intron` / `no_ptc_frameshift` |
| 34 | `category` | 1 not expressed · 2 out of gene or antisense · 3 no NMD · 4 NMD expected |

Columns 17–26 are filled **only** on the two branches that read sequence, and not even fully there:

| | 17–22, 26 (frame, insert, source, positions) | 23–25 (the distance and the threshold) | 27 `note` |
|---|---|---|---|
| `cds_intron` | filled | filled | empty |
| `no_ptc_frameshift` | filled | **empty** — no premature stop, so the rule has nothing to measure | filled |
| `utr` / `last_exon` / `internal_exon` | empty | empty | filled |
| `not_analysed` | empty | empty | empty |

The table never carries a number describing a molecule we cannot pin down. `dist_gt_50` and
`nmd_expected` are the same number wherever both are filled; the duplication is deliberate, so you can
see the threshold did not get overridden elsewhere.

### One row read end to end

```
barcode=TTTCCACCCGAGATTGTA  chrom=chr1  pos=1390629  strand=-  gene_symbol=CCNL2
loc=intronic  feat_idx=6  n_exons=11  cds_last_exon=11  tx_pos=788  region=cds
insert_frame=0  ptc_source=insert_5p  ptc_tx_pos=842  last_junction_pos=1640
dist_to_last_junction=795  dist_gt_50=1  nmd_expected=1  verdict_reason=cds_intron  category=4
```

An integration into intron 6 (of 10) of CCNL2 on the minus strand, at nt 788 of the transcript, inside
the coding part. The ribosome enters the insert in frame 0 and the first stop is at nt 54 of the
insert — check: `842 − 788 = 54`, which is exactly the frame-0 first-stop position documented for our
400-nt insert. The last coding junction is 795 nt further down, so NMD is expected.

---

## Inputs

**Barcodes** — `genuine_barcodes/5ss_noP.txt.retained_barcodes_genuine.bed`, BED6, one 1-nt interval
per integration:

```
chr1	37578	37579	AATACATGATGAATCCCG	1	-
```

**Gene annotation** — `genuine_barcodes/MANE.GRCh38.v1.0.refseq_genomic.bed`, BED12, 19120 transcripts:

```
chr1	65418	71585	NM_001005484.2	0	+	65564	70008	0	3	15,54,2549,	0,101,3618,
```

`chromStart`/`chromEnd` bound **transcription** (the whole transcript); `thickStart`/`thickEnd` bound
**translation** (the CDS only). They differ because the transcript starts before the `ATG` — the 5'UTR
is transcribed but not translated. The last column, `0,101,3618`, holds exon offsets **from
`chromStart`**, not from the CDS. In 37.1% of MANE transcripts the `ATG` is not even in exon 1.

**Transcript sequences** — `MANE.GRCh38.v1.0.refseq_rna.fna`, the **mature mRNA**: introns already
gone, already 5'→3' (so for a minus-strand gene, the reverse complement of the genome). Every record's
length equals the sum of its BED12 exon lengths, which is what makes the coordinate projection exact.

**The insert** — `genuine_barcodes/exon_and_5arm_sequence.txt`, 400 nt on one line: 75 nt of cryptic
exon from the 5' arm plus 325 nt of β-globin exon, with a run of 18 `N` where the real barcode goes.
The script finds that slot itself and reads the length from the file, so nothing about the insert is
hard-coded — it has already been swapped once (325 → 400 nt) with no code change.

**The measurement** — `logFC/CHX_{5ss,CPA}.txt`, the DESeq2 result. We take the four `CHX*` columns and
the accession → symbol lookup, and nothing else: the annotation columns in that file disagree with
MANE v1.0, and `MANE_intron_ind` numbers introns in genomic order, which is inverted on the minus
strand. We join **by barcode**, never by row order.

---

## The decision tree

> **A premature stop more than 50 nt upstream of the last junction BETWEEN CODING EXONS triggers NMD.
> Closer than that, or past it, the transcript escapes.**

Coding exons are 1 … *k*, where exon *k* holds the host's own stop codon. The junction that counts is
the one joining the last two of them — not the transcript's last junction. They differ in the 837 MANE
genes (4.4%) with an intron in the 3'UTR.

Applied in order; first match wins:

| where the integration landed | verdict |
|---|---|
| outside a MANE gene, or on the opposite strand | not analysed |
| **any UTR** — 5' or 3', exonic or intronic | no NMD — no premature stop arises at all. Checked FIRST, before exon vs intron is even asked |
| the last coding exon | no NMD — translation already ends there |
| any earlier coding exon | **NMD** |
| an intron inside the CDS | **computed** — build the modified mRNA and measure |

Exonic hits go by rule rather than calculation because the construct **breaks** the host exon and how
that splices does not follow from the coordinates: the whole construct might stay, half the exon might
leave with a piggyBac ITR, or the exon might be replaced. So we do not invent a stop position.

CDS introns are computed because there the structure is unambiguous: the construct's flanks are
introns, they leave with the host intron, exactly 400 nt reach the mRNA, and one host junction becomes
two.

---

## The four stages

Python 3, standard library only. Each takes `--root DIR` and `--libs a,b`, which is why the toy test
drives the production scripts rather than a copy.

**Stage 1 — `01_build_transcripts.py`.** MANE BED12 + FASTA → one row per transcript in **transcript
coordinates**: exon lengths, junction positions (running sums of exon lengths), CDS bounds, and
`cds_last_exon`. Distances have to be measured in the mRNA, not the genome, or a 20 kb intron wrecks
the arithmetic. The minus-strand mirror lives in exactly one function, `genomic_to_tx` in `common.py`;
every other stage calls it.

Five invariants must hold, and do on all 19120 transcripts: FASTA length equals the block sum, every
transcript has a sequence, no CDS bound falls in an intron, the CDS divides by 3, the CDS ends in a
stop. Two things must *not* fail the run: 33 transcripts with a non-canonical start codon, and 24
selenoproteins with an internal in-frame `TGA` that is read through — their positions are recorded so
stage 3 does not mistake them for a stop.

**Stage 2 — `02_map_barcodes.py`.** Places each barcode: which gene (indexed by chromosome, strand and
100 kb bin), exon or intron, which one **in transcript order**, the insertion point, and the region.
For an intronic hit the insertion point is the **junction**, not a position inside the intron — the
intron is not in the mature mRNA. Where inside the intron the barcode actually sat is deliberately
ignored.

Three distinct reasons a barcode is not analysed: `intergenic` (30553 / 29745), `antisense`
(13761 / 13781), `unplaced_contig` (116 / 109 — calling a `chrUn_*` scaffold "intergenic" would claim
more than we know). Overlapping genes: 379 MANE transcripts (2.0%) overlap another on the same strand;
the longest wins and `n_tx_overlap` flags the row.

For 5ss, 13742 barcodes (23.6%) land inside a gene on their own strand: `intronic/cds` **11006**,
`intronic/5utr` 1981, `intronic/3utr` 31, `exonic/cds` 167, `exonic/5utr` 67, `exonic/3utr` 490. CPA is
similar with `intronic/cds` 10996. Those ~11 000 rows are the computed branch.

**Stage 3 — `03_call_nmd.py`.** The tree. Only the CDS-intron branch reads sequence:

```python
insert   = template[:slot] + barcode + template[slot + bc_len:]
new_seq  = seq[:p] + insert + seq[p:]                  # p is the junction position
new_junc = sorted([j if j < p else j + L for j in coding_junctions] + [p])
stop     = first in-frame TAA/TAG/TGA from the CDS start, skipping selenoprotein TGAs
dist     = max(new_junc) - (stop + 3)
```

One host junction becomes two, at `p` and `p + L`. The frame the insert is read in is set by the
integration point, not the design: `frame = (-(p - cds_start)) % 3`. For our insert the first stop is
at nt 54 in frame 0, nt 23 in frame 2, and **there is none at all in frame 1** — a sequence designed to
be stop-free is stop-free in one phase out of three. The observed frame split is 43.4 / 21.8 / 34.7 %,
not one third each, because junctions are running sums of exon lengths and human exon lengths are
skewed modulo 3.

Two outcomes are not premature stops, and the branch splits accordingly:

- `no_ptc_frameshift` — 400 is not a multiple of 3, so below the insert the host is read out of frame,
  and that can carry the ribosome **past** where its own stop codon sits. Termination is then late, not
  early. 28 rows in 5ss, 42 in CPA, all frame 1. One CPA row finds no stop anywhere at all.
- `cds_intron` — a genuine premature stop; the distance is measured and the threshold applied.

On real data the 50-nt threshold decides 30 rows (5ss) and 33 (CPA) out of ~11 000: the 5th percentile
of the distance is 343 nt and the median ~1290. What decides is the topology, not the threshold.

Where the stop comes from, `ptc_source`, 5ss / CPA: `insert_5p` 8603 / 8681, `host_downstream`
1864 / 1799, `barcode` 519 / 500, `host_insert_seam` 20 / 15. The seam case is a codon straddling the
host|insert boundary: host `TA` + the insert's first base `G` = `TAG`. It only exists in frame 1,
because a codon at `p−2` is in frame exactly when `(p − cds_start) % 3 == 2`.

**Stage 4 — `04_make_table.py`.** Joins placement + verdict + the CHX measurement **by barcode**. The
join is guarded: a barcode appearing twice in any source stops the script rather than silently losing a
row. Categories:

| | condition | 5ss | CPA |
|---|---|---|---|
| 1 | `CHXlog2FoldChange` empty → not expressed | 51871 | 47694 |
| 2 | out of gene, or in a gene but antisense | 2527 | 4803 |
| 3 | in a gene, NMD not expected | 840 | 1022 |
| 4 | in a gene, NMD expected | 2934 | 3784 |

Category 1 beats everything, but the placement columns are still filled — stage 2 knows nothing about
CHX, so even an unexpressed barcode still shows which intron of which gene it sat in.

---

## The toy test

An invented genome small enough to work out in your head, run through **the same production scripts**
with `--root`, and compared field by field against a hand-written expectation. That expectation must
never be generated from the pipeline output, or the test only proves the code agrees with itself.

```bash
bash toy_test/run_toy.sh     # 441 comparisons per variant, both must PASS
```

Two roots, identical except for the insert: `root_A` has 12 nt (a multiple of 3, so a read-through does
not shift the host frame) and `root_B` has 13 nt (not a multiple of 3, so it does). **Everything below
is `root_B`, because that is how the real insert behaves — 400 nt, and `400 mod 3 = 1`.** `root_A`
exercises a branch that never fires on real data.

The toy insert is `GTAA[NNNN]CTAGC`: a stop in two frames out of three (offset 9 in frame 0, offset 1 in
frame 1, none in frame 2 unless the 4-nt barcode brings one) — the same property the real insert has.

Four genes, each carrying one thing that has to be checked:

- **TOY_PLUS**, plus strand, 7 exons, junctions at 6, 9, 22, 36, 98, 106, CDS 9…102, stop at nt 99–101.
  Two 5'UTR introns (one with its junction **exactly on the ATG**), three CDS introns in different
  frames, and a 3'UTR intron — which makes `cds_last_exon = 6` while `n_exons = 7`, the case that
  separates the two.
- **TOY_MINUS**, minus strand, 3 exons — the index mirror: genomic intron 1 is transcript intron 2.
- **TOY_SINGLE**, intronless — its only exon is also the last one.
- **TOY_SEAM**, built so nt 9,10 are `TA` and nt 11 is `C`; drop the construct into that intron and the
  codon becomes `TAG`, straddling the seam.

### What happens to each of the 21 barcodes

The full result is committed alongside the inputs — see
[`toy_test/root_B/output/nmd_table_5ss_noP.tsv`](toy_test/root_B/output/nmd_table_5ss_noP.tsv), one row
per toy barcode with all 34 columns, and the same for `root_A`.

The rows worth understanding, and why they come out that way:

- **`ACCG`** — the junction is at nt 9 and the `ATG` is at nt 9 too. The insert goes in *before* the
  first `A` of the start codon, so it is still above it: `5utr`, no NMD. The boundary case.
- **`ACTG`** — exon 7, which is both 3'UTR and the last exon. Both would give "no", but the UTR check
  sits higher in the tree, so the label is `utr`. The reason is the absence of a stop, not the absence
  of a junction, and the label should say so.
- **`ACTC`** — exon 6 of 7, still coding (nt 100 < 102). Exon 6 is not the transcript's last, but it
  *is* the one holding the stop codon, so `cds_last_exon = 6` and this is `last_exon`, no NMD. Computed
  as `feat_idx == n_exons` it would have gone the opposite way.
- **`CCCC` vs `CTAA`** — the same integration site, opposite verdicts. In frame 2 the insert has no stop
  of its own; `CCCC` brings none either, so the ribosome reads through and the 13-nt frameshift carries
  it past the host's own stop → `no_ptc_frameshift`, no premature stop at all. `CTAA` has `TAA` in the
  barcode itself, giving a stop at offset 5 → a real PTC, dist 81, NMD. This is why the barcode is
  substituted into the slot instead of being ignored.
- **`CCTC`** — frame 1, stop at offset 1 of the insert, dist 9. Here the 50-nt threshold genuinely
  decides: there is a PTC and there are junctions below it, but the last one is too close.
- **`GCCG` and `GCCC`** — the minus-strand check. `GCCG` sits in genomic intron 2, which is
  **transcript intron 1**, and `GCCC` in genomic intron 1, which is transcript intron 2. Taken in
  genomic order the two would swap and both verdicts would flip. That is the mistake baked into
  `MANE_intron_ind`.
- **`TCCG`** — the stop is at nt 9, *below* the insertion point at 11: host `TA` plus the insert's `G`.
  `host_insert_seam`.
- **`AGGC`** — the same genomic position as `CCCC` but on the minus strand, so `antisense`. The whole
  difference is one column.
- **`TCGG`** — the same site as `CCCG`, but with no fold change. Placement and NMD columns are filled
  identically, yet the category is 1: "not expressed" wins without erasing the placement.
- **`TTTT`** — absent from the toy logFC entirely, and still present in the output. That matters
  because in the real data exactly one barcode per library dropped out of logFC.

---

## What this was checked against

**Input invariants.** Five checks in stage 1, all passing on all 19120 transcripts. Any failure exits
with code 1; nothing passes silently.

**The toy test.** 441 comparisons per variant against a hand calculation. Covers both strands, an
intronless gene, all three frames, all three regions, a junction exactly on the `ATG`, the
`cds_last_exon` ≠ `n_exons` case, a stop on the seam, a stop from the barcode, read-through with and
without a frameshift, no stop found, category precedence, and a barcode missing from logFC.

**Cross-check against an independent annotation** — the strongest one, since ours is re-derived from
the BED12 while the other comes from a separate pipeline:

| compared | 5ss (58171 rows) | CPA (57302 rows) |
|---|---|---|
| gene | matched 58111 (**99.90%**) | matched 57262 (**99.93%**) |
| exon/intron | matched all but 5 | matched **all** |
| intron index | exact on the plus strand, an exact mirror (`his = n − ours`) on the minus | same |

Every disagreement is a row with `n_tx_overlap > 1` — the overlapping-gene tie-break, nothing else. On
those rows that annotation contradicts itself: at `chr1:52033516`, `MANE_exonic` names `NM_138417.3`
while `MANE_intron` names `NM_015913.4.5`.

**Format checks.** All barcodes are 18 nt of `ACGT` with no repeats; MANE sequences contain no `N`;
chromosome naming matches on both sides. That matters because an `N` would never match a stop codon,
so the search would quietly run past a real PTC instead of failing. Stage 3 now checks for it
explicitly.

**Two bugs found.** In my code: seam stops were labelled `host_upstream` because the
boundary-straddling check only ran inside the insert — fixed to `host_insert_seam`, no verdict changed.
In the data: the barcode missing from `logFC` is exactly the **first line** of the genuine BED in both
libraries, and logFC's first data row is the BED's second — a headerless BED read with `header=TRUE`,
eating one integration per library as column names. We keep them; they land in category 1.

---

## How to run

The toy test runs as it stands — everything it needs is in this repository:

```bash
cd NMD
bash toy_test/run_toy.sh
```

Running on real data needs input files that are **not** in this repository: the barcode BEDs, the MANE
BED12 and FASTA, and the DESeq2 table. Put them in a directory laid out like the toy roots
(`genuine_barcodes/`, `logFC/`) and point each stage at it with `--root`:

```bash
python3 scripts/01_build_transcripts.py --root DIR   # -> DIR/work/transcripts.tsv
python3 scripts/02_map_barcodes.py      --root DIR   # -> DIR/work/barcodes_<lib>.tsv
python3 scripts/03_call_nmd.py          --root DIR   # -> DIR/work/nmd_<lib>.tsv
python3 scripts/04_make_table.py        --root DIR   # -> DIR/output/nmd_table_<lib>.tsv
```

Each stage prints a report with numbers and exits with code 1 if an invariant breaks. The layering
means a change to the annotation needs a rerun from stage 1, a change to the insert only from stage 3,
and a change to the measurement table only stage 4.

See also: `toy_test/README.md`, and `SUMMARY.md` for the comparison against the CHX measurement.
