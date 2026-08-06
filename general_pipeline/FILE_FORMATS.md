# TRIP base pipeline — file formats at each step

For every step: what the script does (one line), then each **input** and **output** file — the
name pattern with a concrete example, the format, and the real first lines.

Conventions:
- `.gz` files are passed to the scripts **compressed**; the previews below show them **decompressed**.
- Tables are **tab-separated**. Where a header row exists it is shown; most intermediate files have
  **no header** (columns are described instead).
- Example dataset: the single-exon **5ss_noP** chain (map sample `Map_5ss_noP_pool1_rep2_S166`,
  expression sample `Exp_total_5ss_noP_pool2_rep2_S80`).

---

## 01 — cutadapt (`run_cutadapt.sh`)
Trims adapters, moves the 4-nt UMI and 18-nt barcode into the read name, then extracts a
per-read barcode table.

**INPUT — raw paired FASTQ (`.fastq.gz`)**
`{sample}_R1_001.fastq.gz`, `{sample}_R2_001.fastq.gz`
e.g. `Map_5ss_noP_pool1_rep2_S166_R1_001.fastq.gz`
Standard FASTQ, 4 lines/read (header / sequence / `+` / quality):
```
@LH00211:57:222NKVLT1:1:1101:27573:1064 1:N:0:NCAGTTAG+NGCTAAGT
ATTTAGAAGTTCAACTCGGTACCAGAGCCACCGTTNGTGCCCCAGATTTGGGATCAGGCTGGTCTTGAGCTTCTGACCTCAGG...
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII...
```

**OUTPUT — trimmed FASTQ (`.fastq.gz`)**
`{sample}.R1.fastq.gz`, `{sample}.R2.fastq.gz` e.g. `Map_5ss_noP_pool1_rep2_S166.R1.fastq.gz`
Same FASTQ format, but the header now carries `umi=` / `adapter=` / `barcode=` tags:
```
@LH00211:57:222NKVLT1:1:1101:27573:1064 umi=ATTT adapter=1 barcode=AGAAGTTCAACTCGGTACCAGAGCCACCGTTNGTGCCCCAGATTTGGGATC
AGGCTGGTCTTGAGCTTCTGACCTCAGGTGGNCTGCCTGCCTCAGCCTCCCNAAGTGCTGGAATTACAGGTGTGAGTCACT...
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIII#IIIIIIIIIIIIIIIIIII#IIIIIIIIIIIIIIIIIIIIII9IIIIIII...
```

**OUTPUT — barcode table (`.txt.gz`, decompressed below)**
`{sample}.barcodes.txt.gz` e.g. `Map_5ss_noP_pool1_rep2_S166.barcodes.txt.gz`
3 columns, no header — `read_name` · `UMI (4 nt)` · `barcode (18 nt)`:
```
@LH00211:57:222NKVLT1:1:1101:27573:1064   ATTT   GTTNGTGCCCCAGATTTG
@LH00211:57:222NKVLT1:1:1101:41498:1064   CTTC   AGANTCTTGGTCGAGAGG
@LH00211:57:222NKVLT1:1:1101:19501:1080   ATCG   ATTNTCTCGGGGCGCGGA
```

---

## 02 — bowtie2 (`run_bowtie2.sh`)
Aligns the trimmed reads to hg38, keeps the genomic mate (read 2), converts it to BED.

**INPUT — trimmed FASTQ** — `{sample}.R1.fastq.gz`, `{sample}.R2.fastq.gz`
(identical to the trimmed FASTQ from step 01, shown above).

**OUTPUT — genomic mate as BED12 (`.bed.gz`, decompressed below)**
`{sample}.second_read.bed.gz` e.g. `Map_5ss_noP_pool1_rep2_S166.second_read.bed.gz`
12-column BED (chrom · start · end · read/2 · score · strand · thickStart · thickEnd · rgb ·
blockCount · blockSizes · blockStarts), no header:
```
chr1  28464  28575  LH00211:57:222NKVLT1:1:1153:15526:18833/2  1  -  28464  28575  255,0,0  1  111  0
chr1  28465  28575  LH00211:57:222NKVLT1:1:1114:27740:20387/2  1  -  28465  28575  255,0,0  1  110  0
```

---

## 03 — analyzeTRIPMapping_uniform.R
Groups reads by genomic position, keeps 18-nt barcodes, collapses barcode mutants per position,
drops β-globin self-integrations and noisy positions.

**INPUT** — `{sample}.barcodes.txt.gz` (step 01) + `{sample}.second_read.bed.gz` (step 02).

**OUTPUT — per-position table (`.txt`)**
`{sample}.barcodes_poses.txt` e.g. `Map_5ss_noP_pool1_rep2_S166.barcodes_poses.txt`
9 columns, no header — chrom · start · end · `barcode:count` list · score(1000) · strand ·
`counts` · `cleanCounts` · collapsed barcode(s):
```
chr1  28575  28576  AAATACGGCGGAGAAGGC:27              1000  -  1  1  AAATACGGCGGAGAAGGC
chr1  37578  37579  CTAAACGTAAACTGTCAA:73              1000  -  1  1  CTAAACGTAAACTGTCAA
```

**OUTPUT — per-position BED6 (`.bed`)**
`{sample}.barcodes.bed` e.g. `Map_5ss_noP_pool1_rep2_S166.barcodes.bed`
6 columns, no header — chrom · start · end · `barcode:count(,…)` list · #barcodes · strand:
```
chr1  28575  28576  AAATACGGCGGAGAAGGC:27                                    1  -
chr1  30359  30360  TAAATTGACGCAACGACG:48,GCCCGTATCGACTGTTGC:14              2  +
```

---

## 04 — filterMappings.R
Removes near-duplicate positions (a weaker position within 10 bp of a stronger one whose top
barcode is a mutant of it).

**INPUT** — `{sample}.barcodes_poses.txt` + `{sample}.barcodes.bed` (step 03, shown above).

**OUTPUT — filtered per-position table (`.txt`)**
`{sample}.filtered_barcodes_poses.txt` e.g. `Map_5ss_noP_pool1_rep2_S166.filtered_barcodes_poses.txt`
Same as `barcodes_poses.txt` **plus two columns**: `MaxCount` (V10) and `TopBarcode` (V11) — 11 cols:
```
chr1  28575  28576  AAATACGGCGGAGAAGGC:27  1000  -  1  1  AAATACGGCGGAGAAGGC  27  AAATACGGCGGAGAAGGC
chr1  30359  30360  TAAATTGACGCAACGACG:48,GCCCGTATCGACTGTTGC:14  1000  +  2  2  TAAATTGACGCAACGACG,GCCCGTATCGACTGTTGC  48  TAAATTGACGCAACGACG
```

**OUTPUT — filtered BED6 (`.bed`)**
`{sample}.filtered_barcodes.bed` — same 6-column layout as `barcodes.bed`:
```
chr1  28575  28576  AAATACGGCGGAGAAGGC:27                        1  -
chr1  30359  30360  TAAATTGACGCAACGACG:48,GCCCGTATCGACTGTTGC:14  2  +
```

---

## 05 — filterContaminants_uniform.R
Removes cross-library contamination: for each `chr:start:barcode` key, drops rows whose count is below
1/50 of the same key's maximum seen in libraries of **other** experiment groups (a barcode/position not
seen outside its own group is kept).

**INPUT** — `{sample}.filtered_barcodes_poses.txt` (must be 11 columns) + `{sample}.filtered_barcodes.bed` (step 04).
_(Argument-driven: `Rscript filterContaminants_uniform.R <path> <fileListFile>`, where `fileListFile` has
two columns — library name and GROUP (experiment type); the comparison key is `chr:start:barcode` (cols
1,2,9) and `contFactor = 50`. The superseded `filterContaminants.R` — no arguments, 12 hard-coded
libraries, key = barcode only, `factor = 20`, output to `set1_combined_1/` — is what produced the current
tables shown here.)_

**OUTPUT — retained BED6 (`.bed`)**
`{sample}.retained_barcodes.bed` e.g. `Map_5ss_noP_pool1_rep2_S166.retained_barcodes.bed`
Same 6-column layout; UMI counts are now the deduplicated survivors:
```
chr1  28575  28576  AAATACGGCGGAGAAGGC:25                                          1  -
chr1  30359  30360  TAAATTGACGCAACGACG:50,GCCCGTATCGACTGTTGC:16,TAAATTGACGCGACGACG:1  3  +
```

**OUTPUT — retained counts table (`.txt`)**
`{sample}.retained_barcodes_counts.txt` — 11-column layout (like `filtered_barcodes_poses.txt`).

---

## 06 — get_genuine_barcodes_v2.R
Merges replicate maps, collapses barcode mutants (Hamming ≤ 2), and keeps each barcode at its
single genuine position.

**INPUT — manifest (`.txt`)**
`{name}.txt` e.g. `5ss_noP.txt` — plain list of the `retained_barcodes.bed` files to combine:
```
Map_5ss_noP_rep1_S95.retained_barcodes.bed
Map_5ss_noP_rep2_S96.retained_barcodes.bed
Map_5ss_noP_pool1_rep2_S166.retained_barcodes.bed
Map_5ss_noP_pool2_rep2_S167.retained_barcodes.bed
```
**INPUT** — the `*.retained_barcodes.bed` files it lists (step 05 output; BED6, shown above).

**OUTPUT — genuine barcodes BED6 (`.bed`)**
`{name}.txt.retained_barcodes_genuine.bed` e.g. `5ss_noP.txt.retained_barcodes_genuine.bed`
6 columns, no header — chrom · start · end · **genuine barcode (18 nt)** · count · strand:
```
chr1  113509  113510  TATGGCATGGTGCCGTGC  1  -
chr1  115353  115354  GTACAGTATTATATCAAT  1  -
```

---

## 07 — normalizeTRIP_v2.R
Joins, per barcode, the genomic position (map) + copy number (norm) + expression (exp) into one
table, keeping single-position barcodes only.

**INPUT — driver (`combos.txt`)**
5 columns, no header — `mapName(=manifest)` · `normName` · `expName` · `effLength` · `reverseComplement`:
```
5ss_noP.txt   Norm_5ss_noP_rep2_S92   Exp_total_5ss_noP_pool2_rep2_S80   18   FALSE
```
**INPUT** — `{name}.txt.retained_barcodes_genuine.bed` (step 06).
**INPUT — norm & exp counts (`.counts.txt.gz`, decompressed below)**
`{sample}.counts.txt.gz` e.g. `Norm_5ss_noP_rep2_S92.counts.txt.gz`, `Exp_total_5ss_noP_pool2_rep2_S80.counts.txt.gz`
2 columns, no header — `barcode` · `read_count`:
```
ATTGCCCATCCGGCGCTG   5498
CATAGCGGCGATTGAACG   5274
TAATGTTGACGTCGCTCT   4766
```

**OUTPUT — combined counts (`.combined_counts.txt`)**
`{expName}.combined_counts.txt` e.g. `Exp_total_5ss_noP_pool2_rep2_S80.combined_counts.txt`
**Has a header.** The first field of each data row is the **barcode (row name, unlabeled)**,
followed by the 6 header columns `norm exp Chr Pos Strand map`:
```
                    norm  exp  Chr    Pos       Strand  map
AAAAAACGTACTGCGTCC  207   24   chr14  88508786  -       402
AAAAAATAGCCACGGGCA  291   0    chr9   130784479 -       48
```

**OUTPUT — expression BED6 (`.bed`)**
`{expName}.bed` — chrom · start · end · barcode · exp · strand:
```
chr1  113508  113510  TATGGCATGGTGCCGTGC  0   -
chr1  273364  273366  GCTTAGCCCCACTTGAGG  67  +
```

---

## 08 — make_tracks.sh
Intersects the expression BED with genome annotations (MANE genes/exons/introns, gencode
non-coding, repeats), in sense and antisense, producing one file per track.

**INPUT** — `{expName}.bed` (step 07) + `names.exp.txt` (a one-column list of experiment names):
```
Exp_total_5ss_noP_pool2_rep2_S80
```

**OUTPUT — track overlaps (`.txt`)** — one file per annotation, all `bedtools intersect` output
(the 6 barcode-BED columns, then the annotation record, then a trailing overlap flag), no header.

`{expName}.MANE_all.txt` — barcode-BED (6) + MANE transcript BED12 (12) + flag = **19 cols**:
```
chr1  1553068  1553070  TCTGGTGATGAGGTCCCC  587  -  chr1  1541672  1574863  NM_014188.3  0  -  1542065  1574557  0  5  495,119,140,144,386,  0,2196,3190,23100,32805,  2
```
Same 19-column layout: `MANE_exonic.txt`, `MANE_all_antisense.txt`, `MANE_exonic_antisense.txt`,
`MANE_introns_sense.txt`, `MANE_introns_antisense.txt`.

`{expName}.gencode_noncoding_all.txt` — barcode-BED (6) + gencode record + flag:
```
chr1  113508  113510  TATGGCATGGTGCCGTGC  0  -  chr1  89295  133566  ENSG00000238009.7  .  -  lncRNA  ENSG00000238009  2
```

`{expName}.repeats.txt` — barcode-BED (6) + repeat record (chr,start,end,name,0,strand,class,family) + flag = **15 cols**:
```
chr1  3804100  3804102  CCCTCGTAATAAGAACCC  0  -  chr1  3803658  3805070  HERVK13-int  0  -  LTR  ERVK  2
```

---

## 09 — gc_content_single.sh
Computes GC% in a window around each integration and appends it as a column.

**INPUT** — `{expName}.bed` (step 07, shown above).

**OUTPUT — BED with GC (`_with_gc.bed`)**
`{expName}_with_gc.bed` e.g. `Exp_total_5ss_noP_pool2_rep2_S80_with_gc.bed`
8 columns, no header — the expression BED6 + barcode (repeated, from `getfasta`) + **GC%** (last):
```
chr1  113508  113510  TATGGCATGGTGCCGTGC  0  -  TATGGCATGGTGCCGTGC  36.00
chr1  120244  120246  ACTCATGTGTCAGAATCA  1  +  ACTCATGTGTCAGAATCA  34.00
```

**OUTPUT — GC table (`_gc_content.txt`)**
`{expName}_gc_content.txt` — 2 columns, no header — `barcode` · `GC%`:
```
TATGGCATGGTGCCGTGC   36.00
ACTCATGTGTCAGAATCA   34.00
```

---

## 10 — table_annotation.R
Merges combined counts + all track files + GC into the flagship per-barcode table, adding host-gene
expression, distances, chromatin state, GC category and expression classes (~48 columns).

**INPUT** — `{expName}.combined_counts.txt` (step 07) + every track file (step 08) +
`{expName}_with_gc.bed` (step 09) + `names_all_clean.txt` (one-column list of experiment names).

**OUTPUT — annotated table (`.annotated_tables.txt`) — the final result**
`{expName}.annotated_tables.txt` e.g. `Exp_total_5ss_noP_pool2_rep2_S80.annotated_tables.txt`
**Has a header of 48 columns**; each data row is prefixed by the **barcode (row name, unlabeled)**.
Shown below as `column = value` for the first data row (barcode `AAAAAACGTACTGCGTCC`):

```
norm                             = 207
exp                              = 24
Chr                              = chr14
Pos                              = 88508786
Strand                           = -
map                              = 402
exp.norm                         = -3.05658352836637     # log2((1+exp)/(1+norm))
exp.log10                        = 1.39794000867204      # log10(1+exp)
MANE_all                         = NM_007039.4
MANE_exonic                      = NA
MANE_antisense                   = NA
MANE_exonic_antisense            = NA
MANE_intron                      = NM_007039.4.16
MANE_intron_antisense            = NA
gencode_noncoding                = NA
annot_rna_type                   = MANE_mRNA
gencode_noncoding_name           = NA
MANE_all_symbol                  = PTPN21
MANE_all_exp                     = 302
MANE_exonic_exp                  = NA
MANE_antisense_exp               = NA
MANE_exonic_antisense_exp        = NA
MANE_all_exp.norm                = 46.8725748874748
MANE_all_exp.norm.log10          = 1.68008678719258
MANE_exonic_exp.norm             = NA
MANE_antisense_exp.norm          = NA
MANE_exonic_antisense_exp.norm   = NA
MANE_intronic                    = NM_007039.4
MANElength                       = 89230
MANE_intron_length               = 9070
DistToStart                      = 43009
DistToEnd                        = 46221
DistToTSS                        = 46221
DistToTES                        = 43009
RelPos                           = 0.517998431020957
DistToIntronStart                = 765
DistToIntronEnd                  = 8305
DistToIntron5                    = 8305
DistToIntron3                    = 765
RelIntronPos                     = 0.915656008820287
gencodelength                    = NA
ChromHMM                         = Quies
GC                               = 38
GC_category                      = Average GC
RepeatGroup                      = NA
residuals                        = 0.607151975498534
category_5classes                = Class4
category_3classes                = Class3
```
_(This barcode sits in an intron of PTPN21: `MANE_intron`/`MANE_intronic` are set, `MANE_exonic` is NA.)_
