# TRIP — primary data-processing pipeline

Scripts that turn raw TRIP sequencing reads into per-barcode tables, plus a small real
example dataset the pipeline can be run on.

**TRIP** (Thousands of Reporters Integrated in Parallel) is a barcoded-reporter assay. A
reporter construct is integrated at thousands of random genomic positions in K562 cells with
**piggyBac** transposase (which inserts at **TTAA** motifs), and every integration carries a
unique **18 bp barcode**. Sequencing then links each barcode to

- its **genomic position** — from DNA "Map" libraries (inverse PCR after DpnII), and
- its **RNA readout** — expression, splicing, polyadenylation, from the same barcode.

Because the reporter sequence is constant, differences in its behaviour are attributed to the
**local genomic context** (chromatin state, GC content, host-gene transcription, distance to
gene and intron features, proximity to nuclear speckles or the lamina).

This repository holds the **primary processing** part: everything from trimmed reads to the
flagship annotated per-barcode table. Downstream analysis (IntronTRIP splicing, U1/CHX, SON,
cDNA, ML) is out of scope.

## Pipeline

```
RAW PAIRED-END FASTQ   (R1 = UMI + primer + barcode, R2 = genomic DNA / transcript)
        |
   01  cutadapt              trim, move UMI+barcode into the read header
        |                    Map libraries -> .barcodes.txt.gz
        |                    Norm/Exp libraries -> .counts.txt.gz  ------------+
   02  bowtie2               align R2 to hg38, keep read 2 -> .second_read.bed.gz
        |                                                                     |
   03  analyzeTRIPMapping    group reads into integration positions,          |
        |                    collapse barcode mutants per position            |
   04  filterMappings        drop near-duplicate positions (+ optional TTAA)  |
        |                                                                     |
   05  filterContaminants    remove cross-library index hopping               |
        |                                                                     |
   06  get_genuine_barcodes  merge replicates, keep single-position barcodes  |
        |                                                                     |
   07  normalizeTRIP  <------ join Map + Norm + Exp per barcode <-------------+
        |                     -> {exp}.combined_counts.txt
   08  make_tracks           intersect with MANE / gencode / repeats
   09  gc_content            local GC% per integration
        |
   10  table_annotation      -> {exp}.annotated_tables.txt   <-- the flagship table
```

## Contents

| Path | What it is |
|---|---|
| [scripts/](scripts/) | one script per stage, numbered as in `FILE_FORMATS.md` |
| [FILE_FORMATS.md](FILE_FORMATS.md) | every input and output file of every stage: name pattern, column layout, real first lines |
| [LAB_NOTEBOOK.md](LAB_NOTEBOOK.md) | working notebook for the two new stages: the reasoning, the toy tests case by case, and the real-data measurements |
| [figures/](figures/) | the seven figures of those real-data runs, embedded in the notebook |
| [example_data/](example_data/) | small real slice + `run_example.sh` (runs stages 03 → 07 in ~2 min) |

### Scripts

| Stage | Script | Interface |
|---|---|---|
| 01 | [scripts/01_run_cutadapt.sh](scripts/01_run_cutadapt.sh) | `bash 01_run_cutadapt.sh map\|count <name> [threads]` |
| 02 | [scripts/02_run_bowtie2.sh](scripts/02_run_bowtie2.sh) | `bash 02_run_bowtie2.sh <name> [threads] [index]` |
| 03 | [scripts/03_analyzeTRIPMapping_uniform.R](scripts/03_analyzeTRIPMapping_uniform.R) | `Rscript … <path> <path_to_mapping> <name>` |
| 04 | [scripts/04_filterMappings_TTAA.R](scripts/04_filterMappings_TTAA.R) | `Rscript … <path> <name> [ttaaFilter <ttaa_bed>]` |
| 05 | [scripts/05_filterContaminants_uniform.R](scripts/05_filterContaminants_uniform.R) | `Rscript … <path> <fileListFile>` |
| 06 | [scripts/06_get_genuine_barcodes_v2.R](scripts/06_get_genuine_barcodes_v2.R) | `Rscript … <path> <nameMap> <barcodeLen>` |
| 07 | [scripts/07_normalizeTRIP_v2_MM.R](scripts/07_normalizeTRIP_v2_MM.R) | `Rscript … <path> <path_norm> <nameMap> <nameNorm> <nameExp> <effLength> <doReverseComplement> [oneMismatch]` |
| 08 | [scripts/08_make_tracks.sh](scripts/08_make_tracks.sh) | `bash 08_make_tracks.sh <names_file>` |
| 09 | [scripts/09_gc_content_single.sh](scripts/09_gc_content_single.sh) | `bash 09_gc_content_single.sh <name>` |
| 10 | [scripts/10_table_annotation.R](scripts/10_table_annotation.R) | `Rscript … <path> <namesExp_file>` |

Stages 03–07 and 10 are R; 01, 02, 08, 09 are shell. All are argument-driven and process one
library (or one driver row) per call; in production they are submitted as LSF jobs.

## Two stages differ from the lab originals

Both are strict supersets of the script they replace: **without the extra flag they behave
exactly like the original**, so they are drop-in replacements. New code is delimited by
`# ===== NEW BLOCK` comments.

**Stage 04 — optional TTAA filter.** piggyBac integrates at a genomic TTAA and duplicates it,
so a genuine integration position should sit next to a TTAA in the reference. With
`ttaaFilter <hg38_ttaa.bed>` the script keeps only those positions. The check is exact and
strand-aware: on `+` the TTAA occupies `[Start-4, Start)`, on `−` it occupies
`[Start, Start+4)`. On 10 real Map libraries (568,687 positions, ~27.9 M reads) it removed
**29 % of barcodes but only 7.6 % of reads** — the dropped positions are overwhelmingly
low-count noise, as expected. The reference (627 MB, 19.8 M sites) is not in this repository.

**Stage 07 — optional 1-nt mismatch tolerance.** The original joins Map, Norm and Exp by exact
barcode identity, so a single sequencing typo loses the read. With the 8th argument
`oneMismatch`, Norm/Exp barcodes exactly one substitution away from a genuine Map barcode are
snapped onto it (neighbour-hash lookup). Rules: genuine Map barcodes are the reference, exact
matches win over 1-nt ones, and a read that is 1 nt from **more than one** genuine barcode is
dropped rather than assigned arbitrarily. The tolerance is fixed at one substitution. On the
example data it recovers 233 extra barcodes (+3 %) and ~10 % more expression reads.

## Running

```bash
module load R/4.3.1-foss-2022b            # R 4.3 + the packages below
bash example_data/run_example.sh          # stages 03 -> 07 on the bundled example
```

R packages: `DescTools` (barcode edit distances), `dplyr`, `purrr`, `tidyr`, `ggplot2`,
`cowplot`, `gridExtra`. Command-line tools: `cutadapt`, `bowtie2`, `samtools`, `bedtools`, and
`bigWigAverageOverBed` for track extraction.

Two practical notes carried over from production: several stages are heavy (a full Map library
takes tens of minutes and several GB), so they are meant to be submitted as cluster jobs, not
run on a login node; and every script writes its outputs **into the directory given as its
first argument**, next to its inputs.

## Known version mismatch (stage 03 → 05)

The current copy of stage 03 writes a 7th column (alignment `Length`) into
`*.barcodes_poses.txt`, whereas the production files in the data tree have 9 columns without
it. Stage 05 keys on columns 1, 2 and 9 (chrom, start, barcodes) and aborts unless its input
has exactly 11 columns, so the extra column both shifts that key onto `cleanCounts` and stops
the run. `example_data/run_example.sh` therefore drops column 7 between stages 03 and 04, with
a comment at that line. This needs a decision from the authors about which layout is canonical;
until then the adapter keeps the chain faithful to the production files.

## Plan for this repository

Written so far: the ten pipeline scripts, the file-format reference, the lab notebook covering the
two new stages, and a runnable example. Still to add:

- **PIPELINE.md** — per-stage narrative: what each step does biologically, why the thresholds
  are what they are, and how the two branches (position vs counts) merge.
- **PROVENANCE.md** — for every script: which read-only lab path it was copied from, its
  timestamp and checksum, and whether it is verbatim or modified.
- **VERSIONS.md** — the version landscape: which older variants exist, what the two coexisting
  normalization lineages are (6-column genuine BED with the barcode in column 4 vs 9-column
  with it in column 7), and the `exp.log10` ambiguity (log10 expression for single-exon
  reporters, but a 0–1 splicing efficiency in the IntronTRIP branch).
- **Tests** — the toy datasets and oracle checks for the two new stages, including the proof
  that flag-off output is byte-identical to the original. The notebook describes them case by
  case; the datasets and drivers themselves are not here yet.
- **Analysis bundles** — the `summary.csv` tables and the run harness behind the figures in
  [figures/](figures/) (TTAA: 29 % of barcodes / 7.6 % of reads; `oneMismatch`: +2.6 % barcodes,
  +12 % exp reads).
- **Reference bundle** — how to build or fetch `hg38_ttaa.bed`, MANE, gencode non-coding,
  repeats, ChromHMM and the K562 host-expression table, so stages 08–10 become runnable too.
- **Thresholds appendix** — the magic numbers in one table (18-nt barcode, alignment length
  20–300, β-globin `chr11:5,220,000–5,240,000` removed, mutant-collapse edit distance ≤ 2–3,
  proximity filter 10 bp / 0.25, contaminant factor 50, `minNormReads` 2, splicing NA below
  20 reads).

## Provenance and attribution

The pipeline scripts were written in the **Ulitsky lab** (Weizmann Institute of Science) by
Igor Ulitsky and Lidiia (`lidiag`); stages 04 and 07 are extensions written by Artyom Tyukaev
on top of their originals. Copies here are verbatim except for those two stages and the two
shell wrappers for stages 01–02, which were reassembled from the per-sample cluster job
scripts. The example data are a slice of unpublished lab sequencing data — please do not
redistribute them outside this repository.
