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



## Contents

| Path | What it is |
|---|---|
| [scripts/](scripts/) | one script per stage, numbered as in `FILE_FORMATS.md` |
| [FILE_FORMATS.md](FILE_FORMATS.md) | every input and output file of every stage: name pattern, column layout, real first lines |
| [LAB_NOTEBOOK.md](LAB_NOTEBOOK.md) | condensed notebook for the two new stages: the reasoning, what the toy tests cover, and the real-data measurements |
| [figures/](figures/) | the eight figures of those real-data runs, embedded in the notebook |
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

