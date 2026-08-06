# Example data

A small **real** slice of the TRIP data, big enough to run stages 03 → 07 end to end and
small enough for git (~1.5 MB total). Run it with:

```bash
bash example_data/run_example.sh        # ~2 min, writes everything to example_data/out/
```

## What is in here

| File | Stage that reads it | What it is |
|---|---|---|
| `Map_5ss_noP_pool1_rep2_S166.barcodes.txt.gz` | 03 | first 25,000 reads of a real Map library: read name, UMI, 18-nt barcode |
| `Map_5ss_noP_pool1_rep2_S166.second_read.bed.gz` | 03 | the R2 alignments (BED12) of exactly those reads |
| `Map_CPA_noP_pool1_rep2_S164.*` | 03 | the same slice of a second Map library, from a different construct |
| `fileList.txt` | 05 | library name + experiment group; the contaminant filter compares a library only against **other** groups |
| `5ss_noP.txt` | 06 | manifest of the `*.retained_barcodes.bed` files to merge into one position map |
| `combos.txt` | 07 | driver row: `mapName normName expName effLength reverseComplement` |
| `Norm_5ss_noP_rep2_S92.counts.txt.gz` | 07 | gDNA barcode counts (copy-number denominator) |
| `Exp_total_5ss_noP_pool2_rep2_S80.counts.txt.gz` | 07 | RNA barcode counts (expression readout) |

Provenance: the Map slices come from `shared/data/human/TRIP/run5/`, the count tables from
`shared/data/human/TRIP/set1_combined_1/`, and `combos.txt` is the real driver row used for
this sample chain. It is the same library chain used as the worked example throughout
[../FILE_FORMATS.md](../FILE_FORMATS.md), so the file previews there match these files.

The two count tables are filtered to barcodes that are within one substitution of a genuine
barcode of this slice. That keeps them small **and** keeps the sequencing-typo barcodes that
the optional `oneMismatch` mode of stage 07 is meant to rescue.

## What the run prints

```
Stage 03   8,881 / 8,579 integration positions
Stage 04   8,840 / 8,521 positions kept (near-duplicates dropped)
Stage 05   8,840 / 8,521 retained (nothing removed: see note below)
Stage 06   7,998 genuine barcodes, 544 multimappers
Stage 07   7,012 barcodes exact match
           7,245 barcodes with the oneMismatch flag  (+233, +4 % norm reads, +10 % exp reads)
```

Two things to expect on a slice this small:

- **Stage 05 removes nothing.** Index hopping is rare, and a barcode is only judged when the
  same `chr:start:barcode` also shows up in another experiment group. In a 25,000-read slice
  almost every barcode is unseen outside its own group, so it is kept by design (the log says
  `not observed outside group`). On full libraries this filter does remove positions.
- **Stage 04's TTAA mode is off.** Turning it on needs `hg38_ttaa.bed` (627 MB, 19.8 M sites),
  which is not in this repository. See the note in [../README.md](../README.md).

## Stages that cannot run from this folder

- **01 / 02** need the raw FASTQ files and the bowtie2 hg38 index. FASTQs are deliberately not
  included (hundreds of MB each).
- **08 / 09 / 10** need external reference data: MANE annotation, gencode non-coding, hg38
  repeats, ChromHMM, `hg38.fa` and the K562 host-expression table. The scripts are here for
  reading; a reference bundle may be added later.
