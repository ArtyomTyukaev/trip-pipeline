# TRIP

**TRIP** (Thousands of Reporters Integrated in Parallel) is a barcoded-reporter assay. A reporter
construct is integrated at thousands of random positions in the K562 genome with **piggyBac**
transposase, and every integration carries a unique **18 bp barcode**. Sequencing then links each
barcode to its **genomic position** (from DNA "Map" libraries) and to its **RNA readout** — expression,
splicing, polyadenylation. Because the reporter sequence is constant, differences in how it behaves
are attributed to the **local genomic context**.

This repository holds two things.

## [`general_pipeline/`](general_pipeline/) — raw reads to per-barcode tables

The primary processing chain: trimmed sequencing reads in, the flagship annotated per-barcode table
out. Ten numbered stages, plus a small real slice of data the whole thing can be run on end to end.

| | |
|---|---|
| `scripts/` | stages 01–10, one current version of each |
| `example_data/` | a small real dataset and `run_example.sh` — stages 03→07 in about two minutes |
| `FILE_FORMATS.md` | every intermediate file format, column by column |
| `LAB_NOTEBOOK.md` | dated development log; `LAB_NOTEBOOK_EXTENDED.md` is the longer version |
| `figures/` | figures referenced by the notebook |

## [`NMD/`](NMD/) — does an integration make the host transcript an NMD target?

A downstream analysis built on top of that table. For every barcode it predicts, from coordinates
alone, whether the reporter creates a premature stop codon and therefore turns the host transcript
into a substrate for **nonsense-mediated decay** — then holds that prediction against the measured
response to cycloheximide, which switches NMD off.

| | |
|---|---|
| `HOW_IT_WORKS.md` | the pipeline: every column of the output table, the decision tree, the stages, the toy test |
| `scripts/` | the stages, Python 3 standard library only. 01–04 for every construct, 05–06 for IntronTRIP only |
| `toy_test/` | an invented miniature genome with its inputs and outputs, run through the production scripts and compared against a hand-written expectation |
| `single_combined/` | 5ss_noP and CPA_noP against CHX: `SUMMARY.md` with the figures, `analysis/` with the figure script |
| `IntronTRIP/` | the eight IntronTRIP reporters, spliced and unspliced isoforms, against CHX: `SUMMARY.md` with the figures, `pics/`, `make_figures.py` |

Start with `single_combined/SUMMARY.md` or `IntronTRIP/SUMMARY.md` for what was found, or
`HOW_IT_WORKS.md` for how it was computed.
