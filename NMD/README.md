# NMD

Does a TRIP integration turn the host transcript into a target of nonsense-mediated decay? The
pipeline predicts it for every barcode from where the reporter landed and what it inserts into the
mRNA. The prediction is then held against the measured response to cycloheximide (CHX), which
blocks NMD.

| | |
|---|---|
| `HOW_IT_WORKS.md` | the pipeline: output columns, decision tree, stages, toy test, how to run |
| `scripts/` | the stages, Python 3 standard library only. 01–04 for every construct, 05–06 for IntronTRIP only |
| `toy_test/` | an invented miniature genome run through the production scripts |
| `single_combined/` | **5ss_noP and CPA_noP**, the single-exon cassette in two libraries, against CHX. `SUMMARY.md` with the figures, figure script in `analysis/` |
| `IntronTRIP/` | **the eight IntronTRIP reporters** (A1, B1, B3-6, B3-9, C1, D1, E1, H1), a cassette with its own intron, spliced and unspliced isoforms against CHX. `SUMMARY.md` with the figures, `pics/`, `make_figures.py` |

The input and output tables of both datasets are not in the repository, only the figures.
