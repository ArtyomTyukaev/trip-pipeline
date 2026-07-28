#!/bin/bash
# =============================================================================
# Run stages 03 -> 07 of the TRIP pipeline on the small example dataset.
# =============================================================================
# Takes about 2 minutes, needs no cluster, no genome and no annotation files.
# Everything is written to example_data/out/ ; the committed inputs stay untouched.
#
#   bash example_data/run_example.sh
#
# Requires R (>= 4.3) with DescTools, dplyr, purrr, tidyr, ggplot2, cowplot, gridExtra.
# On the Weizmann cluster: module load R/4.3.1-foss-2022b (done automatically below).
# =============================================================================
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SCRIPTS="$HERE/../scripts"
OUT="$HERE/out"

MAP_5SS=Map_5ss_noP_pool1_rep2_S166
MAP_CPA=Map_CPA_noP_pool1_rep2_S164
NORM=Norm_5ss_noP_rep2_S92
EXP=Exp_total_5ss_noP_pool2_rep2_S80

command -v module >/dev/null 2>&1 && module load R/4.3.1-foss-2022b 2>/dev/null || true

mkdir -p "$OUT"
cp -f "$HERE"/*.gz "$HERE"/fileList.txt "$HERE"/5ss_noP.txt "$HERE"/combos.txt "$OUT"/
cd "$OUT"

echo "== Stage 03: barcodes + alignments -> integration positions =="
for LIB in "$MAP_5SS" "$MAP_CPA"; do
  Rscript "$SCRIPTS/03_analyzeTRIPMapping_uniform.R" ./ ./ "$LIB" > "$LIB.stage03.log" 2>&1
  echo "   $LIB: $(wc -l < "$LIB.barcodes_poses.txt") positions"

  # ---------------------------------------------------------------------------
  # Version adapter (see README, "Known version mismatch"). This copy of stage 03
  # writes a 7th column (alignment Length) that the production *.barcodes_poses.txt
  # files do not have. Stage 05 keys on columns 1,2,9 (chrom, start, barcodes) and
  # enforces an 11-column input, so the extra column would both shift that key onto
  # cleanCounts and abort the run. Drop it here to match the production layout.
  # ---------------------------------------------------------------------------
  cut -f1-6,8-10 "$LIB.barcodes_poses.txt" > "$LIB.poses.9col" && mv "$LIB.poses.9col" "$LIB.barcodes_poses.txt"
done

echo "== Stage 04: drop near-duplicate positions (TTAA filter left OFF) =="
# To switch the TTAA filter on, add: ttaaFilter /path/to/hg38_ttaa.bed
# (that reference is 627 MB, so it is not part of this repository.)
for LIB in "$MAP_5SS" "$MAP_CPA"; do
  Rscript "$SCRIPTS/04_filterMappings_TTAA.R" ./ "$LIB" > "$LIB.stage04.log" 2>&1
  echo "   $LIB: $(wc -l < "$LIB.filtered_barcodes_poses.txt") kept"
done

echo "== Stage 05: cross-library contaminant filter (2 groups: 5ss_noP vs CPA_noP) =="
Rscript "$SCRIPTS/05_filterContaminants_uniform.R" ./ fileList.txt > stage05.log 2>&1
for LIB in "$MAP_5SS" "$MAP_CPA"; do
  echo "   $LIB: $(wc -l < "$LIB.retained_barcodes.bed") retained"
done

echo "== Stage 06: genuine (single-position) barcodes =="
Rscript "$SCRIPTS/06_get_genuine_barcodes_v2.R" ./ 5ss_noP.txt 18 > stage06.log 2>&1
echo "   genuine:      $(wc -l < 5ss_noP.txt.retained_barcodes_genuine.bed)"
echo "   multimappers: $(wc -l < 5ss_noP.txt.retained_barcodes_multimappers.bed)"

echo "== Stage 07: join Map + Norm + Exp per barcode =="
while IFS=$'\t' read -r MAP NORM_N EXP_N LEN RC; do
  Rscript "$SCRIPTS/07_normalizeTRIP_v2_MM.R" ./ ./ "$MAP" "$NORM_N" "$EXP_N" "$LEN" "$RC" \
    > stage07_exact.log 2>&1
  cp -f "$EXP_N.combined_counts.txt" "$EXP_N.combined_counts.exact.txt"

  # same call plus the optional 8th argument -> allow exactly one substitution
  Rscript "$SCRIPTS/07_normalizeTRIP_v2_MM.R" ./ ./ "$MAP" "$NORM_N" "$EXP_N" "$LEN" "$RC" oneMismatch \
    > stage07_1mm.log 2>&1
  cp -f "$EXP_N.combined_counts.txt" "$EXP_N.combined_counts.1mm.txt"
done < combos.txt

echo
echo "== Result =="
printf '   exact match     : %s barcodes\n' "$(( $(wc -l < "$EXP.combined_counts.exact.txt") - 1 ))"
printf '   oneMismatch flag: %s barcodes\n' "$(( $(wc -l < "$EXP.combined_counts.1mm.txt") - 1 ))"
awk 'NR>1{n+=$2;e+=$3} END{printf "   exact  : norm=%d exp=%d reads\n", n, e}' "$EXP.combined_counts.exact.txt"
awk 'NR>1{n+=$2;e+=$3} END{printf "   1mm    : norm=%d exp=%d reads\n", n, e}' "$EXP.combined_counts.1mm.txt"
echo
echo "All outputs and per-stage logs are in $OUT"
