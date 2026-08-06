#!/bin/bash
# End-to-end toy test: build the toy genome, run the PRODUCTION scripts over it, compare against
# the hand calculation. Nothing here touches the real data.
set -euo pipefail
cd "$(dirname "$0")/.."

python3 toy_test/make_toy.py

for TAG in A B; do
    ROOT="toy_test/root_${TAG}"
    echo
    echo "================ variant ${TAG} ================"
    python3 scripts/01_build_transcripts.py --root "$ROOT"
    python3 scripts/02_map_barcodes.py      --root "$ROOT" --libs 5ss_noP
    python3 scripts/03_call_nmd.py          --root "$ROOT" --libs 5ss_noP
    python3 scripts/04_make_table.py        --root "$ROOT" --libs 5ss_noP
done

echo
echo "================ comparison ================"
python3 toy_test/check.py
