#!/usr/bin/env bash
# Execute the notebooks in order with the project environment.
# Executed copies (with outputs) go to products/executed/; the notebooks/ sources stay clean.
#   ./run_all.sh            -> all notebooks
#   ./run_all.sh 04 05      -> only the ones whose name starts with these prefixes
set -euo pipefail
cd "$(dirname "$0")"
PY="${PAPER_I_PYTHON:-$HOME/.venvs/solar-corona/bin/python}"
mkdir -p products/executed
if [ $# -eq 0 ]; then set -- 00 01 02 03 04 05 06 07 08 09 10 99; fi
for p in "$@"; do
  for nb in notebooks/${p}_*.ipynb; do
    echo "=== $(date '+%H:%M:%S') running $nb"
    "$PY" -m jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=-1 \
        --output-dir products/executed --output "$(basename "$nb")" "$nb"
  done
done
echo "=== $(date '+%H:%M:%S') done"
