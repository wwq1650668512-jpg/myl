#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${SCRIPT_DIR}/output/demo_run"
mkdir -p "${OUT_DIR}"

if [[ -n "${PYTHON_BIN:-}" ]]; then
  PY="${PYTHON_BIN}"
elif [[ -x "/tmp/microbe_env/bin/python" ]]; then
  PY="/tmp/microbe_env/bin/python"
else
  PY="python3"
fi

"${PY}" "${SCRIPT_DIR}/analyze_aao5774.py" \
  --abundance "${SCRIPT_DIR}/demo/abundance.tsv" \
  --metadata "${SCRIPT_DIR}/demo/metadata.csv" \
  --scfa-list "${SCRIPT_DIR}/demo/scfa_producers.txt" \
  --detrimental-list "${SCRIPT_DIR}/demo/detrimental_taxa.txt" \
  --output-dir "${OUT_DIR}" \
  --sep $'\t' \
  --treatment-label "high_fiber" \
  --control-label "control"

printf '\nDone. Results in: %s\n' "${OUT_DIR}"
