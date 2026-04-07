#!/usr/bin/env bash
set -euo pipefail

ENV_DIR="${1:-/tmp/chemprop_env}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

echo "[step1-chemprop] creating environment at ${ENV_DIR}"
"${PYTHON_BIN}" -m venv "${ENV_DIR}"

echo "[step1-chemprop] upgrading pip"
"${ENV_DIR}/bin/python" -m pip install --upgrade pip

echo "[step1-chemprop] installing CPU-only PyTorch from the official index"
"${ENV_DIR}/bin/python" -m pip install torch --index-url https://download.pytorch.org/whl/cpu

echo "[step1-chemprop] installing Chemprop"
"${ENV_DIR}/bin/python" -m pip install chemprop

echo "[step1-chemprop] done"
echo "[step1-chemprop] python: ${ENV_DIR}/bin/python"
