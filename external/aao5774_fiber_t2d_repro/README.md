# AAO5774 Runnable Analysis Scaffold

This folder contains a **runnable analysis scaffold** inspired by:

- Zhao et al., *Science* (2018)
- "Gut bacteria selectively promoted by dietary fibers alleviate type 2 diabetes"
- DOI: `10.1126/science.aao5774`

## Important note

No official author-released GitHub pipeline was identified from accessible public sources during setup.
This folder is therefore a practical, reusable workflow to run the core analysis logic on your own data.

## What this code does

Given:

- sample-by-taxon abundance table
- sample metadata (including group and HbA1c change)
- SCFA-producer taxa list

it computes:

- sample-level SCFA metrics
- group comparison (`high_fiber` vs `control`) via Mann-Whitney U
- correlation with HbA1c change (Spearman)
- optional figures (if `matplotlib` is installed)

## Folder layout

- `analyze_aao5774.py`: main CLI script
- `run_demo.sh`: one-click demo run
- `demo/`: small synthetic example input files
- `output/`: generated outputs
- `requirements.txt`: dependencies

## Quick start

From this folder:

```bash
python3 -m pip install -r requirements.txt
./run_demo.sh
```

`run_demo.sh` auto-detects Python in this order:

1. `PYTHON_BIN` (if you set it)
2. `/tmp/microbe_env/bin/python` (if present)
3. `python3`

You can force a specific interpreter:

```bash
PYTHON_BIN=/path/to/python ./run_demo.sh
```

Expected outputs under `output/demo_run/`:

- `sample_metrics.csv`
- `group_comparison.csv`
- `hba1c_correlation.csv`
- `run_summary.json`
- PNG plots (when plotting deps are available)

## Run on your own dataset

```bash
python3 analyze_aao5774.py \
  --abundance /path/to/abundance.tsv \
  --metadata /path/to/metadata.csv \
  --scfa-list /path/to/scfa_producers.txt \
  --detrimental-list /path/to/detrimental_taxa.txt \
  --output-dir /path/to/output \
  --sep $'\t' \
  --treatment-label high_fiber \
  --control-label control
```

## Input requirements

### `abundance` table

- Must contain a `sample_id` column.
- Other columns are taxa names (must match names in SCFA taxa list).

### `metadata` table

Must contain columns:

- `sample_id`
- `group`
- `hba1c_change`

### taxa list files

- One taxon name per line.
- Lines starting with `#` are ignored.
