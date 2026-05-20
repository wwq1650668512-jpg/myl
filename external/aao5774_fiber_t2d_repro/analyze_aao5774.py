#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from scipy.stats import mannwhitneyu, spearmanr
except Exception:  # pragma: no cover
    mannwhitneyu = None
    spearmanr = None


def _read_abundance(path: Path, sep: str) -> pd.DataFrame:
    frame = pd.read_csv(path, sep=sep)
    if "sample_id" not in frame.columns:
        raise ValueError("Abundance file must contain a 'sample_id' column.")
    frame = frame.copy()
    frame["sample_id"] = frame["sample_id"].astype(str)
    frame = frame.set_index("sample_id")
    frame = frame.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return frame


def _read_metadata(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = {"sample_id", "group", "hba1c_change"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"Metadata is missing required columns: {', '.join(missing)}")
    frame = frame.copy()
    frame["sample_id"] = frame["sample_id"].astype(str)
    frame["hba1c_change"] = pd.to_numeric(frame["hba1c_change"], errors="coerce")
    return frame


def _read_taxa_list(path: Path) -> list[str]:
    taxa: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        name = raw.strip()
        if not name or name.startswith("#"):
            continue
        taxa.append(name)
    return taxa


def _group_stats(values: pd.Series) -> dict[str, float]:
    numeric = pd.to_numeric(values, errors="coerce").dropna()
    if numeric.empty:
        return {"n": 0, "mean": np.nan, "median": np.nan, "std": np.nan}
    return {
        "n": int(numeric.shape[0]),
        "mean": float(numeric.mean()),
        "median": float(numeric.median()),
        "std": float(numeric.std(ddof=1)) if numeric.shape[0] > 1 else 0.0,
    }


def _safe_spearman(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    merged = pd.concat([x, y], axis=1).dropna()
    if merged.shape[0] < 3:
        return np.nan, np.nan
    if merged.iloc[:, 0].nunique(dropna=True) < 2 or merged.iloc[:, 1].nunique(dropna=True) < 2:
        return np.nan, np.nan
    if spearmanr is not None:
        coef, pvalue = spearmanr(merged.iloc[:, 0], merged.iloc[:, 1])
        return float(coef), float(pvalue)
    coef = merged.iloc[:, 0].corr(merged.iloc[:, 1], method="spearman")
    return float(coef), np.nan


def _safe_mannwhitney(a: pd.Series, b: pd.Series) -> tuple[float, float]:
    a_num = pd.to_numeric(a, errors="coerce").dropna()
    b_num = pd.to_numeric(b, errors="coerce").dropna()
    if a_num.empty or b_num.empty:
        return np.nan, np.nan
    if mannwhitneyu is None:
        return np.nan, np.nan
    stat, pvalue = mannwhitneyu(a_num, b_num, alternative="two-sided")
    return float(stat), float(pvalue)


def _plot_if_possible(data: pd.DataFrame, output_dir: Path, treatment: str, control: str) -> list[str]:
    created: list[str] = []
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return created

    out_box = output_dir / "scfa_relative_abundance_by_group.png"
    sub = data[["group", "scfa_relative_abundance"]].dropna()
    group_order = [control, treatment]
    box_data = [sub.loc[sub["group"] == g, "scfa_relative_abundance"].values for g in group_order]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.boxplot(box_data, labels=group_order)
    ax.set_ylabel("SCFA-producer relative abundance")
    ax.set_title("SCFA abundance by group")
    fig.tight_layout()
    fig.savefig(out_box, dpi=180)
    plt.close(fig)
    created.append(str(out_box))

    out_scatter = output_dir / "scfa_abundance_vs_hba1c_change.png"
    pair = data[["scfa_relative_abundance", "hba1c_change"]].dropna()
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(pair["scfa_relative_abundance"], pair["hba1c_change"], alpha=0.8)
    ax.set_xlabel("SCFA-producer relative abundance")
    ax.set_ylabel("HbA1c change")
    ax.set_title("SCFA abundance vs HbA1c change")
    fig.tight_layout()
    fig.savefig(out_scatter, dpi=180)
    plt.close(fig)
    created.append(str(out_scatter))
    return created


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Runnable analysis scaffold inspired by Science 2018 (10.1126/science.aao5774): "
            "group-wise SCFA-producer enrichment and HbA1c association."
        )
    )
    parser.add_argument("--abundance", type=Path, required=True, help="Sample-by-taxon abundance table (TSV/CSV).")
    parser.add_argument("--metadata", type=Path, required=True, help="Metadata CSV with sample_id/group/hba1c_change.")
    parser.add_argument("--scfa-list", type=Path, required=True, help="Text file of SCFA-producing taxa, one per line.")
    parser.add_argument(
        "--detrimental-list",
        type=Path,
        default=None,
        help="Optional text file of potentially detrimental taxa (for additional metrics).",
    )
    parser.add_argument("--output-dir", type=Path, required=True, help="Output directory.")
    parser.add_argument("--sep", default="\t", help="Delimiter for abundance table (default: tab).")
    parser.add_argument("--treatment-label", default="high_fiber", help="Treatment group label in metadata.group.")
    parser.add_argument("--control-label", default="control", help="Control group label in metadata.group.")
    parser.add_argument("--presence-threshold", type=float, default=1e-4, help="Taxon presence threshold for richness.")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    abundance = _read_abundance(args.abundance, sep=args.sep)
    metadata = _read_metadata(args.metadata)
    scfa_taxa = _read_taxa_list(args.scfa_list)
    detrimental_taxa = _read_taxa_list(args.detrimental_list) if args.detrimental_list else []

    shared = sorted(set(metadata["sample_id"]).intersection(abundance.index.astype(str)))
    if not shared:
        raise ValueError("No overlapping sample_id found between abundance and metadata.")

    abundance = abundance.loc[shared].copy()
    meta = metadata.set_index("sample_id").loc[shared].copy()

    scfa_present = [tax for tax in scfa_taxa if tax in abundance.columns]
    detrimental_present = [tax for tax in detrimental_taxa if tax in abundance.columns]

    if not scfa_present:
        raise ValueError(
            "None of the SCFA taxa were found in abundance columns. "
            "Check taxonomy naming consistency."
        )

    scfa_matrix = abundance[scfa_present]
    sample_metrics = pd.DataFrame(index=abundance.index)
    sample_metrics["group"] = meta["group"]
    sample_metrics["hba1c_change"] = meta["hba1c_change"]
    sample_metrics["scfa_relative_abundance"] = scfa_matrix.sum(axis=1)
    sample_metrics["scfa_richness"] = (scfa_matrix > args.presence_threshold).sum(axis=1)

    if detrimental_present:
        sample_metrics["detrimental_relative_abundance"] = abundance[detrimental_present].sum(axis=1)
    else:
        sample_metrics["detrimental_relative_abundance"] = np.nan

    sample_metrics = sample_metrics.reset_index().rename(columns={"index": "sample_id"})
    sample_metrics.to_csv(args.output_dir / "sample_metrics.csv", index=False)

    comparisons: list[dict[str, object]] = []
    for metric in ["scfa_relative_abundance", "scfa_richness", "detrimental_relative_abundance"]:
        treat_vals = sample_metrics.loc[sample_metrics["group"] == args.treatment_label, metric]
        ctrl_vals = sample_metrics.loc[sample_metrics["group"] == args.control_label, metric]
        stat, pvalue = _safe_mannwhitney(treat_vals, ctrl_vals)
        treat_stat = _group_stats(treat_vals)
        ctrl_stat = _group_stats(ctrl_vals)
        comparisons.append(
            {
                "metric": metric,
                "treatment_label": args.treatment_label,
                "control_label": args.control_label,
                "treatment_n": treat_stat["n"],
                "control_n": ctrl_stat["n"],
                "treatment_mean": treat_stat["mean"],
                "control_mean": ctrl_stat["mean"],
                "treatment_median": treat_stat["median"],
                "control_median": ctrl_stat["median"],
                "mannwhitney_u": stat,
                "p_value": pvalue,
            }
        )

    pd.DataFrame(comparisons).to_csv(args.output_dir / "group_comparison.csv", index=False)

    corr_rows: list[dict[str, object]] = []
    for metric in ["scfa_relative_abundance", "scfa_richness", "detrimental_relative_abundance"]:
        coef, pvalue = _safe_spearman(sample_metrics[metric], sample_metrics["hba1c_change"])
        corr_rows.append(
            {
                "metric": metric,
                "target": "hba1c_change",
                "spearman_rho": coef,
                "p_value": pvalue,
            }
        )
    pd.DataFrame(corr_rows).to_csv(args.output_dir / "hba1c_correlation.csv", index=False)

    created_plots = _plot_if_possible(sample_metrics, args.output_dir, args.treatment_label, args.control_label)

    summary = {
        "n_samples": int(sample_metrics.shape[0]),
        "n_taxa": int(abundance.shape[1]),
        "scfa_taxa_requested": len(scfa_taxa),
        "scfa_taxa_found": len(scfa_present),
        "detrimental_taxa_requested": len(detrimental_taxa),
        "detrimental_taxa_found": len(detrimental_present),
        "treatment_label": args.treatment_label,
        "control_label": args.control_label,
        "plots_created": created_plots,
        "note": (
            "This is a runnable reproduction scaffold inspired by 10.1126/science.aao5774. "
            "It is not an official author-released pipeline."
        ),
    }
    (args.output_dir / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps({"status": "ok", "output_dir": str(args.output_dir), "n_samples": summary["n_samples"]}, indent=2))


if __name__ == "__main__":
    main()
