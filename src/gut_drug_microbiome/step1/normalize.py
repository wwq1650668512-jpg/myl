# 清洗、标准化、整合原始实验数据，生成可用于机器学习建模的药物 - 微生物互作数据集
# 1. 读取原始数据文件，提取药物、微生物、互作等相关信息
# 2. 标准化字段名称，清洗数据格式，处理缺失值
# 3. 根据预设的阈值对互作效果进行打分和分类，生成连续的效果分数和离散的标签
# 4. 将处理后的数据表保存为CSV文件，并生成一个包含统计摘要和输出文件路径的JSON报告

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from .chem_features import enrich_drug_table_with_rdkit
from gut_drug_microbiome.utils.chem import compute_smiles_descriptors as _compute_smiles_descriptors

# 列名标准化函数，将原始列标签转换为稳定的snake_case字段名，去除特殊字符，统一大小写，并处理空值
def _snake_case(text: str) -> str:
    """Normalize a raw column label into a stable snake_case field name."""
    value = str(text).strip().lower()
    replacements = {
        "µ": "u",
        "μ": "u",
        "Å": "a",
        "å": "a",
        "%": "pct",
        "#": "num",
        "/": "_",
        "-": "_",
        ".": "_",
    }
    for source, target in replacements.items():
        value = value.replace(source, target)
    # 非字母数字 → 下划线，多下划线合并，首尾去下划线
    value = re.sub(r"[^a-z0-9]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("_")
    return value or "unnamed"

# 从YAML配置文件中读取Step 1的标签阈值，返回一个包含数值型cutoff的字典
def _read_label_thresholds(config_path: Path) -> dict[str, float]:
    """Load Step 1 labeling thresholds from YAML and return numeric cutoffs."""
    # YAML文件中可能包含字符串形式的数值，使用float()转换为数值类型，并取绝对值以确保正确的比较逻辑
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    step1 = config["step1"]["default_thresholds"]
    return {
        "inhibit_threshold": abs(float(step1["inhibit"]["effect_score_lte"])),
        "promote_threshold": float(step1["promote"]["effect_score_gte"]),
        "inhibit_q_threshold": float(step1["inhibit"]["significance_q_value_lt"]),
        "promote_q_threshold": float(step1["promote"]["significance_q_value_lt"]),
    }

# 下面是一些用于数据聚合和特征计算的辅助函数：
# _collapse_codes：将多个代码值合并为一个以管道分隔的字符串，如果没有有效值则返回NaN
def _collapse_codes(values: pd.Series) -> str | float:
    """Merge multiple code values into one pipe-delimited string, or NaN if empty."""
    cleaned = sorted({str(value).strip() for value in values.dropna() if str(value).strip() and str(value).strip() != "-"})
    if not cleaned:
        return np.nan
    return "|".join(cleaned)

# _first_or_nan：从一系列值中返回第一个非空值，如果没有则返回NaN
def _first_or_nan(values: pd.Series) -> str | float:
    """Return the first non-empty value from a series, otherwise NaN."""
    cleaned = [value for value in values.dropna() if str(value).strip()]
    return cleaned[0] if cleaned else np.nan

# _extract_primary_code：从一个以管道分隔的代码字符串中提取第一个排序的代码，并根据请求的层级宽度进行截断
def _extract_primary_code(code_string: str | float, width: int) -> str | float:
    """Pick the first sorted code and truncate it to the requested hierarchy width."""
    if not isinstance(code_string, str) or not code_string:
        return np.nan
    codes = [code for code in code_string.split("|") if code]
    if not codes:
        return np.nan
    primary = sorted(codes)[0]
    return primary[:width]

# _safe_log10：计算-log10值，先将输入值限制在一个非常小的正数以上，以避免对零或负数取对数时出现无穷大或NaN
def _safe_log10(series: pd.Series) -> pd.Series:
    """Compute -log10 safely by clipping tiny values away from zero first."""
    clipped = series.clip(lower=1e-300)
    return -np.log10(clipped)

# _load_supplementary_drug_metadata：从Maier 2018的补充表格中读取药物元数据，并标准化关键列名以便后续合并和分析
def _load_supplementary_drug_metadata(raw_dir: Path) -> pd.DataFrame:
    """Read Maier supplementary drug metadata and standardize key column names."""
    path = raw_dir / "Supplementary_table_1.xlsx"
    frame = pd.read_excel(path, sheet_name="S1a. Prestwick_Libery")
    frame = frame.rename(columns={column: _snake_case(column) for column in frame.columns})
    frame = frame.rename(
        columns={
            "chemical_name": "chemical_name_supp",
            "stitch4_id": "cid_flat",
            "target_species": "target_species_library",
            "dose_umol": "dose_umol",
            "estimated_intestine_concentration_um": "estimated_intestine_concentration_um",
            "plasma_concentration_um": "plasma_concentration_um",
            "source_for_plasma_concentration": "source_for_plasma_concentration",
            "fraction_excreted_in_feces": "fraction_excreted_in_feces",
            "fraction_excreted_in_urine": "fraction_excreted_in_urine",
            "source_for_excretion_data": "source_for_excretion_data",
            "estimated_colon_concentration_um": "estimated_colon_concentration_um",
            "molecular_weight_g_mol": "molecular_weight_supp",
            "tpsa_a": "tpsa_supp",
            "complexity": "complexity_supp",
            "volume3d_a": "volume3d_supp",
            "screen_conc_20_um_as_ug_ml": "screen_conc_20_um_as_ug_ml",
        }
    )
    return frame

# _load_drug_table：构建一个以药物为行、合并了不同来源的元数据和化学特征的表格，包含药物的基本信息、化学性质、治疗分类、靶标物种等字段，并计算一些基于SMILES的文本描述符
def _load_drug_table(raw_dir: Path) -> pd.DataFrame:
    """Build a one-row-per-drug table with merged metadata and chemistry features."""
    compound = pd.read_csv(raw_dir / "compound_properties.tsv", sep="\t")
    compound = compound.rename(columns={column: _snake_case(column) for column in compound.columns})
    compound = compound.rename(
        columns={
            "molecularformula": "molecular_formula",
            "molecularweight": "molecular_weight",
        }
    )

    atc = pd.read_csv(raw_dir / "prestwick_atc.tsv", sep="\t")
    atc = atc.rename(columns={column: _snake_case(column) for column in atc.columns})
    atc = atc.groupby("prestwick_id", as_index=False).agg(
        {
            "chemical_name": _first_or_nan,
            "chemical_formula": _first_or_nan,
            "molecular_weight": "first",
            "therapeutic_class": _first_or_nan,
            "therapeutic_effect": _first_or_nan,
            "lc_chemical_name": _first_or_nan,
            "cid_flat": _first_or_nan,
            "atc": _collapse_codes,
        }
    )
    atc = atc.rename(columns={"atc": "atc_codes"})
    atc["atc_primary_l1"] = atc["atc_codes"].apply(lambda value: _extract_primary_code(value, 1))
    atc["atc_primary_l3"] = atc["atc_codes"].apply(lambda value: _extract_primary_code(value, 3))
    atc["atc_primary_l4"] = atc["atc_codes"].apply(lambda value: _extract_primary_code(value, 4))

    supplement = _load_supplementary_drug_metadata(raw_dir)

    pair_metadata = pd.read_csv(raw_dir / "combined_pv.tsv", sep="\t")
    pair_metadata = pair_metadata.rename(columns={column: _snake_case(column) for column in pair_metadata.columns})
    pair_metadata = pair_metadata[
        ["prestwick_id", "target_species", "veterinary", "human_use"]
    ].drop_duplicates()

    drugs = (
        compound.merge(atc, on="prestwick_id", how="outer", suffixes=("", "_atc"))
        .merge(supplement, on="prestwick_id", how="left", suffixes=("", "_supp"))
        .merge(pair_metadata, on="prestwick_id", how="left", suffixes=("", "_screen"))
    )

    drugs["chemical_name"] = drugs["chemical_name"].fillna(drugs["chemical_name_supp"])
    drugs["target_species"] = drugs["target_species"].fillna(drugs["target_species_library"])
    drugs["molecular_weight"] = drugs["molecular_weight"].fillna(drugs["molecular_weight_supp"])
    drugs["tpsa"] = drugs["tpsa"].fillna(drugs["tpsa_supp"])
    drugs["complexity"] = drugs["complexity"].fillna(drugs["complexity_supp"])
    drugs["volume3d"] = drugs["volume3d"].fillna(drugs["volume3d_supp"])
    # 以下行列分别是药物ID、化学名称、CID编号、SMILES字符串、分子式、分子量、拓扑极表面积、复杂度、3D体积、治疗分类、治疗效果、ATC代码、靶标物种、人用/兽用标记等关键信息，后续会进行去重和特征计算
    columns_to_keep = [
        "prestwick_id",
        "chemical_name",
        "cid_flat",
        "cid_active",
        "cid_main",
        "main_component_smiles",
        "smiles",
        "molecular_formula",
        "molecular_weight",
        "xlogp",
        "tpsa",
        "complexity",
        "volume3d",
        "therapeutic_class",
        "therapeutic_effect",
        "atc_codes",
        "atc_primary_l1",
        "atc_primary_l3",
        "atc_primary_l4",
        "target_species",
        "human_use",
        "veterinary",
        "dose_umol",
        "estimated_intestine_concentration_um",
        "plasma_concentration_um",
        "fraction_excreted_in_feces",
        "fraction_excreted_in_urine",
        "estimated_colon_concentration_um",
        "screen_conc_20_um_as_ug_ml",
    ]
    # 通过prestwick_id去重，确保每个药物只有一行记录，并计算基于SMILES的文本描述符，然后使用RDKit计算更丰富的化学特征，最终返回一个包含所有关键信息和特征的药物表格
    drugs = drugs[columns_to_keep].drop_duplicates(subset=["prestwick_id"]).reset_index(drop=True)
    drugs = _compute_smiles_descriptors(drugs)
    drugs = enrich_drug_table_with_rdkit(drugs, smiles_columns=["main_component_smiles", "smiles"])
    return drugs

# _load_microbe_table：构建一个以微生物为行、包含分类学信息和实验上下文字段的表格，合并了物种概览和筛选选择表中的相关信息，并生成用于标签显示的字段
def _load_microbe_table(raw_dir: Path) -> pd.DataFrame:
    """Build a one-row-per-microbe table with taxonomy and assay context fields."""
    overview = pd.read_csv(raw_dir / "species_overview.tsv", sep="\t")
    overview = overview.rename(columns={column: _snake_case(column) for column in overview.columns})

    selection = pd.read_excel(raw_dir / "Supplementary_table_2.xlsx", sheet_name="S2. Species selection", header=None)
    header = [_snake_case(value) for value in selection.iloc[2].tolist()]
    selection = selection.iloc[3:].copy()
    selection.columns = header
    selection = selection.rename(columns={"nt_data_base": "nt_code"})
    selection = selection[selection["nt_code"].astype(str).str.startswith("NT")].copy()

    microbes = overview.merge(selection, on="nt_code", how="left", suffixes=("_overview", ""))
    microbes["species_label"] = microbes["species_label"].fillna(microbes["species"])
    microbes["microbe_label"] = microbes["species_label"]
    # 以下字段分别是：nt码、微生物标签、物种标签、物种名称、物种名称（简化）、菌株信息、物种聚类、是否唯一菌株、生物安全等级、门、纲、目、科、属、革兰氏染色结果、培养基偏好、96孔板筛选起始OD值、384孔板筛选起始OD值等，这些字段将用于后续的分析和建模
    columns_to_keep = [
        "nt_code",
        "microbe_label",
        "species_label",
        "species_name",
        "species",
        "strain",
        "speci_cluster",
        "is_unique",
        "biosafety",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "gram_stain",
        "medium_preference",
        "starting_od_96_well_screen",
        "starting_od_384_well_screen",
    ]
    microbes = microbes[columns_to_keep].drop_duplicates(subset=["nt_code"]).reset_index(drop=True)
    return microbes

# _label_effects：根据预设的阈值将实验统计结果转换为连续的效果分数和离散的Step 1标签，计算效果分数、显著性标记、候选促进标签等，并生成最终的标签字段
def _label_effects(frame: pd.DataFrame, thresholds: dict[str, float]) -> pd.DataFrame:
    """Convert assay statistics into continuous effect scores and discrete Step 1 labels."""
    inhibit_q_threshold = thresholds["inhibit_q_threshold"]
    promote_q_threshold = thresholds["promote_q_threshold"]
    inhibit_threshold = thresholds["inhibit_threshold"]
    promote_threshold = thresholds["promote_threshold"]

    labeled = frame.copy()
    # effect_score：基于AUC值计算的连续效果分数，表示相对于中性效应（AUC=1.0）的增强或抑制程度，正值表示促进，负值表示抑制
    labeled["effect_score"] = labeled["auc"] - 1.0
    # effect_score_raw：基于原始AUC值计算的连续效果分数，作为一个未经调整的参考指标，可能与最终的effect_score存在一定的相关性但不完全相同
    labeled["effect_score_raw"] = labeled["mean_norm_auc_raw"] - 1.0
    # effect_score_pct：将effect_score转换为百分比形式，表示相对于中性效应的百分比变化，便于理解和比较不同互作的效果大小
    labeled["effect_score_pct"] = labeled["effect_score"] * 100.0
    # neg_log10_q：基于FDR校正后的p值计算的-log10值，作为一个连续的显著性指标，数值越大表示结果越显著，通常用于可视化和筛选显著互作
    labeled["neg_log10_q"] = _safe_log10(labeled["pv_comb_fdr_bh"])
    # is_significant_inhibit / is_significant_promote：分别记录抑制和促进任务的显著性门槛，允许正负方向使用不同q阈值
    labeled["is_significant_inhibit"] = labeled["pv_comb_fdr_bh"] < inhibit_q_threshold
    labeled["is_significant_promote"] = labeled["pv_comb_fdr_bh"] < promote_q_threshold
    labeled["is_significant"] = labeled["is_significant_inhibit"] | labeled["is_significant_promote"]
    # candidate_promote_no_fdr：一个布尔字段，表示该互作是否具有促进效果且不考虑FDR校正的显著性（即effect_score大于等于promote_threshold），作为一个宽松的候选促进标签，可能包含一些假阳性但有助于捕获潜在的促进互作
    labeled["candidate_promote_no_fdr"] = labeled["effect_score"] >= promote_threshold

    labels = np.full(len(labeled), "no_effect", dtype=object)
    inhibit_mask = labeled["is_significant_inhibit"] & (labeled["effect_score"] <= -inhibit_threshold)
    promote_mask = labeled["is_significant_promote"] & (labeled["effect_score"] >= promote_threshold)
    labels[inhibit_mask] = "inhibit"
    labels[promote_mask] = "promote"
    labeled["effect_label"] = labels
    labeled["binary_effect_label"] = np.where(labeled["effect_label"] == "inhibit", "inhibit", "no_effect")
    labeled["source_dataset"] = "maier_2018"
    labeled["label_tier"] = "gold"
    return labeled

# _load_interaction_table：组装药物-微生物互作表格，应用过滤条件和标签生成逻辑，返回一个包含互作统计结果、标签和相关元数据的DataFrame
def _load_interaction_table(
    raw_dir: Path,
    human_use_only: bool,
    primary_panel_only: bool,
    primary_codes: set[str],
    thresholds: dict[str, float],
) -> pd.DataFrame:
    """Assemble the drug-microbe interaction table and apply filtering plus labels."""
    pv = pd.read_csv(raw_dir / "combined_pv.tsv", sep="\t")
    pv = pv.rename(columns={column: _snake_case(column) for column in pv.columns})

    aucs = pd.read_csv(raw_dir / "aucs.tsv", sep="\t")
    aucs = aucs.rename(columns={column: _snake_case(column) for column in aucs.columns})
    aucs = aucs[aucs["prestwick_id"] != "CONTROL"].copy()

    pairs = pv[["nt_code", "prestwick_id"]].drop_duplicates()
    aucs = aucs.merge(pairs, on=["nt_code", "prestwick_id"], how="inner")
    aucs_summary = (
        aucs.groupby(["nt_code", "prestwick_id"], as_index=False)
        .agg(
            mean_norm_auc_raw=("normauc", "mean"),
            std_norm_auc_raw=("normauc", "std"),
            mean_final_od_raw=("finalod", "mean"),
            n_replicates=("normauc", "size"),
            plate_format=("plate_format", _first_or_nan),
        )
    )

    interactions = pv.merge(aucs_summary, on=["nt_code", "prestwick_id"], how="left")
    interactions["is_primary_panel"] = interactions["nt_code"].isin(primary_codes)

    if human_use_only:
        interactions = interactions[interactions["human_use"]].copy()
    if primary_panel_only:
        interactions = interactions[interactions["is_primary_panel"]].copy()

    interactions = _label_effects(interactions, thresholds)
    interactions["pair_id"] = interactions["prestwick_id"] + "::" + interactions["nt_code"]
    return interactions.reset_index(drop=True)

# build_step1_tables：创建Step 1的药物、微生物、互作和建模表格，并将它们保存到磁盘上，同时生成一个包含统计摘要和输出文件路径的字典作为结果返回
def build_step1_tables(
    raw_dir: str | Path,
    processed_dir: str | Path,
    labels_config_path: str | Path,
    human_use_only: bool = True,
    primary_panel_only: bool = True,
) -> dict:
    """Create the Step 1 drug, microbe, interaction, and modeling tables on disk.

    Args:
        raw_dir: Directory containing the raw Maier 2018 files.
        processed_dir: Output directory for normalized CSV tables and summary JSON.
        labels_config_path: YAML config with Step 1 label thresholds.
        human_use_only: Whether to keep only compounds marked for human use.
        primary_panel_only: Whether to keep only strains in the curated primary panel.

    Returns:
        A summary dictionary describing counts, thresholds, and output file paths.
    """
    raw_dir = Path(raw_dir)
    processed_dir = Path(processed_dir)
    processed_dir.mkdir(parents=True, exist_ok=True)

    thresholds = _read_label_thresholds(Path(labels_config_path))
    drugs = _load_drug_table(raw_dir)
    microbes = _load_microbe_table(raw_dir)
    primary_codes = set(microbes["nt_code"].dropna())
    interactions = _load_interaction_table(
        raw_dir=raw_dir,
        human_use_only=human_use_only,
        primary_panel_only=primary_panel_only,
        primary_codes=primary_codes,
        thresholds=thresholds,
    )
    modeling = (
        interactions.merge(drugs, on="prestwick_id", how="left", suffixes=("", "_drug"))
        .merge(microbes, on="nt_code", how="left", suffixes=("", "_microbe"))
    )

    drugs_path = processed_dir / "step1_drug_table.csv"
    microbes_path = processed_dir / "step1_microbe_table.csv"
    interactions_path = processed_dir / "step1_interactions.csv"
    modeling_path = processed_dir / "step1_modeling_table.csv"

    drugs.to_csv(drugs_path, index=False)
    microbes.to_csv(microbes_path, index=False)
    interactions.to_csv(interactions_path, index=False)
    modeling.to_csv(modeling_path, index=False)

    label_counts = interactions["effect_label"].value_counts(dropna=False).to_dict()
    binary_counts = interactions["binary_effect_label"].value_counts(dropna=False).to_dict()
    summary = {
        "raw_dir": str(raw_dir),
        "processed_dir": str(processed_dir),
        "human_use_only": human_use_only,
        "primary_panel_only": primary_panel_only,
        "thresholds": thresholds,
        "n_drugs": int(modeling["prestwick_id"].nunique()),
        "n_microbes": int(modeling["nt_code"].nunique()),
        "n_interactions": int(len(modeling)),
        "effect_label_counts": {key: int(value) for key, value in label_counts.items()},
        "binary_effect_label_counts": {key: int(value) for key, value in binary_counts.items()},
        "candidate_promote_no_fdr_count": int(interactions["candidate_promote_no_fdr"].sum()),
        "reported_hit_count": int(interactions["hit"].sum()),
        "mean_effect_score": float(interactions["effect_score"].mean()),
        "corr_reported_auc_vs_raw_auc": float(
            interactions[["auc", "mean_norm_auc_raw"]].corr().iloc[0, 1]
        ),
        "output_files": {
            "drug_table": str(drugs_path),
            "microbe_table": str(microbes_path),
            "interactions": str(interactions_path),
            "modeling_table": str(modeling_path),
        },
    }

    summary_path = processed_dir / "step1_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    return summary
