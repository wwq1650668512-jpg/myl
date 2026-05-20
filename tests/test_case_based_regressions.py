from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
FUSION_DIR = ROOT / "predictions" / "evaluation" / "fusion_comparison"
PROFILE_PATH = FUSION_DIR / "profile_evidence.json"


def _require_artifacts() -> None:
    required = [
        PROFILE_PATH,
        FUSION_DIR / "rifaximin.csv",
        FUSION_DIR / "vancomycin.csv",
        FUSION_DIR / "lubiprostone.csv",
        FUSION_DIR / "metronidazole.csv",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        pytest.skip(f"case-based regression artifacts missing: {missing}")


def _profile_payload() -> dict[str, dict[str, object]]:
    _require_artifacts()
    return json.loads(PROFILE_PATH.read_text(encoding="utf-8"))


def _weighted_table(stem: str) -> pd.DataFrame:
    _require_artifacts()
    table = pd.read_csv(FUSION_DIR / f"{stem}.csv")
    sub = table[table["fusion_mode"] == "weighted_0.3_0.7"].copy()
    if sub.empty:
        pytest.skip(f"weighted_0.3_0.7 rows missing for {stem}.csv")
    return sub.sort_values("rank").reset_index(drop=True)


def _inhibit_fraction(profile: dict[str, object]) -> float:
    panel_size = float(profile.get("panel_size", 0) or 0)
    inhibit = float(dict(profile.get("step1_counts", {})).get("inhibit", 0))
    if panel_size <= 0:
        return 0.0
    return inhibit / panel_size


def _strong_butyrate_inhibit_count(profile: dict[str, object]) -> tuple[int, int]:
    rows = profile.get("butyrate_rows", []) or []
    strong = 0
    for row in rows:
        label = str(row.get("predicted_effect_label", "")).strip().lower()
        prob = pd.to_numeric(pd.Series([row.get("predicted_inhibit_probability")]), errors="coerce").fillna(0.0).iloc[0]
        score = pd.to_numeric(pd.Series([row.get("predicted_effect_score")]), errors="coerce").fillna(0.0).iloc[0]
        if label == "inhibit" and (prob >= 0.5 or score < 0.0):
            strong += 1
    return strong, len(rows)


def _best_rank(table: pd.DataFrame, patterns: list[str]) -> int | None:
    matched = table[table["disease_name"].astype(str).str.contains("|".join(patterns), case=False, regex=True, na=False)]
    if matched.empty:
        return None
    return int(matched["rank"].min())


def test_rifaximin_not_indiscriminate_on_core_butyrate_producers() -> None:
    profile = _profile_payload()["Rifaximin"]
    strong, total = _strong_butyrate_inhibit_count(profile)
    assert total > 0
    assert strong / total <= 0.4


def test_lubiprostone_not_broad_spectrum_antimicrobial_like() -> None:
    profile = _profile_payload()["Lubiprostone"]
    assert _inhibit_fraction(profile) <= 0.2
    assert int(profile.get("inhibit_prob_ge_0_7", 0) or 0) <= 10


def test_vancomycin_keeps_high_ecological_risk_signature() -> None:
    profile = _profile_payload()["Vancomycin"]
    eco_risk = float(pd.to_numeric(pd.Series([profile.get("ecological_risk_score")]), errors="coerce").fillna(0.0).iloc[0])
    assert eco_risk >= 0.12
    assert _inhibit_fraction(profile) >= 0.25


def test_metronidazole_shows_benefit_risk_tradeoff() -> None:
    profile = _profile_payload()["Metronidazole"]
    table = _weighted_table("metronidazole")
    infection_rank = _best_rank(
        table,
        [r"abscess", r"脓肿", r"\bibd\b", r"炎症性肠病", r"\bcd\b", r"克罗恩", r"fistula", r"感染"],
    )
    eco_risk = float(pd.to_numeric(pd.Series([profile.get("ecological_risk_score")]), errors="coerce").fillna(0.0).iloc[0])
    assert infection_rank is not None and infection_rank <= 5
    assert (eco_risk >= 0.15) or (_inhibit_fraction(profile) >= 0.3)
