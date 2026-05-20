from __future__ import annotations

import re

import pandas as pd


def normalize_whitespace(value: object) -> str:
    """Collapse repeated whitespace and coerce null-like values to an empty string."""
    if value is None or pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value).replace("\xa0", " ").strip())


def canonicalize_key(value: object, keep_cjk: bool = False) -> str:
    """Build a lowercase matching key from free text."""
    text = normalize_whitespace(value).lower()
    pattern = r"[^a-z0-9\u4e00-\u9fff]+" if keep_cjk else r"[^a-z0-9]+"
    return re.sub(pattern, "", text)
