from __future__ import annotations

import re

import numpy as np
import pandas as pd


def _normalize_token(value: object) -> str:
    """Normalize free-text chemistry labels into compact alphanumeric tokens."""
    if pd.isna(value):
        return ""
    text = str(value).strip().lower()
    return "".join(character for character in text if character.isalnum())


COMPOUND_FAMILY_RULES = [
    {
        "family": "azo_prodrug_sulfonamide",
        "patterns": ["sulfasalazine", "sasp", "azo prodrug sulfonamide"],
        "aliases": ["sulfasalazine", "sasp", "azo prodrug sulfonamide"],
        "keywords": ["sulfasalazine", "sasp", "azo", "sulfapyridine"],
    },
    {
        "family": "sulfonamide_antifolate",
        "patterns": [
            "sulfapyridine",
            "sulfonamide",
            "sulfamethoxazole",
            "sulfadiazine",
            "sulfisoxazole",
            "cotrimoxazole",
            "co trimoxazole",
            "trimethoprim sulfamethoxazole",
            "antifolate antibacterial",
        ],
        "aliases": [
            "sulfapyridine",
            "sulfonamide",
            "sulfamethoxazole",
            "sulfadiazine",
            "sulfisoxazole",
            "cotrimoxazole",
        ],
        "keywords": ["sulfonamide", "sulfapyridine", "dhps", "paba", "antifolate"],
    },
    {
        "family": "flavonoid_glycoside",
        "patterns": ["qrr", "rutin", "rutinoside", "glycoside flavonoid", "flavonoid glycoside"],
        "aliases": ["qrr", "quercetin-3-o-rutinose-7-o-alpha-l-rhamnoside", "rutinoside flavonoid"],
        "keywords": ["qrr", "rutinoside", "glycoside", "flavonoid"],
    },
    {
        "family": "catechin_gallate",
        "patterns": ["egcg", "epigallocatechin gallate", "catechin gallate"],
        "aliases": ["egcg", "epigallocatechin gallate", "epigallocatechin-3-gallate"],
        "keywords": ["egcg", "catechin", "gallate", "green tea polyphenol"],
    },
    {
        "family": "polyphenol",
        "patterns": ["quercetin", "resveratrol", "naringenin", "hesperidin", "polyphenol", "flavonoid"],
        "aliases": ["quercetin", "resveratrol", "naringenin", "hesperidin", "polyphenol"],
        "keywords": ["polyphenol", "flavonoid", "phenol"],
    },
    {
        "family": "plant_polysaccharide",
        "patterns": ["arabinogalactan", "larch arabinogalactan", "galactan"],
        "aliases": ["arabinogalactan", "larch arabinogalactan"],
        "keywords": ["arabinogalactan", "galactan"],
    },
    {
        "family": "mannan_oligosaccharide",
        "patterns": ["galactomannan", "galactomannans", "manno-oligosaccharide", "beta-manno-oligosaccharide"],
        "aliases": ["galactomannan", "galactomannans", "manno-oligosaccharides", "beta-manno-oligosaccharides"],
        "keywords": ["galactomannan", "mannan", "mannooligosaccharide", "oligosaccharide"],
    },
]


def annotate_compound_semantics(frame: pd.DataFrame) -> pd.DataFrame:
    """Attach normalized compound-family hints used by promote and cross-feeding explainers."""
    result = frame.copy()
    chemical_name = result.get("chemical_name", pd.Series("", index=result.index)).fillna("").astype(str)
    therapeutic_class = result.get("therapeutic_class", pd.Series("", index=result.index)).fillna("").astype(str)
    therapeutic_effect = result.get("therapeutic_effect", pd.Series("", index=result.index)).fillna("").astype(str)
    combined_text = (chemical_name + " " + therapeutic_class + " " + therapeutic_effect).str.lower()
    normalized_text = combined_text.map(_normalize_token)

    normalized_names = []
    families = []
    aliases = []
    keywords = []

    for raw_text, compact_text in zip(combined_text, normalized_text, strict=False):
        matched_family = ""
        matched_aliases: list[str] = []
        matched_keywords: list[str] = []
        canonical_name = ""

        for rule in COMPOUND_FAMILY_RULES:
            normalized_patterns = [_normalize_token(value) for value in rule["patterns"]]
            if any(pattern and pattern in compact_text for pattern in normalized_patterns):
                matched_family = rule["family"]
                matched_aliases = list(dict.fromkeys(rule["aliases"]))
                matched_keywords = list(dict.fromkeys(rule["keywords"]))
                canonical_name = _normalize_token(rule["aliases"][0]) if rule["aliases"] else ""
                break

        if not canonical_name:
            canonical_name = _normalize_token(raw_text.split()[0]) if raw_text.strip() else ""
        normalized_names.append(canonical_name or np.nan)
        families.append(matched_family or np.nan)
        aliases.append("|".join(matched_aliases) if matched_aliases else np.nan)
        keywords.append("|".join(matched_keywords) if matched_keywords else np.nan)

    if "compound_name_normalized" not in result.columns:
        result["compound_name_normalized"] = pd.Series(normalized_names, index=result.index)
    else:
        result["compound_name_normalized"] = result["compound_name_normalized"].fillna(pd.Series(normalized_names, index=result.index))
    result["compound_semantic_family"] = pd.Series(families, index=result.index)
    result["compound_semantic_aliases"] = pd.Series(aliases, index=result.index)
    result["compound_semantic_keywords"] = pd.Series(keywords, index=result.index)
    return result
