from __future__ import annotations

import re

import numpy as np
import pandas as pd
try:
    from rdkit import Chem
except Exception:  # pragma: no cover - optional fallback for environments without RDKit
    Chem = None


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


def _pick_smiles_text(row: pd.Series) -> str:
    """Pick the first usable SMILES-like field from a row."""
    for column in ["main_component_smiles", "smiles", "canonical_smiles_rdkit"]:
        value = row.get(column)
        if pd.isna(value):
            continue
        text = str(value).strip()
        if text:
            return text
    return ""


def _structure_semantic_hints(smiles: str) -> tuple[str, list[str], list[str]]:
    """Infer lightweight semantic hints from chemical structure when text labels are unavailable."""
    if not smiles or Chem is None:
        return "", [], []
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "", [], []

    aliases: list[str] = []
    keywords: list[str] = []
    family = ""

    # Keep rules intentionally simple/robust: only add high-level functional-group tags
    # that can be consumed by enzyme-prior keyword matching.
    amide = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    ester = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][#6]")
    nitro = Chem.MolFromSmarts("[N+](=O)[O-]")
    azo = Chem.MolFromSmarts("[N]=[N]")
    sulfate_ester = Chem.MolFromSmarts("O[SX4](=O)(=O)O")
    phosphate_ester = Chem.MolFromSmarts("O[P](=O)(O)O")
    methoxy = Chem.MolFromSmarts("[OX2][CH3]")
    aromatic_hydroxyl = Chem.MolFromSmarts("c[OX2H]")

    if amide is not None and mol.HasSubstructMatch(amide):
        aliases.extend(["amide-containing compound"])
        keywords.extend(["amide", "lactam"])
    if ester is not None and mol.HasSubstructMatch(ester):
        aliases.extend(["ester-containing compound"])
        keywords.extend(["ester", "prodrug ester"])
    if nitro is not None and mol.HasSubstructMatch(nitro):
        aliases.extend(["nitro-containing compound"])
        keywords.extend(["nitro", "nitroaromatic"])
    if azo is not None and mol.HasSubstructMatch(azo):
        family = "azo_prodrug_sulfonamide"
        aliases.extend(["azo compound"])
        keywords.extend(["azo", "diazo"])
    if sulfate_ester is not None and mol.HasSubstructMatch(sulfate_ester):
        aliases.extend(["sulfated compound"])
        keywords.extend(["sulfated", "sulfate ester"])
    if phosphate_ester is not None and mol.HasSubstructMatch(phosphate_ester):
        aliases.extend(["phosphorylated compound"])
        keywords.extend(["phosphate", "phosphorylated"])
    if methoxy is not None and mol.HasSubstructMatch(methoxy):
        keywords.append("methoxy")

    aromatic_oh_matches = (
        mol.GetSubstructMatches(aromatic_hydroxyl) if aromatic_hydroxyl is not None else tuple()
    )
    if len(aromatic_oh_matches) >= 2 and not family:
        family = "polyphenol"
        aliases.append("polyphenol-like aromatic")
        keywords.extend(["polyphenol", "phenol"])

    aliases = list(dict.fromkeys(value for value in aliases if value))
    keywords = list(dict.fromkeys(value for value in keywords if value))
    return family, aliases, keywords


def annotate_compound_semantics(frame: pd.DataFrame) -> pd.DataFrame:
    """Attach normalized compound-family hints used by promote and cross-feeding explainers."""
    result = frame.copy()
    normalized_names = []
    families = []
    aliases = []
    keywords = []

    for _, row in result.iterrows():
        chemical_name = "" if pd.isna(row.get("chemical_name")) else str(row.get("chemical_name"))
        therapeutic_class = "" if pd.isna(row.get("therapeutic_class")) else str(row.get("therapeutic_class"))
        therapeutic_effect = "" if pd.isna(row.get("therapeutic_effect")) else str(row.get("therapeutic_effect"))
        raw_text = f"{chemical_name} {therapeutic_class} {therapeutic_effect}".strip().lower()
        compact_text = _normalize_token(raw_text)

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

        structure_family, structure_aliases, structure_keywords = _structure_semantic_hints(_pick_smiles_text(row))
        if not matched_family and structure_family:
            matched_family = structure_family
        if structure_aliases:
            matched_aliases = list(dict.fromkeys(matched_aliases + structure_aliases))
        if structure_keywords:
            matched_keywords = list(dict.fromkeys(matched_keywords + structure_keywords))

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
