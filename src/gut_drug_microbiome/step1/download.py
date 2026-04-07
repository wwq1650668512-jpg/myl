from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path


REQUIRED_FIGSHARE_FILES = {
    "aucs.tsv": {
        "url": "https://ndownloader.figshare.com/files/9960184",
        "size": 18876579,
        "md5": "f70438771721f9553d298ea9374aa847",
    },
    "combined_hits.tsv": {
        "url": "https://ndownloader.figshare.com/files/9960172",
        "size": 38004,
        "md5": "6b7d73fbf565f706a84d5de1dd5cd42e",
    },
    "combined_pv.tsv": {
        "url": "https://ndownloader.figshare.com/files/9960175",
        "size": 6783957,
        "md5": "e90a348ff373f30871aac70324139924",
    },
    "compound_properties.tsv": {
        "url": "https://ndownloader.figshare.com/files/7970476",
        "size": 191617,
        "md5": "de027b0f4b676686c0b0f78167c7da07",
    },
    "prestwick_atc.tsv": {
        "url": "https://ndownloader.figshare.com/files/7970461",
        "size": 351140,
        "md5": "487a00365de8213f84de7ff29298f575",
    },
    "species_overview.tsv": {
        "url": "https://ndownloader.figshare.com/files/7970428",
        "size": 9369,
        "md5": "89cda76c00635f8b75611fe95a1e0f11",
    },
}

REQUIRED_SPRINGER_FILES = {
    "Supplementary_table_1.xlsx": "https://static-content.springer.com/esm/art%3A10.1038%2Fnature25979/MediaObjects/41586_2018_BFnature25979_MOESM3_ESM.xlsx",
    "Supplementary_table_2.xlsx": "https://static-content.springer.com/esm/art%3A10.1038%2Fnature25979/MediaObjects/41586_2018_BFnature25979_MOESM4_ESM.xlsx",
    "Supplementary_table_3.xlsx": "https://static-content.springer.com/esm/art%3A10.1038%2Fnature25979/MediaObjects/41586_2018_BFnature25979_MOESM5_ESM.xlsx",
    "Source_data_Fig1.xlsx": "https://static-content.springer.com/esm/art%3A10.1038%2Fnature25979/MediaObjects/41586_2018_BFnature25979_MOESM13_ESM.xlsx",
}


def _md5(path: Path) -> str:
    """Compute the MD5 checksum of a local file for integrity checks."""
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download_file(url: str, destination: Path) -> None:
    """Download a file with curl into the requested destination path."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["curl", "-L", "--max-time", "120", "-o", str(destination), url],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def download_step1_data(raw_dir: str | Path, overwrite: bool = False) -> dict:
    """Download the required Step 1 raw files and write a manifest summary.

    Args:
        raw_dir: Destination directory for the Step 1 source files.
        overwrite: Whether to redownload files even if they already exist.

    Returns:
        A manifest dictionary describing downloaded files and local paths.
    """
    raw_dir = Path(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)

    manifest: dict[str, dict] = {"figshare": {}, "springer": {}}

    for name, metadata in sorted(REQUIRED_FIGSHARE_FILES.items()):
        destination = raw_dir / name
        needs_download = overwrite or not destination.exists()
        if destination.exists() and not overwrite:
            needs_download = _md5(destination) != metadata["md5"]
        if needs_download:
            _download_file(metadata["url"], destination)
        manifest["figshare"][name] = {
            "url": metadata["url"],
            "size": metadata["size"],
            "md5": metadata["md5"],
            "local_path": str(destination),
        }

    for name, url in sorted(REQUIRED_SPRINGER_FILES.items()):
        destination = raw_dir / name
        if overwrite or not destination.exists():
            _download_file(url, destination)
        manifest["springer"][name] = {
            "url": url,
            "local_path": str(destination),
            "size": destination.stat().st_size if destination.exists() else None,
        }

    manifest_path = raw_dir / "download_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return manifest
