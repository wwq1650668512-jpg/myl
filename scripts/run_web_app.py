from __future__ import annotations

import argparse
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.web import serve_web_app


DEFAULT_EXPANDED_INTEGRATED_PREDICTIONS = ROOT / "predictions/step2/baseline_scaffold_v1_83/predictions.csv"
DEFAULT_LEGACY_INTEGRATED_PREDICTIONS = ROOT / "predictions/step2/baseline_scaffold_v1/predictions.csv"
DEFAULT_INTEGRATED_PREDICTIONS = (
    DEFAULT_EXPANDED_INTEGRATED_PREDICTIONS
    if DEFAULT_EXPANDED_INTEGRATED_PREDICTIONS.exists()
    else DEFAULT_LEGACY_INTEGRATED_PREDICTIONS
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the local Step 1/2/3 visualization web app.")
    parser.add_argument(
        "--host",
        default="127.0.0.1",
        help="Host to bind the local web server.",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8080,
        help="Port used by the local web server.",
    )
    parser.add_argument(
        "--integrated-predictions",
        default=DEFAULT_INTEGRATED_PREDICTIONS,
        type=Path,
        help="Integrated Step 1 + Step 2 prediction table consumed by the UI.",
    )
    parser.add_argument(
        "--demo-ranking",
        default=ROOT / "predictions/step3/candidate_screen_demo/candidate_ranking.csv",
        type=Path,
        help="Optional Step 3 demo ranking table shown on the landing page.",
    )
    parser.add_argument(
        "--static-dir",
        default=ROOT / "webapp" / "static",
        type=Path,
        help="Static asset directory for the local webpage.",
    )
    args = parser.parse_args()
    demo_ranking = args.demo_ranking if args.demo_ranking.exists() else None
    serve_web_app(
        host=args.host,
        port=args.port,
        integrated_predictions_path=args.integrated_predictions,
        demo_ranking_path=demo_ranking,
        static_dir=args.static_dir,
    )


if __name__ == "__main__":
    main()
