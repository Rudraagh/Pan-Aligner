from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from common import METADATA_DIR, OUTPUTS_DIR
from proposed.workflow import run_proposed_methodology
from proposed.stress_experiment import run_pg_scunk_stress_experiment
from proposed.final_evaluation import run_final_proposed_evaluation


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run the proposed alignment-free methodology stack with minimizer computation, "
            "haplotype-consistency-aware graph construction, intelligent anchors, and PG-SCUnK feedback."
        )
    )
    parser.add_argument(
        "--dataset-manifest",
        type=Path,
        default=METADATA_DIR / "train_manifest.json",
        help="Dataset manifest providing assemblies. Defaults to the existing train manifest.",
    )
    parser.add_argument("--gene", type=str, default=None, help="Optional gene filter.")
    parser.add_argument("--pg-scunk-stress-demo", action="store_true", help="Run the controlled synthetic PG-SCUnK stress demonstration, not a biological benchmark.")
    parser.add_argument("--final-proposed-evaluation", action="store_true", help="Run the complete proposed-only APP/PSEN1/PSEN2 evaluation bundle.")
    parser.add_argument("--demo-quality-threshold", type=float, default=0.99, help="Target quality used only by the controlled PG-SCUnK stress demonstration.")
    parser.add_argument("--backend", type=str, default="auto", help="Minimizer backend: auto, cpu, or gpu.")
    parser.add_argument(
        "--max-assemblies-per-bucket",
        type=int,
        default=8,
        help="Optional cap per healthy/unhealthy bucket to keep the proposed workflow lightweight.",
    )
    parser.add_argument(
        "--quality-threshold",
        type=float,
        default=0.72,
        help="Target quality score for the PG-SCUnK feedback controller.",
    )
    parser.add_argument("--max-rounds", type=int, default=3, help="Maximum PG-SCUnK retuning rounds.")
    parser.add_argument("--improvement-tolerance", type=float, default=1e-4, help="Minimum strict Q improvement needed to replace best-so-far.")
    parser.add_argument("--stagnation-patience", type=int, default=2, help="Consecutive non-improving rounds allowed before stopping.")
    parser.add_argument("--anchor-alpha", type=float, default=0.35, help="Anchor uniqueness weight; weights must sum to 1.")
    parser.add_argument("--anchor-beta", type=float, default=0.35, help="Anchor depth weight; weights must sum to 1.")
    parser.add_argument("--anchor-gamma", type=float, default=0.30, help="Anchor consistency weight; weights must sum to 1.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for proposed-methodology artifacts; the stress demo defaults to outputs/proposed_step5_pgscunk_demo.",
    )
    args = parser.parse_args()

    if args.pg_scunk_stress_demo:
        result = run_pg_scunk_stress_experiment(
            output_dir=args.output_dir.resolve() if args.output_dir else None,
            quality_threshold=args.demo_quality_threshold,
            max_rounds=args.max_rounds,
        )
        print(json.dumps({key: value for key, value in result.items() if key != "feedback"}, indent=2))
        return
    if args.final_proposed_evaluation:
        result = run_final_proposed_evaluation(output_dir=args.output_dir.resolve() if args.output_dir else None)
        print(json.dumps({key: value for key, value in result.items() if key != "details"}, indent=2))
        return

    result = run_proposed_methodology(
        args.dataset_manifest.resolve(),
        output_dir=args.output_dir.resolve() if args.output_dir else None,
        gene=args.gene,
        backend=args.backend,
        max_assemblies_per_bucket=args.max_assemblies_per_bucket,
        quality_threshold=args.quality_threshold,
        max_rounds=args.max_rounds,
        improvement_tolerance=args.improvement_tolerance,
        stagnation_patience=args.stagnation_patience,
        anchor_alpha=args.anchor_alpha,
        anchor_beta=args.anchor_beta,
        anchor_gamma=args.anchor_gamma,
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
