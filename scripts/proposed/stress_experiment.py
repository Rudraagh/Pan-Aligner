from __future__ import annotations

from pathlib import Path

from common import OUTPUTS_DIR, ensure_dir, write_json
from visualize import visualize_graph

from .graph_constructor import write_proposed_graph_bundle
from .quality_feedback import run_pg_scunk_feedback_loop, write_pg_scunk_reports
from .workflow import load_assembly_inputs


DEMO_QUALITY_THRESHOLD = 0.99
DEMO_MAX_ROUNDS = 3


def run_pg_scunk_stress_experiment(
    *,
    output_dir: Path | None = None,
    quality_threshold: float = DEMO_QUALITY_THRESHOLD,
    max_rounds: int = DEMO_MAX_ROUNDS,
) -> dict:
    """Run a controlled synthetic stress experiment, not a biological benchmark.

    The fixed four-haplotype fixture contains shared sequence plus controlled
    substitutions, so it creates shared/private minimizer structure and branching.
    No quality is modified: every PG-SCUnK candidate invokes the real graph builder.
    The stringent target intentionally makes the initial graph fail so retuning,
    reconstruction, comparison, and best-so-far behavior are observable.
    """

    root = Path(__file__).resolve().parents[2]
    manifest = root / "tests" / "data" / "proposed" / "toy_manifest.json"
    destination = ensure_dir(output_dir or (OUTPUTS_DIR / "proposed_step5_pgscunk_demo"))
    assemblies = load_assembly_inputs(manifest)
    feedback = run_pg_scunk_feedback_loop(
        assemblies,
        backend="cpu",
        quality_threshold=quality_threshold,
        max_rounds=max_rounds,
    )
    graph = feedback["best_graph"]
    bundle = write_proposed_graph_bundle(graph, destination)
    visualization = visualize_graph(Path(bundle["graph_gfa"]), destination / "controlled_stress_graph")
    reports = write_pg_scunk_reports(feedback, destination)
    metadata = {
        "experiment_type": "controlled_pg_scunk_stress_experiment",
        "label": "DEMONSTRATION ONLY — not a biological benchmark",
        "dataset_manifest": str(manifest.resolve()),
        "assembly_count": len(assemblies),
        "dataset_design": "Four short haplotypes with controlled substitutions; they retain shared sequence while inducing private/shared minimizer nodes and branching.",
        "initial_configuration": feedback["evaluations"][0]["parameters"],
        "target_quality_threshold": quality_threshold,
        "initial_quality": feedback["initial_quality"],
        "initial_below_target": feedback["initial_quality"] < quality_threshold,
        "all_candidates_rebuilt_with_real_graph_constructor": True,
        "quality_scores_recomputed_from_graph": True,
        "winning_configuration_not_hard_coded": True,
        "best_configuration": feedback["best_configuration"],
        "best_quality": feedback["best_overall_quality"],
        "final_quality": feedback["final_quality"],
        "unique_configurations_evaluated": feedback["unique_configurations_evaluated"],
        "stopping_reason": feedback["stagnation_reason"],
    }
    metadata_path = destination / "experiment_metadata.json"
    write_json(metadata_path, metadata)
    return {
        "metadata": metadata,
        "metadata_json": str(metadata_path.resolve()),
        "feedback": feedback,
        **bundle,
        **reports,
        **visualization,
    }
