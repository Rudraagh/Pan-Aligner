from __future__ import annotations

import csv
import json
import time
from pathlib import Path

from common import OUTPUTS_DIR, ensure_dir, read_json, write_json
from visualize import visualize_graph

from .graph_constructor import validate_proposed_graph, write_proposed_graph_bundle
from .quality_feedback import run_pg_scunk_feedback_loop, write_pg_scunk_reports
from .workflow import load_assembly_inputs


SUMMARY_FIELDS = [
    "gene", "input_assembly_count", "haplotype_path_count", "node_count", "edge_count", "path_count",
    "rescued_node_count", "private_node_count", "shared_node_count", "anchor_candidate_count", "selected_anchor_count",
    "mean_node_consistency", "mean_edge_consistency", "mean_path_consistency", "graph_consistency",
    "anchor_coverage", "overall_quality", "initial_configuration", "initial_quality", "iterations",
    "unique_configurations", "best_configuration", "best_quality", "final_evaluated_quality", "threshold",
    "threshold_reached", "stopping_reason", "backend_used", "validation_passed",
]


def _in_unit_interval(value: float) -> bool:
    return 0.0 <= float(value) <= 1.0


def validate_final_evaluation(graph, feedback: dict, reports: dict) -> dict:
    """Validate serialised proposed-only metrics without changing their algorithms."""
    errors = list(validate_proposed_graph(graph))
    for collection, field in ((graph.nodes, "consistency_score"), (graph.edges, "consistency_score"), (graph.paths, "consistency_score"), (graph.anchors, "final_score")):
        if any(not _in_unit_interval(getattr(item, field)) for item in collection):
            errors.append(f"Out-of-range {field} in {collection.__class__.__name__}.")
    quality_fields = ("mean_node_consistency", "mean_edge_consistency", "mean_path_consistency", "graph_path_consistency", "anchor_coverage_ratio", "overall_quality_score")
    if any(not _in_unit_interval(graph.quality[field]) for field in quality_fields):
        errors.append("Out-of-range graph quality metric.")
    keys = [item["configuration_key"] for item in feedback["evaluations"]]
    if len(keys) != len(set(keys)):
        errors.append("Duplicate PG-SCUnK configuration recorded.")
    spacing = graph.parameters["min_anchor_spacing"]
    selected = sorted(graph.anchors, key=lambda item: item.mean_position)
    if any(right.mean_position - left.mean_position < spacing for left, right in zip(selected, selected[1:])):
        errors.append("Selected anchors violate configured spacing.")
    best = float("-inf")
    for item in feedback["evaluations"]:
        if item["is_best_so_far"]:
            if item["overall_quality"] <= best + feedback["improvement_tolerance"] and best != float("-inf"):
                errors.append("Best-so-far was replaced without meaningful improvement.")
            best = item["overall_quality"]
    if abs(best - feedback["best_overall_quality"]) > 1e-12:
        errors.append("Best-so-far quality disagrees with evaluation trajectory.")
    summary = read_json(Path(reports["pg_scunk_summary_json"]))
    rows = list(csv.DictReader(Path(reports["pg_scunk_iterations_csv"]).open(encoding="utf-8")))
    if len(rows) != len(feedback["evaluations"]) or summary["best_quality"] != feedback["best_overall_quality"]:
        errors.append("PG-SCUnK JSON and CSV reports disagree.")
    return {"passed": not errors, "errors": errors}


def _row(gene: str, graph, feedback: dict, validation: dict) -> dict:
    metadata, quality = graph.metadata, graph.quality
    return {
        "gene": gene, "input_assembly_count": graph.assembly_count, "haplotype_path_count": len(graph.paths),
        "node_count": quality["node_count"], "edge_count": quality["edge_count"], "path_count": len(graph.paths),
        "rescued_node_count": metadata["rescued_node_count"], "private_node_count": metadata["haplotype_specific_node_count"],
        "shared_node_count": metadata["shared_node_count"], "anchor_candidate_count": len(graph.anchor_candidates),
        "selected_anchor_count": len(graph.anchors), "mean_node_consistency": quality["mean_node_consistency"],
        "mean_edge_consistency": quality["mean_edge_consistency"], "mean_path_consistency": quality["mean_path_consistency"],
        "graph_consistency": quality["graph_path_consistency"], "anchor_coverage": quality["anchor_coverage_ratio"],
        "overall_quality": quality["overall_quality_score"], "initial_configuration": json.dumps(feedback["evaluations"][0]["parameters"], sort_keys=True),
        "initial_quality": feedback["initial_quality"], "iterations": feedback["total_rounds"],
        "unique_configurations": feedback["unique_configurations_evaluated"], "best_configuration": json.dumps(feedback["best_configuration"], sort_keys=True),
        "best_quality": feedback["best_overall_quality"], "final_evaluated_quality": feedback["final_quality"],
        "threshold": feedback["quality_threshold"], "threshold_reached": feedback["threshold_reached"],
        "stopping_reason": feedback["stagnation_reason"], "backend_used": graph.backend_used,
        "validation_passed": validation["passed"],
    }


def _write_readme(destination: Path, rows: list[dict], runtime_seconds: float, baseline_exists: bool) -> None:
    table = "\n".join(
        f"| {row['gene']} | {row['input_assembly_count']} | {row['node_count']} | {row['edge_count']} | {row['overall_quality']:.6f} | {row['backend_used']} |"
        for row in rows
    )
    baseline_note = (
        "Existing baseline files were read only and are described separately in `baseline_proposed_comparison.md`."
        if baseline_exists else "No baseline evaluation bundle was found at report time."
    )
    (destination / "README.md").write_text(
        "# Proposed-methodology final evaluation\n\n"
        "This is a proposed-only, alignment-free evaluation bundle. It is not a biological benchmark and does not overwrite or merge with the minigraph/PanAligner baseline.\n\n"
        "## Datasets and results\n\n| Gene | Assemblies / paths | Nodes | Edges | Overall structural quality | Backend |\n|---|---:|---:|---:|---:|---|\n"
        f"{table}\n\n"
        "All assemblies currently listed in `data/metadata/train_manifest.json` were used, without the older eight-per-bucket cap.\n\n"
        "## Pipeline\n\nFASTA inputs → deterministic CPU minimizers → path-preserving proposed graph → Step 2 structural/path coherence → Step 3 anchors → Step 4 PG-SCUnK controller.\n\n"
        "## Interpretation and limitations\n\nGraph/path consistency is structural coherence derived from observed minimizer paths, not independent biological correctness. The proposed method has no GPU acceleration in this evaluation; the recorded backend is CPU. Baseline alignment metrics and proposed structural metrics are not directly interchangeable.\n\n"
        f"Runtime: {runtime_seconds:.2f} seconds. {baseline_note}\n",
        encoding="utf-8",
    )


def run_final_proposed_evaluation(*, output_dir: Path | None = None) -> dict:
    """Run all available APP/PSEN1/PSEN2 inputs with the unchanged proposed stack."""
    root = Path(__file__).resolve().parents[2]
    manifest = root / "data" / "metadata" / "train_manifest.json"
    destination = ensure_dir(output_dir or (OUTPUTS_DIR / "proposed_final_evaluation"))
    started = time.perf_counter()
    assemblies = load_assembly_inputs(manifest, max_assemblies_per_bucket=None)
    grouped: dict[str, list] = {}
    for assembly in assemblies:
        grouped.setdefault(assembly.gene, []).append(assembly)
    rows, details = [], {}
    for gene in ("APP", "PSEN1", "PSEN2"):
        gene_assemblies = grouped[gene]
        gene_dir = ensure_dir(destination / gene)
        feedback = run_pg_scunk_feedback_loop(gene_assemblies, backend="cpu")
        graph = feedback["best_graph"]
        bundle = write_proposed_graph_bundle(graph, gene_dir)
        visualization = visualize_graph(Path(bundle["graph_gfa"]), gene_dir / f"{gene.lower()}.proposed_graph")
        reports = write_pg_scunk_reports(feedback, gene_dir)
        validation = validate_final_evaluation(graph, feedback, reports)
        write_json(gene_dir / "validation.json", validation)
        row = _row(gene, graph, feedback, validation)
        rows.append(row)
        details[gene] = {"summary": row, "validation": validation, "artifacts": {**bundle, **visualization, **reports}}
    runtime = time.perf_counter() - started
    summary_csv = destination / "proposed_final_summary.csv"
    with summary_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS); writer.writeheader(); writer.writerows(rows)
    summary = {"evaluation_type": "proposed_only_structural_evaluation", "dataset_manifest": str(manifest.resolve()), "complete_inputs_used": True, "runtime_seconds": runtime, "rows": rows, "genes": details}
    summary_json = destination / "proposed_final_summary.json"
    write_json(summary_json, summary)
    baseline_metrics = root / "outputs" / "evaluation" / "alignment_metrics.json"
    comparison = destination / "baseline_proposed_comparison.md"
    baseline_text = "Baseline metrics were not available."
    if baseline_metrics.exists():
        baseline = read_json(baseline_metrics)
        baseline_text = "Baseline results exist at `outputs/evaluation/` and were not modified.\n\n```json\n" + json.dumps(baseline, indent=2) + "\n```"
    comparison.write_text(
        "# Pipeline separation\n\n## BASELINE\n\nminigraph/PanAligner workflow; alignment-oriented results retained unchanged.\n\n"
        "## PROPOSED\n\nalignment-free minimizers + haplotype/path structural coherence + intelligent anchors + PG-SCUnK. These structural metrics are not direct alignment-accuracy comparisons, so this report makes no superiority claim.\n\n"
        + baseline_text + "\n",
        encoding="utf-8",
    )
    _write_readme(destination, rows, runtime, baseline_metrics.exists())
    return {"output_dir": str(destination.resolve()), "summary_csv": str(summary_csv.resolve()), "summary_json": str(summary_json.resolve()), "runtime_seconds": runtime, "rows": rows, "details": details}
