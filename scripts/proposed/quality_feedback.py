from __future__ import annotations

import csv
import json
from pathlib import Path

from .anchor_selection import DEFAULT_ALPHA, DEFAULT_BETA, DEFAULT_GAMMA
from .graph_constructor import build_proposed_pangenome_graph
from .models import AssemblyInput, ProposedGraph

CONFIGURATION_FIELDS = ("k", "window_size", "min_node_support", "max_anchors", "min_anchor_spacing", "anchor_alpha", "anchor_beta", "anchor_gamma")
DEFAULT_IMPROVEMENT_TOLERANCE = 1e-4
DEFAULT_STAGNATION_PATIENCE = 2
DIAGNOSIS_MAPPING = {
    "anchor_coverage": {"action": "retune_anchor_selection", "parameters": ["min_anchor_spacing", "max_anchors", "anchor_alpha", "anchor_beta", "anchor_gamma"], "reason": "Selected-anchor coverage is low; vary spacing, anchor count, and ranking weights."},
    "node_consistency": {"action": "retune_minimizer_nodes", "parameters": ["k", "window_size", "min_node_support"], "reason": "Node context coherence is low; vary k, window size, and minimum node support."},
    "node_support": {"action": "retune_node_support", "parameters": ["k", "window_size", "min_node_support"], "reason": "Cross-path node support is low; vary k, window size, and minimum node support."},
    "edge_structural_coherence": {"action": "retune_graph_transitions", "parameters": ["k", "window_size", "min_node_support"], "reason": "Observed-path edge support is weak; vary k, window size, and minimum node support."},
    "path_structural_coherence": {"action": "retune_path_resolution", "parameters": ["k", "window_size"], "reason": "Observed-path structural coherence is weak; vary k and window size."},
    "branch_coherence": {"action": "retune_graph_branching", "parameters": ["k", "window_size", "min_node_support"], "reason": "Branching is high; vary k, window size, and minimum node support."},
}


def _quality_passes(quality: dict, threshold: float, min_anchor_coverage: float) -> bool:
    return quality["overall_quality_score"] >= threshold and quality["anchor_coverage_ratio"] >= min_anchor_coverage and quality["mean_node_consistency"] >= 0.45 and quality["mean_path_consistency"] >= 0.45


def _adaptive_default_parameters(assemblies: list[AssemblyInput]) -> dict:
    shortest = min(len(item.sequence) for item in assemblies)
    k = max(3, min(11, shortest // 2))
    return {"k": k, "window_size": max(2, min(6, k - 1 if k > 3 else k)), "min_node_support": 2 if len(assemblies) >= 3 else 1, "max_anchors": min(12, max(4, len(assemblies) * 2)), "min_anchor_spacing": 4.0 if shortest >= 40 else 2.0, "anchor_alpha": DEFAULT_ALPHA, "anchor_beta": DEFAULT_BETA, "anchor_gamma": DEFAULT_GAMMA}


def canonical_configuration(parameters: dict) -> tuple[tuple[str, int | float], ...]:
    """Immutable key containing every active graph/anchor construction setting."""
    return tuple((field, round(float(parameters[field]), 12) if isinstance(parameters[field], float) else int(parameters[field])) for field in CONFIGURATION_FIELDS)


def configuration_key_text(key: tuple[tuple[str, int | float], ...]) -> str:
    return json.dumps(dict(key), sort_keys=True, separators=(",", ":"))


def diagnose_quality(quality: dict) -> dict:
    """Deterministically diagnose a weak signal without altering Step 2 formulas.

    Edge/path scores are observed-path structural coherence, not independent
    biological correctness measures, because graph edges originate from paths.
    """
    signals = {"anchor_coverage": quality.get("anchor_coverage_ratio", 0.0), "node_consistency": quality.get("mean_node_consistency", 0.0), "node_support": quality.get("mean_node_support_ratio", 0.0), "edge_structural_coherence": quality.get("mean_edge_consistency", 0.0), "path_structural_coherence": quality.get("mean_path_consistency", 0.0), "branch_coherence": 1.0 - min(1.0, quality.get("branching_node_ratio", 1.0))}
    priority = ("anchor_coverage", "node_consistency", "node_support", "edge_structural_coherence", "path_structural_coherence", "branch_coherence")
    weakest = min(priority, key=lambda name: (signals[name], priority.index(name)))
    mapping = DIAGNOSIS_MAPPING[weakest]
    return {"weakest_component": weakest, "signals": signals, "action": mapping["action"], "reason": mapping["reason"], "parameters": mapping["parameters"]}


def _boost_alpha(parameters: dict) -> tuple[float, float, float]:
    alpha, beta, gamma = parameters["anchor_alpha"] + 0.10, parameters["anchor_beta"], parameters["anchor_gamma"]
    total = alpha + beta + gamma
    return tuple(round(value / total, 12) for value in (alpha, beta, gamma))


def generate_feedback_candidates(parameters: dict, diagnosis: dict, visited: set[tuple[tuple[str, int | float], ...]]) -> list[dict]:
    """Generate a bounded, deterministic neighbourhood and remove seen configs."""
    action = diagnosis["action"]
    if action == "retune_anchor_selection":
        alpha, beta, gamma = _boost_alpha(parameters)
        raw = [{**parameters, "min_anchor_spacing": max(1.0, float(parameters["min_anchor_spacing"]) / 2.0)}, {**parameters, "max_anchors": int(parameters["max_anchors"]) + 2}, {**parameters, "anchor_alpha": alpha, "anchor_beta": beta, "anchor_gamma": gamma}]
    elif action == "retune_path_resolution":
        raw = [{**parameters, "k": max(3, int(parameters["k"]) - 2)}, {**parameters, "window_size": max(2, int(parameters["window_size"]) - 1)}, {**parameters, "k": max(3, int(parameters["k"]) - 2), "window_size": max(2, int(parameters["window_size"]) - 1)}]
    else:
        raw = [{**parameters, "k": max(3, int(parameters["k"]) - 2)}, {**parameters, "window_size": max(2, int(parameters["window_size"]) - 1)}, {**parameters, "min_node_support": max(1, int(parameters["min_node_support"]) - 1)}]
    candidates, local_seen = [], set()
    for candidate in raw:
        key = canonical_configuration(candidate)
        if key not in visited and key not in local_seen:
            local_seen.add(key)
            candidates.append({"parameters": candidate, "configuration_key": key, "action": action, "reason": diagnosis["reason"]})
    return candidates


def _evaluation(round_index: int, candidate_id: str, candidate: dict, graph: ProposedGraph, quality_delta: float | None, improved: bool) -> dict:
    quality, diagnosis = graph.quality, diagnose_quality(graph.quality)
    return {"round": round_index, "candidate_id": candidate_id, "configuration_key": configuration_key_text(candidate["configuration_key"]), "parameters": dict(candidate["parameters"]), "overall_quality": quality["overall_quality_score"], "mean_node_consistency": quality["mean_node_consistency"], "mean_edge_consistency": quality["mean_edge_consistency"], "mean_path_consistency": quality["mean_path_consistency"], "anchor_coverage": quality["anchor_coverage_ratio"], "node_count": quality["node_count"], "edge_count": quality["edge_count"], "path_count": len(graph.paths), "anchor_count": quality["anchor_count"], "accepted": False, "is_best_so_far": improved, "quality_delta": quality_delta, "weakest_component": diagnosis["weakest_component"], "action": candidate["action"], "reason": candidate["reason"]}


def run_pg_scunk_feedback_loop(assemblies: list[AssemblyInput], *, base_parameters: dict | None = None, backend: str = "auto", quality_threshold: float = 0.72, min_anchor_coverage: float = 0.70, max_rounds: int = 3, improvement_tolerance: float = DEFAULT_IMPROVEMENT_TOLERANCE, stagnation_patience: int = DEFAULT_STAGNATION_PATIENCE) -> dict:
    """Run deterministic feedback search; always return the protected best-so-far graph.

    Meaningful improvement is Q > best_Q + tolerance. Stops on threshold success,
    max rounds, no unseen configurations, or stagnation patience.
    """
    if not assemblies:
        raise ValueError("At least one assembly is required.")
    if max_rounds < 1 or stagnation_patience < 1 or improvement_tolerance < 0:
        raise ValueError("max_rounds and stagnation_patience must be positive; tolerance must be non-negative.")
    initial = _adaptive_default_parameters(assemblies)
    if base_parameters:
        initial.update(base_parameters)
    initial_key = canonical_configuration(initial)
    pending = [{"parameters": initial, "configuration_key": initial_key, "action": "initial_evaluation", "reason": "Evaluate adaptive/default configuration."}]
    visited, evaluations = set(), []
    best_graph = best_configuration = best_evaluation = final_graph = final_evaluation = None
    best_quality, threshold_reached, stagnation_rounds, reason, round_index, sequence = float("-inf"), False, 0, "max_rounds_reached", 0, 0
    while pending and round_index < max_rounds:
        round_index += 1
        round_improved, round_records = False, []
        for candidate in pending:
            if candidate["configuration_key"] in visited:
                continue
            visited.add(candidate["configuration_key"]); sequence += 1
            graph = build_proposed_pangenome_graph(assemblies, backend=backend, **candidate["parameters"])
            delta = None if best_graph is None else graph.quality["overall_quality_score"] - best_quality
            improved = best_graph is None or graph.quality["overall_quality_score"] > best_quality + improvement_tolerance
            record = _evaluation(round_index, f"c{sequence:03d}", candidate, graph, delta, improved)
            record["accepted"] = _quality_passes(graph.quality, quality_threshold, min_anchor_coverage)
            evaluations.append(record); round_records.append(record); final_graph, final_evaluation = graph, record
            if improved:
                best_graph, best_configuration, best_quality, best_evaluation, round_improved = graph, dict(candidate["parameters"]), graph.quality["overall_quality_score"], record, True
            if record["accepted"]:
                threshold_reached, reason = True, "quality_threshold_reached"
                break
        if threshold_reached:
            break
        if not round_records:
            reason = "no_unseen_configurations"; break
        stagnation_rounds = 0 if round_improved else stagnation_rounds + 1
        if stagnation_rounds >= stagnation_patience:
            reason = "stagnation_patience_exceeded"; break
        diagnosis = diagnose_quality(best_graph.quality)
        pending = generate_feedback_candidates(best_configuration, diagnosis, visited)
        if not pending:
            reason = "no_unseen_configurations"; break
    if best_graph is None or best_evaluation is None or final_evaluation is None:
        raise RuntimeError("PG-SCUnK feedback loop could not build any proposed graph.")
    return {"accepted": threshold_reached, "threshold_reached": threshold_reached, "initial_quality": evaluations[0]["overall_quality"], "best_quality": best_graph.quality, "best_overall_quality": best_quality, "final_quality": final_evaluation["overall_quality"], "quality_threshold": quality_threshold, "min_anchor_coverage": min_anchor_coverage, "improvement_tolerance": improvement_tolerance, "stagnation_patience": stagnation_patience, "stagnation_reason": reason, "best_round": best_evaluation["round"], "total_rounds": round_index, "unique_configurations_evaluated": len(visited), "best_configuration": best_configuration, "final_configuration": final_evaluation["parameters"], "final_evaluation": final_evaluation, "best_evaluation": best_evaluation, "evaluations": evaluations, "final_parameters": best_configuration, "round_summaries": [{"round_index": number, "evaluations": [item for item in evaluations if item["round"] == number]} for number in range(1, round_index + 1)], "graph": best_graph, "best_graph": best_graph, "final_graph": final_graph}


def write_pg_scunk_reports(result: dict, output_dir: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    iterations_path, summary_path, plot_path = output_dir / "pg_scunk_iterations.csv", output_dir / "pg_scunk_summary.json", output_dir / "pg_scunk_quality.png"
    fields = ["round", "candidate_id", "configuration_key", "parameters", "overall_quality", "mean_node_consistency", "mean_edge_consistency", "mean_path_consistency", "anchor_coverage", "node_count", "edge_count", "path_count", "anchor_count", "accepted", "is_best_so_far", "quality_delta", "weakest_component", "action", "reason"]
    with iterations_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader()
        for evaluation in result["evaluations"]:
            row = dict(evaluation); row["parameters"] = json.dumps(row["parameters"], sort_keys=True); writer.writerow(row)
    summary = {"initial_quality": result["initial_quality"], "best_quality": result["best_overall_quality"], "final_quality": result["final_quality"], "threshold": result["quality_threshold"], "threshold_reached": result["threshold_reached"], "best_round": result["best_round"], "total_rounds": result["total_rounds"], "unique_configurations_evaluated": result["unique_configurations_evaluated"], "stagnation_reason": result["stagnation_reason"], "best_configuration": result["best_configuration"], "final_configuration": result["final_configuration"], "best_evaluation": result["best_evaluation"], "final_evaluation": result["final_evaluation"], "diagnostic_mapping": DIAGNOSIS_MAPPING, "quality_interpretation": "Path and edge consistency are observed-path structural coherence measures; because edges are constructed from observed paths, they are not independent biological correctness evidence."}
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    try:
        import matplotlib.pyplot as plt
        values = [item["overall_quality"] for item in result["evaluations"]]
        figure, axis = plt.subplots(figsize=(7, 4)); axis.plot(range(1, len(values) + 1), values, marker="o", label="candidate quality"); axis.axhline(result["quality_threshold"], color="tab:red", linestyle="--", label="threshold"); axis.set(xlabel="Evaluation", ylabel="Overall quality", title="PG-SCUnK quality trajectory"); axis.set_ylim(0, 1.02); axis.legend(); figure.tight_layout(); figure.savefig(plot_path, dpi=160); plt.close(figure)
    except Exception:
        plot_path = None
    reports = {"pg_scunk_iterations_csv": str(iterations_path.resolve()), "pg_scunk_summary_json": str(summary_path.resolve())}
    if plot_path is not None:
        reports["pg_scunk_quality_plot"] = str(plot_path.resolve())
    return reports
