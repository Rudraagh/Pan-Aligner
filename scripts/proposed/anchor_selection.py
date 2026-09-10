from __future__ import annotations

import csv
import json
import math
from collections import Counter
from pathlib import Path

from .models import AnchorCandidate, ProposedGraph


DEFAULT_ALPHA = 0.35
DEFAULT_BETA = 0.35
DEFAULT_GAMMA = 0.30


def validate_anchor_weights(alpha: float, beta: float, gamma: float) -> None:
    if min(alpha, beta, gamma) < 0:
        raise ValueError("Anchor weights must be non-negative.")
    if not math.isclose(alpha + beta + gamma, 1.0, rel_tol=0.0, abs_tol=1e-9):
        raise ValueError("Anchor weights must sum to 1.0.")


def compute_anchor_score(uniqueness: float, depth: float, consistency: float, alpha: float, beta: float, gamma: float) -> float:
    """Return Score(a) = alpha*U(a) + beta*D(a) + gamma*C(a)."""

    validate_anchor_weights(alpha, beta, gamma)
    if any(value < 0.0 or value > 1.0 for value in (uniqueness, depth, consistency)):
        raise ValueError("Anchor components must be normalized to [0, 1].")
    return alpha * uniqueness + beta * depth + gamma * consistency


def score_anchor_candidates(
    candidates: list[AnchorCandidate],
    *,
    alpha: float = DEFAULT_ALPHA,
    beta: float = DEFAULT_BETA,
    gamma: float = DEFAULT_GAMMA,
) -> list[AnchorCandidate]:
    """Apply the proposed score to already-derived graph candidates."""

    validate_anchor_weights(alpha, beta, gamma)
    for candidate in candidates:
        candidate.alpha = alpha
        candidate.beta = beta
        candidate.gamma = gamma
        candidate.final_score = compute_anchor_score(
            candidate.uniqueness, candidate.depth, candidate.consistency_score, alpha, beta, gamma
        )
        candidate.informativeness = candidate.final_score
        candidate.selected = False
    return candidates


def build_graph_anchor_candidates(
    graph: ProposedGraph,
    *,
    alpha: float = DEFAULT_ALPHA,
    beta: float = DEFAULT_BETA,
    gamma: float = DEFAULT_GAMMA,
) -> list[AnchorCandidate]:
    """Derive U, D, and C directly from graph nodes and HaplotypePath walks.

    U(a) = 1 - node path coverage. D(a) is node path occurrences divided by
    the maximum node occurrence count in this graph. C(a) is Step 2's
    graph-derived node consistency C(v).
    """

    occurrence_counts = Counter(node_id for path in graph.paths for node_id in path.node_ids)
    max_occurrence = max(occurrence_counts.values(), default=1)
    candidates = [
        AnchorCandidate(
            anchor_id=f"a{index + 1}",
            node_id=node.node_id,
            minimizer=node.sequence,
            support=node.support,
            mean_position=node.mean_position,
            consistency_score=node.consistency_score,
            informativeness=0.0,
            haplotypes=node.haplotypes,
            uniqueness=1.0 - node.path_coverage,
            depth=occurrence_counts[node.node_id] / max_occurrence,
        )
        for index, node in enumerate(graph.nodes)
    ]
    return score_anchor_candidates(candidates, alpha=alpha, beta=beta, gamma=gamma)


def anchor_summary(candidates: list[AnchorCandidate]) -> dict:
    def mean(values: list[float]) -> float:
        return sum(values) / len(values) if values else 0.0

    selected = [candidate for candidate in candidates if candidate.selected]
    return {
        "formula": "Score(a) = alpha*U(a) + beta*D(a) + gamma*C(a)",
        "uniqueness_formula": "U(a) = 1 - node_path_coverage",
        "depth_formula": "D(a) = node_path_occurrences / maximum_node_path_occurrences",
        "consistency_source": "Step 2 graph-derived node consistency C(v)",
        "weights": {
            "alpha": candidates[0].alpha if candidates else DEFAULT_ALPHA,
            "beta": candidates[0].beta if candidates else DEFAULT_BETA,
            "gamma": candidates[0].gamma if candidates else DEFAULT_GAMMA,
        },
        "candidate_anchor_count": len(candidates),
        "selected_anchor_count": len(selected),
        "mean_uniqueness": mean([candidate.uniqueness for candidate in candidates]),
        "mean_depth": mean([candidate.depth for candidate in candidates]),
        "mean_consistency": mean([candidate.consistency_score for candidate in candidates]),
        "mean_final_score": mean([candidate.final_score for candidate in candidates]),
        "selected_anchor_mean_score": mean([candidate.final_score for candidate in selected]),
    }


def select_informative_anchors(
    candidates: list[AnchorCandidate],
    max_anchors: int,
    min_spacing: float,
) -> list[AnchorCandidate]:
    """Select highest-scoring candidates, using spacing only as a coexistence rule."""

    for candidate in candidates:
        candidate.selected = False
    ranked = sorted(
        candidates,
        key=lambda item: (-item.final_score, -item.consistency_score, -item.depth, item.node_id, item.anchor_id),
    )
    selected: list[AnchorCandidate] = []
    for candidate in ranked:
        if any(abs(candidate.mean_position - accepted.mean_position) < min_spacing for accepted in selected):
            continue
        candidate.selected = True
        selected.append(candidate)
        if len(selected) >= max_anchors:
            break
    return sorted(selected, key=lambda item: (item.mean_position, item.node_id, item.anchor_id))


def write_anchor_reports(candidates: list[AnchorCandidate], output_dir: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    scores_path = output_dir / "anchor_scores.csv"
    summary_path = output_dir / "anchor_summary.json"
    with scores_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "anchor_id", "node_id", "haplotype_count", "U", "D", "C",
            "alpha", "beta", "gamma", "final_score", "selected",
        ])
        writer.writeheader()
        for candidate in candidates:
            writer.writerow({
                "anchor_id": candidate.anchor_id,
                "node_id": candidate.node_id,
                "haplotype_count": len(candidate.haplotypes),
                "U": candidate.uniqueness,
                "D": candidate.depth,
                "C": candidate.consistency_score,
                "alpha": candidate.alpha,
                "beta": candidate.beta,
                "gamma": candidate.gamma,
                "final_score": candidate.final_score,
                "selected": candidate.selected,
            })
    summary_path.write_text(json.dumps(anchor_summary(candidates), indent=2) + "\n", encoding="utf-8")
    return {
        "anchor_scores_csv": str(scores_path.resolve()),
        "anchor_summary_json": str(summary_path.resolve()),
    }
