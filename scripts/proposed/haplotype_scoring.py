from __future__ import annotations

from .models import MinimizerRecord


def _minimizer_positions(records: list[MinimizerRecord]) -> dict[str, list[int]]:
    positions: dict[str, list[int]] = {}
    for record in records:
        positions.setdefault(record.minimizer, []).append(record.position)
    return positions


def legacy_pairwise_consistency(left: list[MinimizerRecord], right: list[MinimizerRecord]) -> float:
    """Legacy minimizer-sketch comparator retained for experimental baselines only."""
    if not left or not right:
        return 0.0

    left_positions = _minimizer_positions(left)
    right_positions = _minimizer_positions(right)
    shared = sorted(set(left_positions).intersection(right_positions))
    if not shared:
        return 0.0

    overlap_ratio = len(shared) / max(len(set(left_positions)), len(set(right_positions)))

    ordered_matches = 0
    order_checks = 0
    for first, second in zip(shared, shared[1:]):
        left_delta = min(left_positions[second]) - min(left_positions[first])
        right_delta = min(right_positions[second]) - min(right_positions[first])
        order_checks += 1
        if left_delta >= 0 and right_delta >= 0:
            ordered_matches += 1

    order_ratio = 1.0 if order_checks == 0 else ordered_matches / order_checks
    spacing_penalty = 0.0
    for minimizer in shared:
        spacing_penalty += abs(min(left_positions[minimizer]) - min(right_positions[minimizer]))
    normalized_spacing = spacing_penalty / max(1, len(shared))
    spacing_score = 1.0 / (1.0 + normalized_spacing)
    return 0.5 * overlap_ratio + 0.3 * order_ratio + 0.2 * spacing_score


def pairwise_haplotype_consistency(left: list[MinimizerRecord], right: list[MinimizerRecord]) -> float:
    """Backward-compatible alias for the isolated legacy comparator."""

    return legacy_pairwise_consistency(left, right)
