from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class AssemblyInput:
    assembly_id: str
    gene: str
    haplotype_label: str
    sequence: str


@dataclass(frozen=True)
class MinimizerRecord:
    minimizer: str
    hash_value: int
    position: int
    window_start: int
    source_assembly: str


@dataclass
class AnchorCandidate:
    anchor_id: str
    minimizer: str
    support: int
    mean_position: float
    consistency_score: float
    informativeness: float
    haplotypes: list[str] = field(default_factory=list)
    node_id: str = ""
    uniqueness: float = 0.0
    depth: float = 0.0
    alpha: float = 0.0
    beta: float = 0.0
    gamma: float = 0.0
    final_score: float = 0.0
    selected: bool = False


@dataclass
class ProposedGraphNode:
    node_id: str
    sequence: str
    support: int
    mean_position: float
    consistency_score: float
    haplotypes: list[str] = field(default_factory=list)
    is_rescued: bool = False
    path_coverage: float = 0.0
    predecessor_consistency: float = 0.0
    successor_consistency: float = 0.0
    local_path_consistency: float = 0.0


@dataclass
class ProposedGraphEdge:
    source: str
    target: str
    support: int
    consistency_score: float
    haplotypes: list[str] = field(default_factory=list)
    path_coverage: float = 0.0
    transition_consistency: float = 0.0


@dataclass
class HaplotypePath:
    """Observed ordered minimizer walk for one input assembly."""

    assembly_id: str
    haplotype_label: str
    node_ids: list[str]
    positions: list[int]
    valid_transition_count: int = 0
    transition_fraction: float = 0.0
    mean_node_consistency: float = 0.0
    mean_edge_consistency: float = 0.0
    rescued_node_count: int = 0
    questionable_node_count: int = 0
    consistency_score: float = 0.0


@dataclass
class ProposedGraph:
    gene: str
    assembly_count: int
    parameters: dict
    backend_used: str
    nodes: list[ProposedGraphNode]
    edges: list[ProposedGraphEdge]
    anchors: list[AnchorCandidate]
    quality: dict
    paths: list[HaplotypePath] = field(default_factory=list)
    metadata: dict = field(default_factory=dict)
    anchor_candidates: list[AnchorCandidate] = field(default_factory=list)
