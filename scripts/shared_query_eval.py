"""Shared held-out-query comparison of the proposed graph and the minigraph graph.

The proposed graph stores minimizer k-mers as nodes (GFA links are 0M and carry no
intervening sequence), so reads cannot be aligned to it and alignment accuracy is
not comparable. This evaluation therefore uses one alignment-free measure on both:

  graphs     built from the same sequences (reference + train haplotypes):
             minigraph  graphs/train/<gene>.combined.gfa (existing build)
             proposed   build_proposed_pangenome_graph(...) at a fixed configuration
             linear     the reference alone (control)
  queries    the held-out test haplotypes (not used to build either graph)
  scheme     the proposed minimizer engine, same k and window, forward strand

  node recall        fraction of the query's minimizer walk present in the graph
                     (proposed: as a node; minigraph/linear: as a spelled k-mer)
  transition recall  fraction of consecutive query minimizer pairs supported
                     proposed: an edge between the two nodes exists (no check of the
                               bases in between, so this is weaker evidence)
                     minigraph/linear: some graph walk spells the exact query
                               substring from the first minimizer to the second

This measures how well each graph represents unseen haplotypes, not alignment
accuracy; alignment accuracy for minigraph/PanAligner is in outputs/truth_eval.
"""

from __future__ import annotations

import argparse
import csv
import platform
import statistics
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from common import METADATA_DIR, OUTPUTS_DIR, ROOT, ensure_dir, project_relative, write_json
from proposed.graph_constructor import build_proposed_pangenome_graph
from proposed.minimizer_engine import compute_minimizers
from proposed.models import AssemblyInput
from proposed.workflow import load_assembly_inputs
from simulate_truth_reads import read_fasta, sha256_file


GENES = ("APP", "PSEN1", "PSEN2")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


class SequenceGraph:
    """Oriented sequence graph from GFA S/L lines (0M overlaps only), for substring queries."""

    def __init__(self, segments: dict[str, str], links: list[tuple[str, str, str, str]]):
        if not segments:
            raise ValueError("graph has no segments")
        self.segments = segments
        self.successors: dict[tuple[str, str], set[tuple[str, str]]] = defaultdict(set)
        flip = {"+": "-", "-": "+"}
        for source, source_orient, target, target_orient in links:
            if source not in segments or target not in segments:
                raise ValueError(f"link {source}->{target} references a missing segment")
            self.successors[(source, source_orient)].add((target, target_orient))
            self.successors[(target, flip[target_orient])].add((source, flip[source_orient]))
        self._oriented: dict[tuple[str, str], str] = {}
        self._continuations: dict[tuple[tuple[str, str], int], frozenset[str]] = {}
        self._index: dict[int, set[str]] = {}

    @classmethod
    def from_gfa(cls, path: Path) -> "SequenceGraph":
        segments: dict[str, str] = {}
        links: list[tuple[str, str, str, str]] = []
        with Path(path).open("r", encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if fields[0] == "S":
                    segments[fields[1]] = fields[2].upper()
                elif fields[0] == "L":
                    if fields[5] not in ("0M", "*"):
                        raise ValueError(f"unsupported link overlap {fields[5]} in {path}")
                    links.append((fields[1], fields[2], fields[3], fields[4]))
        return cls(segments, links)

    @classmethod
    def linear(cls, sequence: str) -> "SequenceGraph":
        return cls({"ref": sequence.upper()}, [])

    def total_bases(self) -> int:
        return sum(len(sequence) for sequence in self.segments.values())

    def _sequence(self, state: tuple[str, str]) -> str:
        if state not in self._oriented:
            raw = self.segments[state[0]]
            self._oriented[state] = raw if state[1] == "+" else reverse_complement(raw)
        return self._oriented[state]

    def _prefixes(self, state: tuple[str, str], length: int) -> frozenset[str]:
        """All strings of `length` spelled by walks that start at the beginning of `state`."""
        key = (state, length)
        if key not in self._continuations:
            sequence = self._sequence(state)
            if len(sequence) >= length:
                result = frozenset({sequence[:length]})
            else:
                result = frozenset(
                    sequence + tail
                    for successor in self.successors.get(state, ())
                    for tail in self._prefixes(successor, length - len(sequence))
                )
            self._continuations[key] = result
        return self._continuations[key]

    def substrings(self, length: int) -> set[str]:
        """Every string of `length` spelled by a walk starting on a forward segment."""
        if length <= 0:
            raise ValueError("length must be positive")
        if length not in self._index:
            found: set[str] = set()
            for name in self.segments:
                state = (name, "+")
                sequence = self._sequence(state)
                for start in range(len(sequence)):
                    if start + length <= len(sequence):
                        found.add(sequence[start:start + length])
                    else:
                        head = sequence[start:]
                        for successor in self.successors.get(state, ()):
                            found.update(head + tail for tail in self._prefixes(successor, length - len(head)))
            self._index[length] = found
        return self._index[length]

    def contains(self, text: str) -> bool:
        return text.upper() in self.substrings(len(text))


@dataclass(frozen=True)
class WalkStep:
    minimizer: str
    position: int


def minimizer_walk(sequence: str, k: int, window_size: int, backend: str = "cpu") -> list[WalkStep]:
    """Query minimizers in order, with consecutive repeats of one k-mer collapsed as in the proposed graph."""
    records = compute_minimizers(sequence, k, window_size, "query", backend=backend)["records"]
    walk: list[WalkStep] = []
    for record in records:
        if not walk or walk[-1].minimizer != record["minimizer"]:
            walk.append(WalkStep(record["minimizer"], record["position"]))
    return walk


def score_sequence_graph(walk: list[WalkStep], query: str, graph: SequenceGraph, k: int) -> dict:
    node_hits = sum(graph.contains(step.minimizer) for step in walk)
    pairs = list(zip(walk, walk[1:]))
    transition_hits = sum(graph.contains(query[left.position:right.position + k]) for left, right in pairs)
    return _scores(len(walk), node_hits, len(pairs), transition_hits)


def score_proposed_graph(walk: list[WalkStep], node_by_sequence: dict[str, str], edges: set[tuple[str, str]]) -> dict:
    node_hits = sum(step.minimizer in node_by_sequence for step in walk)
    pairs = list(zip(walk, walk[1:]))
    transition_hits = sum(
        (node_by_sequence.get(left.minimizer), node_by_sequence.get(right.minimizer)) in edges for left, right in pairs
    )
    return _scores(len(walk), node_hits, len(pairs), transition_hits)


def _scores(nodes: int, node_hits: int, transitions: int, transition_hits: int) -> dict:
    return {
        "query_minimizers": nodes, "node_hits": node_hits,
        "node_recall": node_hits / nodes if nodes else 0.0,
        "query_transitions": transitions, "transition_hits": transition_hits,
        "transition_recall": transition_hits / transitions if transitions else 0.0,
    }


def graph_inputs(gene: str, train_manifest: Path) -> tuple[list[AssemblyInput], str, str]:
    """Reference + train haplotypes: the sequences minigraph used for <gene>.combined.gfa."""
    reference_name, reference = read_fasta(ROOT / "data" / "raw" / f"{gene.lower()}_reference.fa")[0]
    train = load_assembly_inputs(train_manifest, gene=gene)
    if not train:
        raise ValueError(f"no training assemblies for {gene} in {train_manifest}")
    return [AssemblyInput("REFERENCE", gene, "REFERENCE", reference), *train], reference_name, reference


def evaluate_gene(gene: str, args: argparse.Namespace) -> tuple[list[dict], list[dict]]:
    assemblies, _, reference = graph_inputs(gene, args.train_manifest)
    if args.query_set == "heldout":
        queries = load_assembly_inputs(args.test_manifest, gene=gene)
        overlap = {item.assembly_id for item in queries} & {item.assembly_id for item in assemblies}
        if overlap:
            raise ValueError(f"held-out queries also used to build the graphs: {sorted(overlap)}")
    else:  # in-sample positive control: graphs should contain these, the linear reference should not
        queries = [item for item in assemblies if item.assembly_id != "REFERENCE"]
    if not queries:
        raise ValueError(f"no {args.query_set} queries for {gene}")

    started = time.perf_counter()
    proposed = build_proposed_pangenome_graph(
        assemblies, k=args.k, window_size=args.window_size, min_node_support=args.min_node_support, backend=args.backend,
    )
    proposed_seconds = time.perf_counter() - started
    node_by_sequence = {node.sequence: node.node_id for node in proposed.nodes}
    proposed_edges = {(edge.source, edge.target) for edge in proposed.edges}

    minigraph_gfa = ROOT / "graphs" / "train" / f"{gene.lower()}.combined.gfa"
    minigraph = SequenceGraph.from_gfa(minigraph_gfa)
    linear = SequenceGraph.linear(reference)

    construction = [
        {"gene": gene, "graph": "proposed", "inputs": len(assemblies), "nodes": len(proposed.nodes), "edges": len(proposed.edges),
         "stored_bases": sum(len(node.sequence) for node in proposed.nodes), "stores_full_sequence": False,
         "build_seconds": round(proposed_seconds, 3), "backend_used": proposed.backend_used, "source": "built here"},
        {"gene": gene, "graph": "minigraph", "inputs": len(assemblies), "nodes": len(minigraph.segments),
         "edges": sum(len(targets) for targets in minigraph.successors.values()) // 2,
         "stored_bases": minigraph.total_bases(), "stores_full_sequence": True,
         "build_seconds": "", "backend_used": "cpu", "source": project_relative(minigraph_gfa)},
        {"gene": gene, "graph": "linear", "inputs": 1, "nodes": 1, "edges": 0, "stored_bases": linear.total_bases(),
         "stores_full_sequence": True, "build_seconds": "", "backend_used": "", "source": "reference FASTA"},
    ]

    rows = []
    for query in queries:
        walk = minimizer_walk(query.sequence, args.k, args.window_size, args.backend)
        results = {
            "proposed": score_proposed_graph(walk, node_by_sequence, proposed_edges),
            "minigraph": score_sequence_graph(walk, query.sequence, minigraph, args.k),
            "linear": score_sequence_graph(walk, query.sequence, linear, args.k),
        }
        for graph_name, scores in results.items():
            rows.append({"gene": gene, "query": query.assembly_id, "label": query.haplotype_label, "graph": graph_name, **scores})
        print(f"{gene} {query.assembly_id}: " + ", ".join(
            f"{name} node {scores['node_recall']:.5f} trans {scores['transition_recall']:.5f}" for name, scores in results.items()
        ), flush=True)
    return rows, construction


def aggregate(rows: list[dict]) -> list[dict]:
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in rows:
        grouped[(row["gene"], row["graph"])].append(row)
        grouped[("ALL", row["graph"])].append(row)
    summary = []
    for (gene, graph), members in sorted(grouped.items()):
        entry = {"gene": gene, "graph": graph, "queries": len(members)}
        for metric in ("node", "transition"):
            hits = sum(m[f"{metric}_hits"] for m in members)
            total = sum(m[f"query_{metric}s" if metric == "transition" else "query_minimizers"] for m in members)
            per_query = [m[f"{metric}_recall"] for m in members]
            entry[f"{metric}_recall_pooled"] = hits / total if total else 0.0
            entry[f"{metric}_recall_query_mean"] = statistics.mean(per_query)
            entry[f"{metric}_recall_query_min"] = min(per_query)
            entry[f"{metric}_misses"] = total - hits
        summary.append(entry)
    return summary


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Compare the proposed and minigraph graphs on shared held-out queries.")
    parser.add_argument("--genes", nargs="+", choices=GENES, default=list(GENES))
    parser.add_argument("--k", type=int, default=11, help="Minimizer k (proposed final-evaluation best configuration: 11).")
    parser.add_argument("--window-size", type=int, default=6, help="Minimizer window (best configuration: 6).")
    parser.add_argument("--min-node-support", type=int, default=2)
    parser.add_argument("--backend", choices=("cpu", "gpu", "auto"), default="cpu", help="Minimizer backend for graph and queries.")
    parser.add_argument("--query-set", choices=("heldout", "train"), default="heldout",
                        help="heldout: test haplotypes (the comparison); train: in-sample positive control.")
    parser.add_argument("--train-manifest", type=Path, default=METADATA_DIR / "train_manifest.json")
    parser.add_argument("--test-manifest", type=Path, default=METADATA_DIR / "test_manifest.json")
    parser.add_argument("--output-dir", type=Path, default=OUTPUTS_DIR / "shared_query_eval")
    args = parser.parse_args(argv)
    if args.k > 55 and args.backend != "cpu":
        parser.error("the CUDA kernel supports k <= 55")

    output_dir = ensure_dir(args.output_dir.resolve() / args.query_set)
    started = time.perf_counter()
    rows: list[dict] = []
    construction: list[dict] = []
    for gene in args.genes:
        gene_rows, gene_construction = evaluate_gene(gene, args)
        rows.extend(gene_rows)
        construction.extend(gene_construction)
    write_csv(output_dir / "per_query.csv", rows)
    write_csv(output_dir / "by_gene.csv", aggregate(rows))
    write_csv(output_dir / "graph_construction.csv", construction)
    write_json(output_dir / "run_manifest.json", {
        "command": "python " + " ".join([project_relative(Path(__file__)), *(argv if argv is not None else sys.argv[1:])]),
        "parameters": {"k": args.k, "window_size": args.window_size, "min_node_support": args.min_node_support, "backend": args.backend},
        "graph_inputs": "reference + train haplotypes (train_manifest), matching graphs/train/<gene>.combined.gfa",
        "queries": "held-out haplotypes (test_manifest)" if args.query_set == "heldout"
                   else "train haplotypes (in-sample positive control; not an evaluation of generalisation)",
        "input_sha256": {
            project_relative(path): sha256_file(path)
            for path in [args.train_manifest.resolve(), args.test_manifest.resolve(),
                         *(ROOT / "graphs" / "train" / f"{gene.lower()}.combined.gfa" for gene in args.genes)]
        },
        "host": platform.platform(),
        "python": platform.python_version(),
        "runtime_seconds": round(time.perf_counter() - started, 3),
        "comparison_type": "graph representation of unseen haplotypes (alignment-free); not alignment accuracy",
        "notes": [
            "Proposed transition recall checks only that an edge exists; minigraph/linear transition recall checks the exact bases.",
            "Node recall: the proposed graph stores only k-mers selected as minimizers in its inputs (low-support ones are "
            "kept as rescued nodes), while minigraph/linear are checked for the k-mer anywhere in their sequence.",
            "Forward strand only: the proposed minimizers are not canonical.",
            "minigraph build time is not measured here; graphs/train/<gene>.combined.gfa is prebuilt by scripts/build_graph.py.",
        ],
    })
    print(f"wrote {project_relative(output_dir)}")


if __name__ == "__main__":
    main()
