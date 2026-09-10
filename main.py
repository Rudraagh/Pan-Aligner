from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS_DIR = ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from build_graph import build_graphs
from common import METADATA_DIR, OUTPUTS_DIR, ROOT as PROJECT_ROOT, ensure_dir, resolve_binary
from paper_evaluation import evaluate_panaligner_workflow
from preprocess import preprocess_fasta_files
from predictor import analyze_custom_query, prepare_query_fasta
from split_dataset import split_dataset
from theory.anchors_demo import run_anchors_demo
from theory.chaining_demo import run_chaining_demo
from theory.dag_converter_demo import run_dag_converter_demo
from theory.graph_analysis_demo import run_graph_analysis
from theory.path_cover_demo import run_path_cover_demo
from theory.precedence_demo import run_precedence_demo
from theory.scc_demo import run_scc_demo
from visualize import visualize_graph


REPORTS_DIR = OUTPUTS_DIR / "reports"


def default_fasta_inputs() -> list[Path]:
    return sorted(PROJECT_ROOT.glob("*_combined.fasta"))


def build_training_artifacts(
    input_fastas: list[Path],
    minigraph_bin: Path,
    threads: int,
    split_seed: int,
    test_fraction: float,
) -> tuple[dict, Path]:
    preprocess_fasta_files(input_fastas)
    split_dataset(METADATA_DIR / "preprocess_manifest.json", test_fraction, split_seed)
    train_graph_dir = ensure_dir(PROJECT_ROOT / "graphs" / "train")
    train_graph_manifest_path = METADATA_DIR / "train_graph_manifest.json"
    graph_manifest = build_graphs(
        METADATA_DIR / "train_manifest.json",
        minigraph_bin,
        threads,
        train_graph_dir,
        train_graph_manifest_path,
    )

    graph_output_dir = ensure_dir(OUTPUTS_DIR / "graphs" / "train")
    for gene_name, graph_entry in graph_manifest["graphs"].items():
        visualize_graph(Path(graph_entry["combined_graph"]), graph_output_dir / f"{gene_name.lower()}.combined")
        visualize_graph(Path(graph_entry["healthy_graph"]), graph_output_dir / f"{gene_name.lower()}.healthy")
        visualize_graph(Path(graph_entry["unhealthy_graph"]), graph_output_dir / f"{gene_name.lower()}.unhealthy")
    return graph_manifest, train_graph_manifest_path


def choose_theory_graph(graph_manifest: dict) -> Path:
    if "APP" in graph_manifest["graphs"]:
        return Path(graph_manifest["graphs"]["APP"]["combined_graph"]).resolve()
    first_gene = sorted(graph_manifest["graphs"])[0]
    return Path(graph_manifest["graphs"][first_gene]["combined_graph"]).resolve()


def run_theory_suite(graph_path: Path) -> dict:
    graph_context = run_graph_analysis(graph_path)
    scc_result = run_scc_demo(graph_path)
    dag_result = run_dag_converter_demo(graph_path)
    path_cover_result = run_path_cover_demo(graph_path)
    anchors_result = run_anchors_demo(graph_path)
    precedence_result = run_precedence_demo(graph_path)
    chaining_result = run_chaining_demo(graph_path)
    return {
        "graph_context": graph_context,
        "scc_result": scc_result,
        "dag_result": dag_result,
        "path_cover_result": path_cover_result,
        "anchors_result": anchors_result,
        "precedence_result": precedence_result,
        "chaining_result": chaining_result,
    }


def run_alignment_evaluation(
    train_graph_manifest_path: Path,
    panaligner_bin: Path,
    minigraph_bin: Path | None,
    threads: int,
) -> dict:
    return evaluate_panaligner_workflow(
        METADATA_DIR / "test_manifest.json",
        train_graph_manifest_path,
        panaligner_bin,
        minigraph_bin,
        threads,
        OUTPUTS_DIR / "evaluation",
    )


def existing_evaluation_inputs_ready() -> bool:
    required = [
        METADATA_DIR / "test_manifest.json",
        METADATA_DIR / "train_graph_manifest.json",
    ]
    return all(path.exists() for path in required)


def write_theoretical_report(theory_result: dict) -> Path:
    ensure_dir(REPORTS_DIR)
    context = theory_result["graph_context"]
    scc_count = len(theory_result["scc_result"]["components"])
    path_count = len(theory_result["path_cover_result"]["paths"])
    chain = theory_result["chaining_result"]["best_chain"]
    lines = [
        "Theoretical Reproduction Report",
        "===============================",
        "This report documents the simplified educational reproductions of PanAligner paper concepts implemented alongside the real PanAligner integration.",
        "",
        "1. SCC logic",
        f"- Tarjan SCC detection was applied to the educational cyclic graph derived from {context.source_label}.",
        f"- SCC count in the demo graph: {scc_count}.",
        "",
        "2. DAG conversion",
        "- A DFS traversal identifies back edges and removes them only in the educational copy of the graph.",
        "- This mirrors the paper's idea of removing back edges inside SCCs to obtain a DAG-like structure for downstream analysis.",
        "",
        "3. Path cover concept",
        f"- A simplified greedy path cover was generated on the DAG approximation with {path_count} paths.",
        "- This is an educational analogue of the paper's width-aware path-cover preprocessing.",
        "",
        "4. Anchor representation",
        "- Anchors are represented as (vertex, [x..y], [c..d], weight).",
        "- The demo anchor set includes walk-consistent anchors plus distractor anchors to illustrate chaining decisions.",
        "",
        "5. Precedence relation",
        "- A precedes B when query order increases and graph reachability exists.",
        "- Same-vertex precedence is additionally allowed in cyclic situations in the simplified model.",
        "",
        "6. Chaining DP",
        "- The project includes a simplified dynamic programming recurrence for best chain score ending at each anchor.",
        f"- Best chain recovered in the current run: {' -> '.join(anchor.anchor_id for anchor in chain)}.",
        "",
        "7. Iterative convergence",
        "- An iterative score-propagation demonstration is included to mirror the paper's convergence-oriented algorithms conceptually.",
        "",
        "8. Relation to the PanAligner paper",
        "- The real PanAligner binary remains the production alignment engine.",
        "- The theory layer is a pedagogical reproduction of concepts from the paper, not a reimplementation of PanAligner's optimized internals.",
    ]
    report_path = REPORTS_DIR / "theoretical_reproduction_report.txt"
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path


def write_final_summary(theory_result: dict | None = None, evaluation_metrics: dict | None = None) -> Path:
    ensure_dir(REPORTS_DIR)
    lines = [
        "Final Project Summary",
        "=====================",
        "Architecture",
        "- Real PanAligner implementation for graph alignment.",
        "- Simplified theoretical reproductions for educational explanation of the paper.",
        "- Train/test graph-alignment evaluation pipeline with plots and reports.",
        "",
        "Implementation summary",
        "- FASTA preprocessing, stratified train/test split, graph construction, GAF parsing, theory demonstrations, evaluation, and plotting remain intact.",
        "- A new main.py now acts as the single project entry point.",
        "",
        "Reproduced paper components",
        "- SCC detection",
        "- DFS-based back-edge removal",
        "- Cyclic graph to DAG approximation",
        "- Path cover generation",
        "- Anchor representation",
        "- Precedence relation",
        "- Simplified co-linear chaining DP",
        "- Gap cost demonstration",
        "- Iterative chaining convergence demo",
        "- Reachability analysis",
    ]
    if evaluation_metrics:
        overall = evaluation_metrics["overall"]
        lines.extend(
            [
                "",
                "Evaluation pipeline",
                f"- query count: {overall['query_count']}",
                f"- aligned query count: {overall['aligned_query_count']}",
                f"- alignment rate: {overall['alignment_rate']:.4f}",
                f"- mean identity: {overall['mean_identity']:.4f}",
                f"- mean coverage: {overall['mean_coverage']:.4f}",
                f"- mean MAPQ: {overall['mean_mapq']:.4f}",
            ]
        )
    lines.extend(
        [
            "",
            "Strengths",
            "- The project now separates production alignment from theoretical understanding clearly.",
            "- The theory modules make the PanAligner paper easier to defend during viva without interfering with the real binary.",
            "",
            "Limitations",
            "- The theory modules are simplified educational reproductions and do not match PanAligner's optimized internal data structures or asymptotic performance.",
            "- Some current locus graphs are DAG-like, so cyclic behavior is demonstrated on a representative educational augmentation of the real graph.",
            "",
            "Future work",
            "- Add richer anchor generation from real seed hits.",
            "- Compare educational chaining outputs directly against PanAligner alignment traces for selected reads.",
            "- Extend evaluation to larger cyclic pangenome datasets closer to the full paper-scale setup.",
        ]
    )
    report_path = REPORTS_DIR / "final_project_summary.txt"
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path


def print_json_summary(payload: dict) -> None:
    print(json.dumps(payload, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description="Single entry point for the PanAligner paper reproduction project.")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--full-pipeline", action="store_true", help="Run preprocessing, graph construction, theory modules, PanAligner execution, evaluation, and report generation.")
    mode.add_argument("--theory-only", action="store_true", help="Run only the educational PanAligner theory reproductions.")
    mode.add_argument("--evaluate", action="store_true", help="Run the held-out PanAligner evaluation workflow on the current train/test split.")
    mode.add_argument("--custom-query-analysis", action="store_true", help="Analyze a custom query sequence or FASTA against the current graphs and generate score outputs.")
    parser.add_argument("--theory-graph", type=Path, default=None, help="Optional graph GFA for theory-only mode.")
    parser.add_argument("--input-fastas", nargs="+", type=Path, default=default_fasta_inputs(), help="Input FASTA files.")
    parser.add_argument("--threads", type=int, default=4, help="Thread count for minigraph and PanAligner.")
    parser.add_argument("--split-seed", type=int, default=42, help="Random seed for reproducible splitting.")
    parser.add_argument("--test-fraction", type=float, default=0.20, help="Fraction reserved for test data.")
    parser.add_argument("--minigraph-bin", type=str, default=None, help="Path to the minigraph binary.")
    parser.add_argument("--panaligner-bin", type=str, default=None, help="Path to the PanAligner binary.")
    parser.add_argument("--hybrid-approximation", action="store_true", help="Use local Minigraph-first heuristic screening before PanAligner during evaluation.")
    parser.add_argument("--graph-manifest", type=Path, default=METADATA_DIR / "graph_manifest.json", help="Graph manifest used by custom query analysis.")
    parser.add_argument("--query-fasta", type=Path, default=None, help="Custom query FASTA for custom-query-analysis mode.")
    parser.add_argument("--query-sequence", type=str, default=None, help="Inline DNA query for custom-query-analysis mode.")
    parser.add_argument("--query-argument", type=str, default=None, help="Optional note or argument to include in custom query outputs.")
    parser.add_argument("--gene", type=str, default=None, help="Optional gene override for custom query analysis.")
    parser.add_argument("--output-dir", type=Path, default=None, help="Output directory for custom query analysis artifacts.")
    args = parser.parse_args()

    if args.theory_only:
        theory_result = run_theory_suite(args.theory_graph.resolve() if args.theory_graph else None)
        theory_report = write_theoretical_report(theory_result)
        final_report = write_final_summary(theory_result=theory_result)
        print_json_summary(
            {
                "mode": "theory-only",
                "theory_report": str(theory_report.resolve()),
                "final_report": str(final_report.resolve()),
            }
        )
        return

    if args.custom_query_analysis:
        if args.query_fasta is None and args.query_sequence is None:
            raise ValueError("Pass --query-fasta or --query-sequence with --custom-query-analysis.")
        panaligner_bin = resolve_binary(args.panaligner_bin, [Path("PanAligner/PanAligner")])
        query_fasta = prepare_query_fasta(
            args.query_sequence,
            args.query_fasta,
            ROOT / "data" / "queries" / "custom_query.fa",
        )
        prediction = analyze_custom_query(
            query_fasta,
            args.graph_manifest.resolve(),
            panaligner_bin,
            args.threads,
            args.gene,
            args.output_dir.resolve() if args.output_dir else None,
            args.query_argument,
        )
        print_json_summary(
            {
                "mode": "custom-query-analysis",
                "query_fasta": prediction["query_fasta"],
                "selected_gene": prediction["selected_gene"],
                "prediction": prediction["prediction"],
                "confidence": prediction["confidence"],
                "alignment_detected": prediction["alignment_detected"],
                "score_plot": prediction["score_plot"],
                "prediction_json": str((Path(prediction["score_plot"]).parent / "prediction.json").resolve()),
            }
        )
        return

    panaligner_bin = resolve_binary(args.panaligner_bin, [Path("PanAligner/PanAligner")])
    minigraph_bin = resolve_binary(args.minigraph_bin, [Path("minigraph/minigraph")]) if args.hybrid_approximation else None
    theory_result = None
    evaluation_metrics = None
    graph_manifest = None
    train_graph_manifest_path = METADATA_DIR / "train_graph_manifest.json"

    if args.evaluate and existing_evaluation_inputs_ready():
        evaluation_metrics = run_alignment_evaluation(
            train_graph_manifest_path,
            panaligner_bin,
            minigraph_bin,
            args.threads,
        )
        final_report = write_final_summary(None, evaluation_metrics)
        print_json_summary(
            {
                "mode": "evaluate",
                "final_report": str(final_report.resolve()),
                "alignment_rate": evaluation_metrics["overall"]["alignment_rate"],
                "aligned_query_count": evaluation_metrics["overall"]["aligned_query_count"],
            }
        )
        return

    minigraph_bin = resolve_binary(args.minigraph_bin, [Path("minigraph/minigraph")])
    graph_manifest, train_graph_manifest_path = build_training_artifacts(
        [path.resolve() for path in args.input_fastas],
        minigraph_bin,
        args.threads,
        args.split_seed,
        args.test_fraction,
    )

    if args.full_pipeline:
        theory_result = run_theory_suite(choose_theory_graph(graph_manifest))
        evaluation_metrics = run_alignment_evaluation(
            train_graph_manifest_path,
            panaligner_bin,
            minigraph_bin,
            args.threads,
        )
        theory_report = write_theoretical_report(theory_result)
        final_report = write_final_summary(theory_result, evaluation_metrics)
        print_json_summary(
            {
                "mode": "full-pipeline",
                "theory_report": str(theory_report.resolve()),
                "final_report": str(final_report.resolve()),
                "alignment_rate": evaluation_metrics["overall"]["alignment_rate"],
            }
        )
        return

    if args.evaluate:
        evaluation_metrics = run_alignment_evaluation(
            train_graph_manifest_path,
            panaligner_bin,
            minigraph_bin,
            args.threads,
        )
        final_report = write_final_summary(None, evaluation_metrics)
        print_json_summary(
            {
                "mode": "evaluate",
                "final_report": str(final_report.resolve()),
                "alignment_rate": evaluation_metrics["overall"]["alignment_rate"],
                "aligned_query_count": evaluation_metrics["overall"]["aligned_query_count"],
            }
        )


if __name__ == "__main__":
    main()
