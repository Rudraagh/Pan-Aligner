from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from common import OUTPUTS_DIR, ensure_dir, read_json, validate_dna, write_json
from visualize import visualize_graph

from .anchor_selection import DEFAULT_ALPHA, DEFAULT_BETA, DEFAULT_GAMMA
from .graph_constructor import write_proposed_graph_bundle
from .models import AssemblyInput
from .quality_feedback import (
    DEFAULT_IMPROVEMENT_TOLERANCE,
    DEFAULT_STAGNATION_PATIENCE,
    run_pg_scunk_feedback_loop,
    write_pg_scunk_reports,
)


def _read_single_fasta_sequence(path: Path) -> str:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"FASTA file is empty: {path}")
    sequence = "".join(line for line in lines if not line.startswith(">")).upper()
    if not validate_dna(sequence):
        raise ValueError(f"Invalid DNA sequence found in {path}")
    return sequence


def load_assembly_inputs(
    dataset_manifest_path: Path,
    *,
    gene: str | None = None,
    max_assemblies_per_bucket: int | None = None,
) -> list[AssemblyInput]:
    payload = read_json(dataset_manifest_path)

    if "assemblies" in payload:
        assemblies = []
        for item in payload["assemblies"]:
            if gene and item["gene"].upper() != gene.upper():
                continue
            assemblies.append(
                AssemblyInput(
                    assembly_id=item["id"],
                    gene=item["gene"].upper(),
                    haplotype_label=item.get("label", "UNKNOWN").upper(),
                    sequence=_read_single_fasta_sequence(Path(item["path"])),
                )
            )
        return assemblies

    assemblies: list[AssemblyInput] = []
    for gene_name, gene_payload in payload.get("genes", {}).items():
        if gene and gene_name.upper() != gene.upper():
            continue
        for bucket_name in ("healthy", "unhealthy"):
            entries = gene_payload.get(bucket_name, [])
            if max_assemblies_per_bucket is not None:
                entries = entries[:max_assemblies_per_bucket]
            for entry in entries:
                assemblies.append(
                    AssemblyInput(
                        assembly_id=entry["id"],
                        gene=gene_name.upper(),
                        haplotype_label=bucket_name.upper(),
                        sequence=_read_single_fasta_sequence(Path(entry["path"])),
                    )
                )
    return assemblies


def run_proposed_methodology(
    dataset_manifest_path: Path,
    *,
    output_dir: Path | None = None,
    gene: str | None = None,
    backend: str = "auto",
    max_assemblies_per_bucket: int | None = 8,
    quality_threshold: float = 0.72,
    max_rounds: int = 3,
    improvement_tolerance: float = DEFAULT_IMPROVEMENT_TOLERANCE,
    stagnation_patience: int = DEFAULT_STAGNATION_PATIENCE,
    anchor_alpha: float = DEFAULT_ALPHA,
    anchor_beta: float = DEFAULT_BETA,
    anchor_gamma: float = DEFAULT_GAMMA,
) -> dict:
    assemblies = load_assembly_inputs(
        dataset_manifest_path,
        gene=gene,
        max_assemblies_per_bucket=max_assemblies_per_bucket,
    )
    if not assemblies:
        raise ValueError("No assemblies matched the requested manifest and gene filters.")

    grouped: dict[str, list[AssemblyInput]] = {}
    for assembly in assemblies:
        grouped.setdefault(assembly.gene, []).append(assembly)

    root_output_dir = ensure_dir(output_dir or (OUTPUTS_DIR / "proposed"))
    summary = {
        "dataset_manifest": str(dataset_manifest_path.resolve()),
        "backend_requested": backend,
        "genes": {},
    }

    for gene_name, gene_assemblies in sorted(grouped.items()):
        gene_output_dir = ensure_dir(root_output_dir / gene_name.lower())
        feedback = run_pg_scunk_feedback_loop(
            gene_assemblies,
            backend=backend,
            quality_threshold=quality_threshold,
            max_rounds=max_rounds,
            improvement_tolerance=improvement_tolerance,
            stagnation_patience=stagnation_patience,
            base_parameters={
                "anchor_alpha": anchor_alpha,
                "anchor_beta": anchor_beta,
                "anchor_gamma": anchor_gamma,
            },
        )
        graph = feedback.pop("graph")
        bundle = write_proposed_graph_bundle(graph, gene_output_dir)
        visualize_graph(Path(bundle["graph_gfa"]), gene_output_dir / f"{gene_name.lower()}.proposed_graph")
        feedback_reports = write_pg_scunk_reports(feedback, gene_output_dir)
        feedback_json = gene_output_dir / f"{gene_name.lower()}.pg_scunk_feedback.json"
        write_json(
            feedback_json,
            {
                key: feedback[key]
                for key in (
                    "accepted", "threshold_reached", "initial_quality", "best_quality", "best_overall_quality", "final_quality",
                    "quality_threshold", "min_anchor_coverage", "improvement_tolerance", "stagnation_patience",
                    "stagnation_reason", "best_round", "total_rounds", "unique_configurations_evaluated",
                    "best_configuration", "final_configuration", "best_evaluation", "final_evaluation",
                    "evaluations", "final_parameters", "round_summaries",
                )
            },
        )
        summary["genes"][gene_name] = {
            "assembly_count": len(gene_assemblies),
            "graph_quality": graph.quality,
            "final_parameters": graph.parameters,
            "backend_used": graph.backend_used,
            "accepted": feedback["accepted"],
            **bundle,
            **feedback_reports,
            "feedback_json": str(feedback_json.resolve()),
        }

    summary_path = root_output_dir / "proposed_methodology_summary.json"
    write_json(summary_path, summary)
    summary["summary_json"] = str(summary_path.resolve())
    return summary
