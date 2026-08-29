# 14. Proposed Methodology Incremental Implementation

## Purpose

This document explains the newly added **proposed methodology** that now coexists with the existing PanAligner-centered workflow.

It is intentionally implemented as a **modular research prototype**, not as a replacement for:

- the current `main.py` pipeline
- the real `minigraph` graph construction workflow
- the real `PanAligner` alignment workflow
- the educational theory/demo modules already present in the repository

## Existing architecture we preserved

Before adding the new methodology, the working project already had these active components:

1. FASTA preprocessing
2. reproducible train/test split
3. `minigraph`-based train graph construction
4. graph visualization and graph statistics
5. educational theory demonstrations
6. held-out evaluation using real alignment
7. JSON, TXT, CSV, and PNG reporting
8. Tkinter GUI support

These remain intact.

## Integration strategy

To avoid breaking the stable workflow, the new methodology was added as a **parallel component stack** under:

- `scripts/proposed/`
- `proposed_methodology.py`

This means:

- the current PanAligner pipeline still runs exactly as before
- the GUI remains focused on the stable existing workflow
- the new proposed implementation can be tested independently
- the new modules can later be integrated more deeply only if they prove useful

## What was implemented

### 1. GPU-ready k-mer / minimizer engine

Implemented in:

- `scripts/proposed/minimizer_engine.py`

Current behavior:

- computes deterministic k-mers and minimizers
- exposes backend selection: `auto`, `cpu`, or `gpu`
- honestly reports which backend was actually used
- falls back to CPU when a GPU backend such as CuPy is unavailable

Important honesty note:

The current implementation is a **GPU-ready backend interface with deterministic CPU fallback**, not a full CUDA-optimized production kernel.

### 2. Haplotype-consistency scorer

Implemented in:

- `scripts/proposed/haplotype_scoring.py`

Current behavior:

- compares minimizer sketches between assemblies
- scores shared minimizer overlap
- scores order preservation
- applies a spacing-based similarity factor
- produces a bounded pairwise consistency score

This gives the proposed graph builder a way to prefer stable shared structure across related assemblies.

### 3. Alignment-free graph constructor

Implemented in:

- `scripts/proposed/graph_constructor.py`

Current behavior:

- consumes assemblies directly from a dataset manifest
- computes minimizers for each assembly
- merges shared minimizers into graph nodes
- creates directed edges from within-assembly minimizer order
- stores support counts per node and edge
- computes node and edge consistency statistics
- exports the proposed graph to:
  - JSON
  - GFA

Important distinction:

This graph is **our alignment-free proposed graph abstraction**. It is not the same thing as the repository's existing `minigraph` construction pipeline.

### 4. Intelligent anchor selection

Implemented in:

- `scripts/proposed/anchor_selection.py`

Current behavior:

- ranks candidate anchors by support and consistency
- computes an informativeness score
- prevents overly dense nearby anchor choices using spacing constraints
- returns a compact selected anchor set for downstream use

### 5. PG-SCUnK-style quality feedback controller

Implemented in:

- `scripts/proposed/quality_feedback.py`

Current behavior:

- evaluates proposed graph quality after construction
- measures:
  - node consistency
  - edge consistency
  - support ratio
  - anchor coverage
  - branching ratio
  - overall quality score
- decides whether the graph meets threshold
- retunes parameters if quality is insufficient
- reconstructs the graph using updated parameters
- records each round's decision

Important honesty note:

This is an **in-loop quality-feedback controller inspired by the proposal**, not a claim of reproducing any external tool's exact internals.

### 6. Workflow runner

Implemented in:

- `scripts/proposed/workflow.py`
- `proposed_methodology.py`

Current behavior:

- loads assemblies from either:
  - the existing project manifests such as `data/metadata/train_manifest.json`
  - a small synthetic assembly manifest for testing
- runs the proposed graph workflow per gene
- writes separate artifacts under `outputs/proposed/` or a user-specified output directory
- writes summary JSON describing quality, chosen parameters, and generated files

## Output artifacts produced by the proposed stack

For each gene, the workflow writes artifacts like:

- `*.proposed_graph.json`
- `*.proposed_graph.gfa`
- `*.proposed_graph.png`
- `*.proposed_graph.stats.json`
- `*.selected_anchors.json`
- `*.quality.json`
- `*.pg_scunk_feedback.json`

And a run-level summary:

- `proposed_methodology_summary.json`

## Tests added

New unit tests were added under:

- `tests/test_proposed_minimizer_engine.py`
- `tests/test_proposed_haplotype_scoring.py`
- `tests/test_proposed_anchor_selection.py`
- `tests/test_proposed_graph_constructor.py`
- `tests/test_proposed_quality_feedback.py`

Synthetic test data was added under:

- `tests/data/proposed/`

This includes:

- small toy FASTA files
- a toy manifest for end-to-end proposed-methodology checks

## Validation performed

After implementation, the following checks were run successfully:

- `python -m unittest discover -s tests -p "test_*.py"`
- `python main.py --help`
- `python proposed_methodology.py --help`
- import smoke check for `main.py` and `gui_app.py`

A real-data smoke check was also run against the existing `train_manifest.json` on a small `PSEN2` slice using the proposed PG-SCUnK loop. It completed successfully after caching pairwise consistency values inside the proposed graph constructor.

## Current limitations

### GPU limitation

The GPU contribution is currently an incremental backend-ready implementation. It does not yet include a full low-level CUDA or ROCm production kernel.

### Biological limitation

The current project dataset is not a true whole-genome haplotype-resolved benchmark suite. In this implementation, available sequence records are treated as assembly-like haplotype inputs for research prototyping.

### Graph-construction limitation

The proposed graph builder is a lightweight alignment-free abstraction based on minimizer sharing and order, not a replacement for established production graph constructors.

### Quality-feedback limitation

The PG-SCUnK controller currently retunes a practical subset of parameters:

- `k`
- `window_size`
- `min_node_support`
- `min_anchor_spacing`

This is enough to demonstrate the feedback-loop idea honestly, but it is not yet a full graph-optimization framework.

## Recommended current framing

The safest and most accurate way to present the project now is:

**The repository contains a stable PanAligner-centered reproduction workflow, plus a new modular prototype of our proposed alignment-free methodology with GPU-ready minimizer computation, haplotype-consistency-aware graph construction and anchoring, and PG-SCUnK-style quality feedback.**
