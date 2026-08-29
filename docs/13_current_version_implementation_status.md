# 13. Current Version Implementation Status

## Document purpose

This document summarizes everything implemented in the **current version** of the project as of **2026-08-29**.

It is intended to answer:

- what is actually working now
- which parts are part of the main workflow
- which outputs are already generated
- which modules remain in the repository as legacy or archival work

## Current project identity

The current version of the repository is best described as:

**Implementation and partial practical reproduction of the PanAligner paper, with an added educational theory layer and a desktop GUI for running and inspecting the workflow.**

This is no longer primarily framed as a disease-classification or ML-first project.

## What is fully implemented in the current version

### 1. Single main entry point

The project now has a unified top-level runner:

- `main.py`

Implemented run modes:

- `--full-pipeline`
- `--theory-only`
- `--evaluate`
- `--custom-query-analysis`

This means the project can now be run from one central script instead of requiring the user to execute each stage manually.

### 2. FASTA preprocessing

Implemented in:

- `scripts/preprocess.py`

Current capabilities:

- reads the root-level combined FASTA inputs
- normalizes sequence formatting
- validates DNA characters
- separates reference, healthy, and unhealthy records
- supports deduplication
- writes structured processed outputs
- writes `data/metadata/preprocess_manifest.json`

### 3. Reproducible train/test splitting

Implemented in:

- `scripts/split_dataset.py`

Current capabilities:

- reproducible splitting with configurable seed
- class-aware splitting inside each gene
- writes combined and per-sample FASTA files
- keeps manifest-based bookkeeping for later stages

Generated metadata includes:

- `data/metadata/train_manifest.json`
- `data/metadata/test_manifest.json`
- `data/metadata/split_manifest.json`
- `data/metadata/split_summary.txt`

### 4. Pangenome graph construction

Implemented in:

- `scripts/build_graph.py`

Current capabilities:

- builds graphs using the real `minigraph` workflow
- builds three graphs per gene:
  - healthy
  - unhealthy
  - combined
- stores train graphs under `graphs/train/`
- stores graph manifest metadata in `data/metadata/train_graph_manifest.json`

### 5. Graph parsing, statistics, and visualization

Implemented in:

- `scripts/visualize.py`

Current capabilities:

- parses GFA graphs into directed graph form
- computes graph statistics
- detects whether cycles exist
- renders graph PNGs
- supports alignment-aware graph highlighting

Generated outputs include:

- `outputs/graphs/*.png`
- `outputs/graphs/*.stats.json`
- `outputs/graphs/train/*.png`
- `outputs/graphs/train/*.stats.json`

### 6. Real PanAligner alignment wrapper

Implemented in:

- `scripts/align.py`

Current capabilities:

- runs the real `PanAligner` binary
- writes GAF alignment files
- also supports local `minigraph` alignment runs for hybrid evaluation logic

This is a real tool integration, not a mock implementation.

### 7. GAF parsing and best-alignment selection

Implemented in:

- `scripts/parse_gaf.py`

Current capabilities:

- parses GAF lines into structured alignment objects
- extracts alignment score, normalized score, identity, coverage, MAPQ, path span, and traversed nodes
- selects the best alignment using score-driven ranking

### 8. Educational theory reproduction layer

Implemented in:

- `scripts/theory/common.py`
- `scripts/theory/graph_analysis_demo.py`
- `scripts/theory/scc_demo.py`
- `scripts/theory/dag_converter_demo.py`
- `scripts/theory/path_cover_demo.py`
- `scripts/theory/anchors_demo.py`
- `scripts/theory/precedence_demo.py`
- `scripts/theory/chaining_demo.py`
- `scripts/theory/visualization_demo.py`

Current theory capabilities:

- representative graph extraction
- reachability context generation
- synthetic cyclic demo preparation when real graphs are DAG-like
- Tarjan SCC detection
- DFS back-edge detection
- DAG approximation by back-edge removal
- simplified greedy path cover
- simplified anchor generation
- precedence relation checking
- simplified co-linear chaining dynamic programming
- gap-cost analysis
- iterative convergence demonstration
- theory visual generation and text reporting

Generated outputs include:

- `outputs/theory/scc_analysis.txt`
- `outputs/theory/scc_graph.png`
- `outputs/theory/back_edges.txt`
- `outputs/theory/dag_graph.png`
- `outputs/theory/path_cover.txt`
- `outputs/theory/path_cover_graph.png`
- `outputs/theory/anchors.csv`
- `outputs/theory/anchor_examples.txt`
- `outputs/theory/precedence_analysis.txt`
- `outputs/theory/chaining_results.txt`
- `outputs/theory/chaining_graph.png`
- `outputs/theory/gap_cost_analysis.txt`
- `outputs/theory/iteration_log.txt`
- `outputs/theory/convergence_plot.png`
- `outputs/theory/theory_context.json`

### 9. Held-out evaluation workflow

Implemented in:

- `scripts/paper_evaluation.py`

Current capabilities:

- aligns held-out test sequences to train-derived combined graphs
- records per-query alignment outputs
- computes overall metrics
- computes per-gene metrics
- computes per-bucket metrics for healthy and unhealthy groups
- creates MQ cutoff sweep analysis
- builds local graph property tables
- creates example alignment visualizations
- writes paper-style local summary files

Generated evaluation outputs include:

- `outputs/evaluation/alignment_metrics.json`
- `outputs/evaluation/alignment_records.json`
- `outputs/evaluation/evaluation_report.txt`
- `outputs/evaluation/paper_style_results.json`
- `outputs/evaluation/mq_cutoff_sweep.csv`
- `outputs/evaluation/mq_cutoff_sweep.png`
- `outputs/evaluation/app.mq_cutoff_sweep.png`
- `outputs/evaluation/psen1.mq_cutoff_sweep.png`
- `outputs/evaluation/psen2.mq_cutoff_sweep.png`
- `outputs/evaluation/local_graph_properties.csv`
- `outputs/evaluation/dataset_specific_reproduction_report.txt`
- `outputs/evaluation/*.example_alignment.png`

### 10. Local hybrid Minigraph-first evaluation logic

Implemented in:

- `scripts/paper_evaluation.py`

Current capabilities:

- runs a local hybrid approximation
- uses `minigraph` first for easy cases
- falls back to `PanAligner` when needed
- tracks which aligner produced the selected result

This is implemented as a local dataset-specific approximation, not as the paper's full benchmarked hybrid pipeline.

### 11. Custom query analysis

Implemented in:

- `scripts/predictor.py`
- `main.py --custom-query-analysis`

Current capabilities:

- accepts inline DNA sequence or external FASTA
- can infer the most likely gene by aligning against combined graphs
- aligns the query against healthy, unhealthy, and combined graphs for the chosen gene
- computes a composite score from identity, coverage, MAPQ, and normalized score
- predicts healthy vs unhealthy state
- produces a confidence value
- stores an explanatory summary
- generates a score plot

Generated custom-query outputs include:

- `outputs/alignments/custom_query/prediction.json`
- `outputs/alignments/custom_query/custom_query_scores.png`
- custom GAF files for the query

### 12. Desktop GUI

Implemented in:

- `gui_app.py`
- `run_gui.bat`

Current GUI capabilities:

- launches the current workflow from a desktop app
- runs `full-pipeline`, `theory-only`, and `evaluate`
- runs custom query analysis
- accepts inline DNA sequence or FASTA input
- supports binary-path configuration for PanAligner and minigraph
- shows live console output
- shows a run summary
- previews generated PNG graphs
- displays evaluation metrics in dedicated panels
- shows a metric cheatsheet
- shows local reproduction summaries
- shows paper coverage audit text
- browses generated files such as TXT, JSON, CSV, and GAF
- opens output files and folders directly

### 13. Reporting layer

Implemented in:

- `main.py`
- `gui_app.py`
- `scripts/paper_evaluation.py`

Current report outputs include:

- `outputs/reports/theoretical_reproduction_report.txt`
- `outputs/reports/final_project_summary.txt`
- `docs/paper_coverage_audit.md`
- `outputs/evaluation/dataset_specific_reproduction_report.txt`
- `outputs/evaluation/paftools_readiness_report.txt`

### 14. Tooling and platform helpers

Implemented support files include:

- `tools/minigraph_wsl.cmd`
- `tools/panaligner_wsl.cmd`
- `tools/paftools_wsl.cmd`
- `tools/wsl_exec.py`
- `scripts/setup_tools.sh`
- `run_pipeline.sh`

These make it easier to run the Linux-oriented toolchain from the current mixed Windows/WSL setup.

## Current output snapshot already present in the repo

From the current generated evaluation artifacts:

- total held-out queries evaluated: `46`
- aligned queries: `46`
- alignment rate: `1.0000`
- mean identity: `0.9999767438`
- mean coverage: `0.9999854009`
- mean MAPQ: `60.0`
- mean normalized score: `0.9999621425`

Per gene:

- `APP`: 10 queries, 10 aligned
- `PSEN1`: 32 queries, 32 aligned
- `PSEN2`: 4 queries, 4 aligned

Hybrid summary from current outputs:

- hybrid mode enabled: `true`
- `minigraph` selected: `45`
- `PanAligner` selected: `1`
- easy cases detected by `minigraph`: `45`

Current local graph property snapshot:

- `APP` combined train graph: 96 nodes, 130 edges, no cycle detected
- `PSEN1` combined train graph: 323 nodes, 428 edges, no cycle detected
- `PSEN2` combined train graph: 34 nodes, 43 edges, no cycle detected

## Main workflow in the current version

The active workflow is:

1. preprocess the root FASTA inputs
2. split into train and test sets
3. build train graphs with `minigraph`
4. visualize the graphs and compute graph statistics
5. run the theory suite on a representative graph
6. evaluate held-out sequences using real alignment
7. generate reports, JSON summaries, plots, and GUI-viewable artifacts

## What is implemented but not part of the current main story

The repository still contains earlier or archival modules from the project's older ML-oriented phase.

These include:

- `scripts/pipeline.py`
- `scripts/extract_features.py`
- `scripts/train_ml.py`
- `scripts/predict_ml.py`
- `scripts/evaluate.py`

What they do:

- extract alignment-derived features
- train ML models
- evaluate classification metrics
- support ML-based prediction on query sequences

Status in the current version:

- still implemented in code
- still part of the repository history
- not the main workflow described by the current README and `main.py`
- best treated as legacy or archival functionality unless the project scope is intentionally expanded back toward ML

## What is only partially reproduced

The following are implemented as simplified educational reproductions rather than exact low-level PanAligner internals:

- SCC-based theory preprocessing
- DAG approximation for cyclic graphs
- path cover logic
- anchor examples
- precedence relation
- co-linear chaining DP
- convergence demonstration

These are suitable for explanation, viva defense, and conceptual understanding, but they should not be described as a full reimplementation of PanAligner's optimized internal algorithm.

## What is not implemented yet

Compared with the full PanAligner paper benchmark suite, the current version does **not** implement:

- the exact large-scale 10H, 40H, 80H, and 95H benchmark reproduction
- full multi-aligner runtime and memory comparison tables
- the complete paftools-based truth-scoring workflow from simulated reads
- all paper figures and all reported benchmark tables
- full correctness benchmarking against the exact published external datasets

## Best one-line description for viva or report

**The current version implements a real PanAligner-based graph-alignment workflow, a theory/demo layer for SCC-to-chaining concepts, a held-out local evaluation pipeline with paper-style artifacts, a custom query analysis module, and a desktop GUI for running and inspecting all of it.**

## Best honest conclusion

As of the current version, the project is **substantially implemented** as a working local reproduction and explanation platform for the PanAligner paper.

It includes:

- real graph construction
- real alignment execution
- real held-out evaluation
- theory demonstrations
- output visualization
- reporting
- GUI-based usability support

It should be presented as a **practical partial reproduction plus educational implementation**, not as a complete reproduction of every benchmark and table from the original paper.
