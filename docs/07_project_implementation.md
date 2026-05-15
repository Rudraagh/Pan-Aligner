# 07. Project Implementation

## What it is

This file explains how the repository is implemented in practice.

## Big picture

The implementation has two layers:

1. `Practical pipeline`
2. `Educational theory pipeline`

Both are connected, but their purposes are different.

## Layer 1: Practical pipeline

### `scripts/preprocess.py`

#### What it does

- reads input FASTA files
- parses headers
- normalizes sequences
- validates DNA
- separates healthy, unhealthy, and reference sequences
- writes organized FASTA files and manifest JSON

#### Why it is needed

Raw combined FASTA files are not convenient for graph building and evaluation. Preprocessing organizes the data into a clean project structure.

#### Important implementation details

- deduplication is supported
- invalid sequences are skipped
- each gene must have exactly one reference sequence

#### Output

- `data/healthy/`
- `data/unhealthy/`
- `data/raw/`
- `data/metadata/preprocess_manifest.json`

### `scripts/split_dataset.py`

#### What it does

- creates reproducible train/test splits
- splits healthy and unhealthy separately
- writes both per-sample and combined FASTA outputs

#### Why it is needed

To evaluate alignment on held-out sequences rather than sequences already used to build the graph.

#### Important implementation details

- random seed based
- class-aware within each gene
- always keeps at least one train and one test sequence per class when possible

#### Output

- `data/train/...`
- `data/test/...`
- `data/metadata/train_manifest.json`
- `data/metadata/test_manifest.json`
- `data/metadata/split_manifest.json`
- `data/metadata/split_summary.txt`

### `scripts/build_graph.py`

#### What it does

- reads manifest data
- resolves reference and sample FASTA paths
- calls `minigraph`
- builds three graphs per gene

#### Graphs built

- healthy graph
- unhealthy graph
- combined graph

#### Output

- `graphs/train/*.gfa`
- `data/metadata/train_graph_manifest.json`

### `scripts/align.py`

#### What it does

- runs the real `PanAligner` binary
- saves GAF output
- parses GAF immediately for convenience

### `scripts/parse_gaf.py`

#### What it does

- parses GAF into structured Python objects
- computes identity, coverage, normalized score, and traversed nodes

### `scripts/visualize.py`

#### What it does

- parses GFA into a directed graph
- computes graph stats
- plots graphs
- optionally highlights nodes touched by an alignment

### `scripts/paper_evaluation.py`

#### What it does

- loops over test sequences
- aligns each query to the appropriate train-derived combined graph
- selects best alignment
- computes per-gene and overall evaluation metrics

## Layer 2: Educational theory pipeline

This layer lives inside `scripts/theory/`.

### `common.py`

Builds the teaching context:

- selects a source graph
- extracts a representative subgraph
- may add synthetic back edges
- creates a demo query walk
- computes reachability

### `graph_analysis_demo.py`

Creates a summary of:

- original graph
- representative graph
- educational cyclic graph

### `scc_demo.py`

Implements Tarjan SCC detection.

### `dag_converter_demo.py`

Finds DFS back edges and removes them to create a DAG approximation.

### `path_cover_demo.py`

Creates a simplified greedy path cover on the DAG.

### `anchors_demo.py`

Creates demo anchors for the query walk plus a few distractor anchors.

### `precedence_demo.py`

Decides which anchors can validly come before others.

### `chaining_demo.py`

Runs simplified co-linear chaining DP and iterative convergence demo.

## Implementation design strengths

### Clear separation of concerns

- preprocessing is isolated
- graph building is isolated
- alignment is isolated
- parsing is isolated
- theory explanation is isolated

This makes the project easier to explain and debug.

### Manifest-driven workflow

The project uses JSON manifests to track:

- what data exists
- where files are stored
- which graphs belong to which genes

This is better than hardcoding file paths everywhere.

### Reproducibility

The split uses a fixed seed by default:

- seed: `42`
- test fraction: `0.20`

That helps repeat the same experiment.

## Important implementation distinction

### Real implementation

- uses external binaries
- produces actual graph/alignment outputs
- evaluated on held-out sequences

### Educational demos

- use `networkx`
- use smaller representative graphs
- use simplified scoring and chaining
- exist mainly for explanation

## How our project handles cyclic teaching honestly

This is a likely viva question.

The current real graphs are often mostly DAG-like. To still demonstrate the cyclic case from the paper:

1. the project extracts a representative subgraph
2. synthetic back edges may be added
3. SCC and DAG conversion are demonstrated on that educational copy

That means the theoretical concept is shown without altering the real alignment graphs.

## Implementation flow in one table

| Stage | Script | Main output |
|---|---|---|
| Preprocess | `preprocess.py` | clean FASTA + manifest |
| Split | `split_dataset.py` | train/test FASTA + manifests |
| Build graph | `build_graph.py` | train `.gfa` graphs |
| Visualize | `visualize.py` | `.png` and stats JSON |
| Theory | `scripts/theory/*` | educational reports and plots |
| Align | `align.py` | `.gaf` alignments |
| Evaluate | `paper_evaluation.py` | alignment metrics and report |

## Important outputs in the repo

### Graph outputs

- `outputs/graphs/*.png`
- `outputs/graphs/*.stats.json`

### Alignment outputs

- `outputs/alignments/*.gaf`

### Evaluation outputs

- `outputs/evaluation/alignment_metrics.json`
- `outputs/evaluation/alignment_records.json`
- `outputs/evaluation/evaluation_report.txt`

### Theory outputs

- `outputs/theory/scc_analysis.txt`
- `outputs/theory/back_edges.txt`
- `outputs/theory/path_cover.txt`
- `outputs/theory/anchor_examples.txt`
- `outputs/theory/precedence_analysis.txt`
- `outputs/theory/chaining_results.txt`

## Common viva questions

### Why did you split the project into practical and theory layers?

So the project can both run real alignment experiments and clearly teach the paper's harder concepts.

### Why use manifests?

They make the workflow organized, reproducible, and easy to connect across scripts.

### Why use train/test split in a graph alignment project?

To test graph alignment on unseen sequences, which is more meaningful than aligning only sequences used to construct the graph.

### Why is the theory layer written in Python?

Because Python with `networkx` is easier to read, visualize, and defend in viva than deeply optimized aligner internals.

## Common mistakes and confusions

- Assuming every script is equally important
  `main.py` is the main orchestrator; other scripts are components.

- Forgetting that some older ML files exist in the repo
  They are archival and not the main paper-reproduction path now.

- Saying the theory modules modify PanAligner
  They do not.

## Fast viva answer

`The implementation is modular: preprocessing organizes FASTA data, splitting creates held-out queries, minigraph builds train graphs, PanAligner performs real alignment, GAF parsing summarizes results, and a separate theory layer reproduces SCC, DAG, path cover, precedence, and chaining concepts in a small educational form.`
