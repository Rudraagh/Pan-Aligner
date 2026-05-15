# 01. Project Overview

## One-line idea

This project reproduces the PanAligner paper in two ways:

1. `Real implementation layer`: it actually runs `minigraph` and `PanAligner`.
2. `Educational theory layer`: it explains the paper's graph ideas with simple Python demos.

That is the most important sentence for your viva.

## What the project does

The project takes genomic sequence data for three genes:

- `APP`
- `PSEN1`
- `PSEN2`

and builds pangenome graphs from them. Then it aligns held-out test sequences to those graphs using the real `PanAligner` tool. Alongside that, it demonstrates the paper concepts such as:

- cyclic graph handling
- SCC detection
- DAG conversion
- path cover
- precedence relation
- simplified co-linear chaining DP

## Why this project matters

A normal linear reference genome is like using only one road map for all people.

But real biological variation means different people may have slightly different sequence paths. A pangenome graph stores many valid paths together. That makes graph alignment more flexible and biologically meaningful than plain reference alignment.

PanAligner is a tool that aligns sequences to pangenome graphs. The paper studies how to do co-linear chaining even when the graph has cycles.

## Project goal in simple words

The goal is not to invent a new aligner.

The goal is to:

- understand the PanAligner paper
- reproduce its workflow in a project
- show real graph construction and real alignment
- explain difficult theory in a beginner-friendly way
- defend the project confidently in viva

## Real PanAligner vs our educational demos

| Part | What it is | Why it exists |
|---|---|---|
| Real PanAligner implementation | Actual binary from the official repository | Shows true graph alignment on real data |
| Simplified educational demos | Python scripts in `scripts/theory/` | Helps explain the paper concepts clearly in viva |

### Very important viva sentence

`We did not replace PanAligner. We used the real PanAligner for alignment, and we built separate educational demos to understand and explain the paper concepts.`

## Main workflow

```text
Input FASTA
   ->
Preprocessing
   ->
Train/Test Split
   ->
Graph Construction with minigraph
   ->
Graph Visualization and Analysis
   ->
Theory Demos (SCC, DAG, chaining, etc.)
   ->
Real PanAligner Alignment
   ->
GAF Parsing
   ->
Held-out Evaluation
   ->
Reports
```

## Repository components

### Core files

- `main.py`: master entry point
- `README.md`: overall explanation
- `run_pipeline.sh`: convenience runner

### Practical pipeline scripts

- `scripts/preprocess.py`
- `scripts/split_dataset.py`
- `scripts/build_graph.py`
- `scripts/align.py`
- `scripts/parse_gaf.py`
- `scripts/visualize.py`
- `scripts/paper_evaluation.py`

### Theory scripts

- `scripts/theory/graph_analysis_demo.py`
- `scripts/theory/scc_demo.py`
- `scripts/theory/dag_converter_demo.py`
- `scripts/theory/path_cover_demo.py`
- `scripts/theory/anchors_demo.py`
- `scripts/theory/precedence_demo.py`
- `scripts/theory/chaining_demo.py`

## What the current implementation includes

The current repo already includes:

- FASTA preprocessing
- validation and deduplication
- train/test split
- graph construction using `minigraph`
- graph visualization
- real `PanAligner` execution
- GAF parsing
- held-out evaluation
- SCC and cyclic graph demos
- DAG conversion demo
- path cover demo
- precedence demo
- simplified chaining DP demo

## Current dataset summary

From the repository metadata:

| Gene | Total sequences | Train healthy | Test healthy | Train unhealthy | Test unhealthy |
|---|---:|---:|---:|---:|---:|
| APP | 46 | 18 | 5 | 18 | 5 |
| PSEN1 | 156 | 62 | 16 | 62 | 16 |
| PSEN2 | 18 | 7 | 2 | 7 | 2 |
| Overall | 220 | 87 | 23 | 87 | 23 |

There are also `3` reference sequences, one per gene.

## Current evaluation summary

From `outputs/evaluation/alignment_metrics.json`:

- query count: `46`
- aligned query count: `46`
- alignment rate: `1.0`
- mean identity: approximately `0.99955`
- mean coverage: approximately `1.00007`
- mean MAPQ: `60`

### Simple interpretation

This means the held-out test sequences aligned successfully to the train-derived combined graphs, and the alignments were extremely strong.

## Beginner intuition

Think of the whole project like this:

- FASTA files are raw travel routes.
- `minigraph` builds a road network from many routes.
- `PanAligner` checks where a new route fits in that network.
- the theory scripts explain how the aligner can reason about path order even if the network contains loops.

## Technical explanation

The practical layer:

- preprocesses genomic sequence data
- separates healthy, unhealthy, and reference sequences
- performs a reproducible train/test split
- builds pangenome graphs using `minigraph`
- aligns held-out sequences using the real `PanAligner`
- parses GAF output and summarizes alignment quality

The theory layer:

- extracts a representative graph
- adds synthetic back edges for teaching cyclic behavior when needed
- runs SCC analysis
- removes DFS back edges to create a DAG approximation
- computes a simplified path cover
- builds anchors
- computes precedence relations
- runs a simplified co-linear chaining dynamic program

## How the paper relates to the implementation

The paper is about co-linear chaining on pangenome graphs, especially in cyclic settings.

Your project reflects that in two forms:

1. `Practical reproduction`
   Real graphs are built and real alignments are performed.
2. `Conceptual reproduction`
   Hard algorithmic ideas are re-created in a small, teachable form.

## Common viva questions

### What is this project in one sentence?

It is a reproduction of the PanAligner paper that combines real graph alignment with beginner-friendly educational demos of the paper's main algorithms.

### Why did you include theory demos if PanAligner already exists?

Because PanAligner's internal implementation is optimized and harder to explain directly. The demos help us understand and defend the paper concepts clearly.

### Is this a machine learning project?

No. The current main project framing is pangenome graph construction, graph alignment, and paper reproduction.

### Did you reimplement the full PanAligner algorithm?

No. We used the real PanAligner for alignment and created simplified educational reproductions for explanation.

## Common mistakes and confusions

- Confusing `PanAligner` with `minigraph`
  `minigraph` builds graphs, `PanAligner` aligns sequences to graphs.

- Saying the theory demo is the full aligner
  It is not. It is only an educational approximation.

- Saying cycles were found directly in all real graphs
  In this repo, many real locus graphs are mostly DAG-like, so synthetic back edges are added only in the educational copy to teach cyclic handling.

- Confusing alignment evaluation with disease prediction
  The current pipeline evaluates sequence-to-graph alignment quality, not clinical diagnosis.

## Fast viva answer

`Our project reproduces the PanAligner paper by combining real pangenome graph alignment with simplified theory demos for SCCs, DAG conversion, path cover, precedence, and co-linear chaining, so we can both implement the workflow and explain it confidently in viva.`
