# 12. Complete One-shot Guide

## Start Here

If you study only this file, you should be able to explain the project from start to finish in viva.

## 1. Project in One Sentence

This project reproduces the PanAligner paper by building pangenome graphs from FASTA data, aligning held-out sequences to those graphs using the real PanAligner, and separately explaining difficult paper concepts like SCCs and co-linear chaining using simplified educational demos.

## 2. The Most Important Distinction

### Real PanAligner implementation

- uses the official `PanAligner` binary
- performs actual sequence-to-graph alignment
- produces real `GAF` outputs
- is used in held-out evaluation

### Our simplified educational demos

- are written in Python
- live in `scripts/theory/`
- explain paper concepts
- do not replace PanAligner internals

### Best viva-safe sentence

`We used the real PanAligner for alignment, and we built separate simplified demos only to explain the paper concepts clearly.`

## 3. What the Project Does End-to-End

```text
Combined FASTA input
   ->
Preprocess and separate reference/healthy/unhealthy
   ->
Create train/test split
   ->
Build pangenome graphs from training sequences with minigraph
   ->
Visualize graph structure
   ->
Run educational theory demos
   ->
Align held-out test queries to train-derived combined graphs with PanAligner
   ->
Parse GAF alignments
   ->
Compute evaluation metrics
   ->
Write reports
```

## 4. Inputs and Data

The project uses three genes:

- `APP`
- `PSEN1`
- `PSEN2`

For each gene, the data contains:

- one reference sequence
- healthy sample sequences
- unhealthy sample sequences

### Current counts

- healthy sequences: `110`
- unhealthy sequences: `110`
- references: `3`

### Why reference is needed

The reference acts as the initial backbone for graph construction.

## 5. FASTA Preprocessing

### What it is

The first data-cleaning stage.

### Why it is needed

Raw combined FASTA files are not organized enough for graph building or held-out evaluation.

### What it does

- reads root FASTA files
- parses sequence metadata from headers
- normalizes sequences
- validates DNA
- separates:
  - reference
  - healthy
  - unhealthy
- writes one-sample FASTA files
- writes combined FASTA files
- writes metadata manifest

### Important file

`scripts/preprocess.py`

### Important output

`data/metadata/preprocess_manifest.json`

## 6. Train/Test Split

### Why split?

To evaluate graph alignment on unseen sequences.

### What it does

- healthy and unhealthy are split separately for each gene
- default seed: `42`
- default test fraction: `0.20`

### Current split summary

| Gene | Total | Train H | Test H | Train U | Test U |
|---|---:|---:|---:|---:|---:|
| APP | 46 | 18 | 5 | 18 | 5 |
| PSEN1 | 156 | 62 | 16 | 62 | 16 |
| PSEN2 | 18 | 7 | 2 | 7 | 2 |

### Why this is good

It gives a fairer evaluation than aligning only sequences already used in graph construction.

### Important file

`scripts/split_dataset.py`

## 7. Pangenome Graph Construction

### What is a pangenome graph?

A graph representation of multiple related sequences where:

- nodes store sequence segments
- edges connect possible next segments
- different paths represent different variants

### Simple analogy

A pangenome graph is like a road network, while a linear reference is like one single road.

### Why use it?

Because many samples share most sequence but differ at some points. A graph stores shared parts once and branches where variation happens.

### Which tool builds the graph?

`minigraph`

### What graphs are built per gene?

- healthy graph
- unhealthy graph
- combined graph

### Why combined graph matters

It captures the broadest variation and is used for held-out alignment evaluation.

### Important file

`scripts/build_graph.py`

## 8. Graph Representation Formats

### GFA

Graph format.

- `S` lines define segments
- `L` lines define links

### GAF

Alignment format.

Stores how a query aligns to a graph path.

### Memory trick

- `GFA` = `A` for assembly/graph structure
- `GAF` = `A` for alignment format in practice here

Better simple memory:

- `GFA` is graph
- `GAF` is alignment

## 9. Graph Visualization

### Why visualize?

To understand graph topology and show it clearly in viva.

### What the project computes

- node count
- edge count
- connected components
- SCC count
- cycle detection

### Example from current APP graph

- nodes: `115`
- edges: `154`
- cycle detected: `false`

### Meaning

The current APP graph is mostly DAG-like, not strongly cyclic.

## 10. Why Cyclic Graphs Matter Anyway

The PanAligner paper is specifically interested in co-linear chaining on pangenome graphs, including cyclic cases.

Even if our current real graph is often DAG-like, the theory still matters because:

- the paper studies general cyclic handling
- viva questions may focus on the paper, not only the current dataset

So the project includes a teaching layer for cyclic graph concepts.

## 11. SCC: Strongly Connected Components

### What is an SCC?

A set of nodes where every node can reach every other node.

### Why important?

It identifies cyclic regions.

### Simple analogy

An SCC is like a neighborhood where every road can eventually lead to every other road.

### Algorithm used

Tarjan's algorithm.

### How to explain Tarjan simply

It performs DFS, keeps a stack, and uses lowlink values to detect when a strongly connected group is complete.

### Important file

`scripts/theory/scc_demo.py`

## 12. Educational Cyclic Graph Construction

### Honest explanation

The real locus graphs in this repo are often mostly acyclic. So the theory layer:

1. picks a representative subgraph from a real graph
2. optionally adds synthetic back edges
3. uses that as the educational cyclic graph

### Why this is acceptable

Because it does not alter the real alignment workflow. It only creates a small teaching graph for explaining the paper's cyclic case.

### Very useful viva sentence

`The synthetic back edges are added only in the educational copy, not in the real production graph used for alignment.`

## 13. DAG Conversion

### What is a DAG?

A directed acyclic graph.

### Why convert cyclic graph to DAG approximation?

Because ordered processing becomes easier when cycles are removed.

### How our demo does it

- run DFS
- detect back edges
- remove those back edges in the educational graph copy

### Important file

`scripts/theory/dag_converter_demo.py`

### Important honesty note

This is a simplified educational approach, not a claim about exact PanAligner internals.

## 14. Path Cover

### What is path cover?

A set of paths that together covers all nodes in a DAG.

### Why useful?

It gives structure to the DAG and helps explain ordered traversal ideas.

### How our project does it

The demo uses a greedy path cover algorithm.

### Important file

`scripts/theory/path_cover_demo.py`

### Current demo report

The current theory report mentions `3` paths in the simplified path cover.

## 15. Anchors

### What is an anchor?

A local match between the query and a graph location.

### Representation in our demo

```text
(vertex, [x..y], [c..d], weight)
```

Where:

- `vertex` = graph node
- `[x..y]` = graph interval
- `[c..d]` = query interval
- `weight` = anchor strength

### Why anchors matter

Chaining works on anchors, not on the whole query at once.

### Important file

`scripts/theory/anchors_demo.py`

## 16. Precedence Relation

### What is it?

A rule deciding whether anchor `A` can validly come before anchor `B`.

### Simplified rule in this project

Anchor `A` precedes anchor `B` if:

- query order increases
- graph reachability exists

### Why it matters

Without precedence, chaining could connect impossible anchor pairs.

### Important file

`scripts/theory/precedence_demo.py`

## 17. Co-linear Chaining

### Core idea

Choose the best ordered sequence of compatible anchors.

### Why it is needed

A good alignment is usually formed by combining multiple local matches, not by looking at one anchor alone.

### Simple analogy

Like picking the best sequence of matching puzzle pieces that fit in the correct order.

### DP recurrence

```text
dp[j] = best score of a chain ending at anchor j
```

For each earlier compatible anchor `i`:

```text
candidate = dp[i] + weight(j) - gap_cost(i, j)
```

### Gap cost

The project uses:

- query gap
- graph gap
- total gap

### Important file

`scripts/theory/chaining_demo.py`

### Current output

- best chain: `A6 -> A7 -> A8`
- best score: `12.00`

## 18. Real PanAligner Workflow

### What PanAligner does

It aligns a query sequence to a graph.

### Input

- graph `GFA`
- query `FASTA`

### Output

- alignment `GAF`

### Command shape

```text
PanAligner -t 4 -cx lr graph.gfa query.fa > out.gaf
```

### Important file

`scripts/align.py`

## 19. GAF Parsing

### What the parser extracts

- query name
- query span
- path span
- identity
- coverage
- mapping quality
- alignment score
- traversed graph nodes

### Why parse GAF?

So that raw alignments become measurable and explainable.

### Important file

`scripts/parse_gaf.py`

## 20. Held-out Evaluation Pipeline

### Why held-out evaluation matters

It shows whether the graph can explain unseen sequences.

### How the project evaluates

1. build graphs using train data
2. use test data as held-out queries
3. align each query to the corresponding combined train graph
4. choose best GAF alignment
5. compute summary metrics

### Important file

`scripts/paper_evaluation.py`

## 21. Current Evaluation Results

### Overall

- query count: `46`
- aligned query count: `46`
- alignment rate: `1.0`
- mean identity: `0.9995539276605313`
- mean coverage: `1.0000659121670492`
- mean MAPQ: `60.0`
- mean alignment score: `128845.0`

### Per gene

| Gene | Queries | Aligned | Mean identity |
|---|---:|---:|---:|
| APP | 10 | 10 | 0.999756 |
| PSEN1 | 32 | 32 | 0.999437 |
| PSEN2 | 4 | 4 | 0.999983 |

### Easy interpretation

The held-out sequences aligned extremely well to the combined train graphs.

## 22. What `main.py` Does

### `--full-pipeline`

Runs everything end to end.

### `--theory-only`

Runs only the educational explanation modules.

### `--evaluate`

Runs held-out PanAligner evaluation.

### Why this is good

It gives one clean entry point for the whole project.

## 23. Practical Outputs You Should Mention

### Graph outputs

- graph PNGs
- graph stats JSON

### Theory outputs

- SCC analysis
- back edge report
- path cover report
- anchor examples
- precedence analysis
- chaining results

### Evaluation outputs

- GAF alignments
- alignment metrics JSON
- alignment records JSON
- evaluation report

## 24. Common Viva Questions and Best Answers

### What is the project doing overall?

It constructs pangenome graphs and evaluates sequence-to-graph alignment while explaining the paper's core algorithms in a simplified educational form.

### Why use pangenome graphs?

Because they represent multiple sequence variants in one structure rather than forcing everything onto a single linear reference.

### Why split into train and test?

To evaluate on unseen queries instead of sequences already used to build the graph.

### What is the role of SCC?

To identify cyclic regions in directed graphs.

### Why convert to DAG in the demo?

To simplify ordering and graph reasoning for downstream explanation like path cover and chaining.

### What is co-linear chaining?

Selecting the best ordered compatible anchors while balancing weights and gap penalties.

### Did you implement the exact PanAligner internals?

No. We used the real PanAligner binary for practical alignment and separate simplified demos for explanation.

## 25. Common Confusions to Avoid

### Confusion 1

`PanAligner builds graphs`

Wrong. `minigraph` builds graphs.

### Confusion 2

`The theory demo is the actual PanAligner algorithm`

Wrong. It is an educational approximation.

### Confusion 3

`Current real graphs are strongly cyclic everywhere`

Wrong. Many are mostly DAG-like in this dataset.

### Confusion 4

`This project proves disease prediction`

Wrong. The current main focus is graph alignment and paper reproduction.

## 26. Best Diagrams to Remember

### Full workflow

```text
FASTA
  ->
Preprocess
  ->
Split
  ->
Train Graph Build
  ->
Theory Demos
  ->
PanAligner Alignment
  ->
GAF Parsing
  ->
Evaluation
```

### Graph concept

```text
Linear:
A -> B -> C -> D

Graph:
      -> C1 ->
A -> B         -> D
      -> C2 ->
```

### Cycle concept

```text
A -> B -> C
     ^    |
     |    v
     E <- D
```

### Chaining concept

```text
Anchors in query order:
A1 -> A2 -> A3

Graph-compatible route:
v1 -> v2 -> v3

Then:
A1 -> A2 -> A3 is a valid chain
```

## 27. 30-second Viva Summary

`Our project reproduces the PanAligner paper by taking FASTA data for APP, PSEN1, and PSEN2, preprocessing and splitting it into train and test sets, building pangenome graphs with minigraph, and aligning held-out test sequences to train-derived combined graphs using the real PanAligner. To explain the paper clearly, we also built a separate theory layer that demonstrates SCCs, DAG conversion, path cover, anchor precedence, and simplified co-linear chaining in an easy-to-understand way.`

## 28. Final Quick Memory Sheet

- `minigraph` builds graphs
- `PanAligner` aligns queries
- `GFA` is graph format
- `GAF` is alignment format
- `SCC` finds cyclic regions
- `DAG` has no cycles
- `path cover` covers DAG nodes with paths
- `anchor` is a local match
- `precedence` says which anchor can come first
- `chaining` picks the best valid anchor sequence
- real aligner and theory demo are separate

## 29. Final One-line Answer

`This project is a practical and educational reproduction of the PanAligner paper: it runs real graph alignment with PanAligner and separately teaches the cyclic graph chaining concepts needed to explain the paper confidently in viva.`
