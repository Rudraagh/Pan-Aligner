# 11. Complete Revision Notes

## 1-minute revision

- Project uses `minigraph` to build pangenome graphs and `PanAligner` to align held-out sequences.
- Main genes: `APP`, `PSEN1`, `PSEN2`.
- Main entry point: `main.py`.
- Main modes: `--full-pipeline`, `--theory-only`, `--evaluate`.
- Data pipeline: preprocess -> split -> build graphs -> align -> parse GAF -> evaluate.
- Theory pipeline: graph analysis -> SCC -> DAG conversion -> path cover -> anchors -> precedence -> chaining.
- Biggest distinction: real PanAligner alignment versus simplified educational demos.

## 3-minute revision

### Project aim

Reproduce the PanAligner paper practically and conceptually.

### Practical tools

- `minigraph` builds graphs
- `PanAligner` aligns queries to graphs

### Formats

- `GFA` stores graph
- `GAF` stores graph alignment

### Why train/test split?

To build graphs on train data and evaluate on unseen test queries.

### Current split

- total sample sequences: `220`
- train sample sequences: `174`
- test sample sequences: `46`
- references: `3`

### Current evaluation

- aligned queries: `46/46`
- alignment rate: `1.0`
- mean identity: about `0.99955`
- mean coverage: about `1.00007`
- mean MAPQ: `60`

## 5-minute revision

### Pangenome graph

A graph that stores multiple related sequences together using nodes and edges.

### Why graph instead of linear reference?

Because graphs capture variation using alternate branches.

### Cyclic graphs

Graphs may contain loops, which make ordering harder.

### SCC

A strongly connected component is a set of nodes that can all reach each other.

### Why SCC matters

It identifies cyclic regions.

### DAG conversion

In the educational demo, DFS back edges are removed to produce a DAG approximation for simpler downstream explanation.

### Path cover

A set of paths covering all nodes in the DAG.

### Anchor

A local match between query and graph location.

### Precedence

One anchor can come before another if query order and graph reachability both make sense.

### Chaining

Choose the best ordered sequence of compatible anchors using dynamic programming and gap costs.

## Absolute must-remember distinctions

### Distinction 1

`minigraph` builds graphs.

### Distinction 2

`PanAligner` aligns queries to graphs.

### Distinction 3

`GFA` is graph format.

### Distinction 4

`GAF` is alignment format.

### Distinction 5

Real PanAligner alignment is not the same as the simplified theory demo.

## Best viva-safe sentences

### Sentence 1

`We used the real PanAligner binary for sequence-to-graph alignment.`

### Sentence 2

`We created separate educational Python demos to explain the paper's concepts like SCCs, DAG conversion, path cover, precedence, and co-linear chaining.`

### Sentence 3

`The theory layer is pedagogical and simplified; it does not replace PanAligner internals.`

### Sentence 4

`The held-out evaluation aligns test sequences to train-derived combined graphs, which gives a meaningful reproduction-style test.`

## Common viva traps

### Trap 1

Saying the theory DP is the actual PanAligner algorithm.

Correct response:

`No, it is a simplified educational version.`

### Trap 2

Saying PanAligner builds graphs.

Correct response:

`Graph construction is done by minigraph.`

### Trap 3

Saying current real graphs are strongly cyclic everywhere.

Correct response:

`Many current real graphs are mostly DAG-like, so the cyclic case is demonstrated on an educational copy with synthetic back edges.`

### Trap 4

Confusing evaluation with disease prediction.

Correct response:

`The current project focus is graph alignment behavior and paper reproduction, not diagnosis.`

## Quick command revision

```text
python main.py --full-pipeline
python main.py --theory-only
python main.py --evaluate
```

## Quick file revision

```text
main.py
scripts/preprocess.py
scripts/split_dataset.py
scripts/build_graph.py
scripts/align.py
scripts/parse_gaf.py
scripts/paper_evaluation.py
scripts/theory/scc_demo.py
scripts/theory/dag_converter_demo.py
scripts/theory/path_cover_demo.py
scripts/theory/precedence_demo.py
scripts/theory/chaining_demo.py
```

## Quick numbers revision

- genes: `3`
- references: `3`
- healthy sequences: `110`
- unhealthy sequences: `110`
- test fraction: `0.20`
- split seed: `42`
- held-out test queries evaluated: `46`
- aligned queries: `46`
- best current theory chain: `A6 -> A7 -> A8`

## Quick oral answers

### What is a pangenome graph?

A graph representation of multiple related sequences with alternate paths for variation.

### Why use PanAligner?

To align query sequences to pangenome graphs.

### Why use SCCs?

To identify cyclic regions.

### Why convert to DAG in the demo?

To simplify ordering and teaching of downstream graph algorithms.

### Why chaining?

To combine local matches into a strong global alignment explanation.

## 10-second final summary

`This project builds pangenome graphs with minigraph, aligns held-out sequences using the real PanAligner, and explains cyclic graph chaining concepts with separate educational demos for viva clarity.`
