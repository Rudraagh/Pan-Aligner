# 06. PanAligner Workflow

## What it is

This file explains how the real `PanAligner` part of the project works from input graph to output alignment.

## One-sentence summary

`PanAligner takes a graph and a query sequence, aligns the query to the graph, and produces GAF output describing the alignment path and quality.`

## Why it is needed

Building a graph is not enough. We also need to test whether unseen sequences can align well to that graph.

That is the main practical purpose of PanAligner in this project.

## Simple intuition

If `minigraph` builds the road network, `PanAligner` drives a new vehicle through that network and reports the best route.

## Inputs to PanAligner

The real aligner uses:

- one graph in `GFA`
- one query sequence in `FASTA`
- optional runtime parameters like threads

## How our project wraps PanAligner

The script `scripts/align.py` contains:

- `run_panaligner(...)`

It executes a command of the form:

```text
PanAligner -t 4 -cx lr graph.gfa query.fa > output.gaf
```

### Meaning of main pieces

- `-t 4`: use 4 threads
- `-cx lr`: use long-read style preset
- `graph.gfa`: graph to align against
- `query.fa`: query sequence
- `output.gaf`: saved alignment output

## End-to-end practical flow

```text
Train FASTA
   ->
minigraph builds train graph
   ->
held-out test FASTA chosen as query
   ->
PanAligner aligns query to combined graph
   ->
GAF file produced
   ->
parser extracts best alignment
   ->
evaluation metrics recorded
```

## How held-out evaluation uses PanAligner

The script `scripts/paper_evaluation.py`:

1. reads test manifest
2. reads train graph manifest
3. for each test sequence:
   aligns it to the corresponding train-derived combined graph
4. parses the best GAF alignment
5. stores per-sequence record
6. computes overall and per-gene metrics

## Why train-derived graphs are used

This is important for viva.

We do not want the query sequence to already be part of the graph during evaluation, because then evaluation would be less meaningful.

So the project:

- builds graphs from train sequences
- aligns held-out test sequences later

That creates a cleaner reproduction-style evaluation.

## Real repo evaluation behavior

For each gene:

- query comes from `data/test/...`
- graph comes from `graphs/train/...combined.gfa`

That is good scientific workflow because it separates building and testing.

## Current real results

From the project outputs:

- total held-out queries: `46`
- aligned queries: `46`
- alignment rate: `100%`
- mean identity: about `0.99955`
- mean coverage: about `1.00007`
- mean MAPQ: `60`

### Interpretation

The queries align extremely well to the train-derived combined graphs. This suggests the graph captures the relevant sequence structure effectively for these loci.

## Per-gene results

### APP

- query count: `10`
- aligned: `10`
- mean identity: about `0.99976`

### PSEN1

- query count: `32`
- aligned: `32`
- mean identity: about `0.99944`

### PSEN2

- query count: `4`
- aligned: `4`
- mean identity: about `0.99998`

## Visualization support

The project also creates alignment visuals by:

- parsing traversed nodes from GAF
- highlighting those nodes on the graph

This is done by `visualize_alignment(...)` in `scripts/visualize.py`.

## What PanAligner is not doing here

PanAligner is not:

- building the graph
- splitting the dataset
- teaching SCC theory

Those jobs are handled by other parts of the repo.

## Real PanAligner vs theory layer

| Task | Real PanAligner | Theory scripts |
|---|---|---|
| sequence-to-graph alignment | Yes | No |
| SCC demo | No | Yes |
| DAG conversion demo | No | Yes |
| simplified chaining explanation | No | Yes |
| output GAF | Yes | No |

## Common viva explanation

If asked, say:

`PanAligner is the real alignment engine in our project. The theory modules explain the paper concepts, but they do not perform the actual graph alignment.`

## Common viva questions

### What is the role of PanAligner in your project?

It performs real query-to-graph alignment on the graphs built by minigraph.

### What does PanAligner take as input?

A graph in GFA format and a query sequence in FASTA format.

### What does it output?

A GAF alignment file.

### Why use held-out queries?

To evaluate alignment on sequences that were not used to build the train graph.

### How do you choose the best alignment?

The project parses all GAF records and chooses the best one using normalized score, identity, coverage, and mapping quality.

### Did you modify PanAligner internals?

No. We integrate and use the official PanAligner binary.

## Common mistakes and confusions

- Saying PanAligner builds pangenome graphs
  That is `minigraph`'s role.

- Saying PanAligner output is JSON
  It outputs `GAF`; the project later converts summaries into JSON.

- Forgetting the train/test separation
  That is a key strength of the evaluation design.

- Confusing theory demo output with PanAligner output
  PanAligner produces GAF; theory demos produce educational plots and text reports.

## Fast viva answer

`In our project, PanAligner is the real sequence-to-graph aligner. It takes a train-derived combined GFA graph and a held-out query FASTA, produces GAF alignments, and those alignments are parsed to compute identity, coverage, MAPQ, and alignment rate.`
