# 02. Pangenome Graphs

## What it is

A pangenome graph is a graph-based representation of many related genome sequences together instead of storing only one linear reference.

### Simple definition

- nodes store sequence segments
- edges connect segments that can come next
- different paths represent different sequence variants

## Why it is needed

A single reference genome is limited. Real biological data contains variation:

- substitutions
- insertions
- deletions
- alternate sequence paths

If we align everything to only one linear reference, some variation is hard to represent cleanly. A graph can represent many valid alternatives.

## Simple intuition

Think of a pangenome graph like a city road map.

- a straight reference genome is one road
- a pangenome graph is the full road network
- different people may travel through slightly different streets
- the graph stores all these possible paths

## Technical explanation

In this project, graphs are stored in `GFA` format.

Each graph contains:

- `S` lines for sequence segments
- `L` lines for links between segments

Graph nodes correspond to sequence pieces. Graph edges represent adjacency or allowed traversal between those sequence pieces.

For each gene in this project, three graphs are built:

- healthy-only graph
- unhealthy-only graph
- combined graph

The combined graph is especially important because it holds variation from both groups.

## How our project builds them

The script `scripts/build_graph.py` uses `minigraph`.

### Core logic

1. read the reference FASTA
2. collect all sample FASTA files for one gene
3. run `minigraph`
4. write the graph as `.gfa`

### Command shape

```text
minigraph -t <threads> -cxggs -l1k -L1 reference.fa sample1.fa sample2.fa ...
```

You do not need to memorize every parameter, but you should know the role:

- `minigraph` constructs the graph
- reference is the starting backbone
- samples add variation branches

## How the repo uses pangenome graphs

For each gene:

- `APP`
- `PSEN1`
- `PSEN2`

the project builds:

- `graphs/train/<gene>.healthy.gfa`
- `graphs/train/<gene>.unhealthy.gfa`
- `graphs/train/<gene>.combined.gfa`

These train graphs are then used for held-out alignment evaluation.

## What the graph means biologically

The graph is not just a picture. It is a compact representation of sequence variation around the gene region.

Different paths through the graph can represent:

- reference-like sequence
- healthy sample variations
- unhealthy sample variations
- shared subsequences

## Example intuition

Suppose three sequences are:

```text
Reference: A C G T A
Sample 1 : A C G T A
Sample 2 : A C T T A
```

Instead of storing two full copies, a graph can store:

```text
A -> C -> G -> T -> A
         \
          -> T -> T -> A
```

This is simplified, but the idea is correct: shared parts are reused, and differences become alternate branches.

## Why this matters for alignment

When a new query comes in, it may match one branch better than another.

A graph aligner does not just ask:

`Where does this sequence fit on one line?`

It asks:

`Which path through this graph best explains the sequence?`

That is more flexible and better suited for variable genomic regions.

## Real graph stats from this repo

For `outputs/graphs/app.combined.stats.json`:

- nodes: `115`
- edges: `154`
- weakly connected components: `1`
- strongly connected components: `115`
- largest SCC size: `1`
- cycle detected: `false`

### What this means

This APP graph is connected and mostly DAG-like in the current dataset. That is why the theory layer uses a synthetic educational cyclic example when teaching cyclic handling.

## Beginner-friendly comparison

| Linear reference | Pangenome graph |
|---|---|
| One path | Many possible paths |
| Easier to visualize | Richer representation |
| Can miss variation patterns | Better captures variation |
| Good baseline | Better for diverse genomes |

## ASCII diagram

```text
Linear:
A -> C -> G -> T -> A

Graph:
        -> G ->
A -> C -        - > T -> A
        -> T ->
```

## How our project visualizes graphs

The script `scripts/visualize.py`:

- parses the GFA
- builds a directed graph with `networkx`
- computes graph statistics
- saves a `.png` plot
- saves a `.stats.json` file

## Common viva questions

### What is a pangenome graph?

A graph representation of multiple related genomes or sequences where nodes hold sequence segments and different paths represent different variants.

### Why use a graph instead of a linear reference?

Because graphs can represent multiple sequence alternatives and shared structure at the same time.

### What is the role of minigraph?

It constructs pangenome graphs from the reference and sample sequences.

### Why build healthy, unhealthy, and combined graphs?

To compare variation structure and to create a broader combined graph for alignment evaluation.

### What is stored in a graph node?

A sequence segment.

### What is stored in a graph edge?

A valid connection or transition between segments.

## Common mistakes and confusions

- Saying nodes store whole genomes
  Usually they store sequence segments, not entire genomes.

- Saying every graph must contain cycles
  No. Some graphs are acyclic or mostly DAG-like.

- Confusing graph construction with graph alignment
  Construction creates the graph; alignment maps a query onto it.

- Assuming more nodes always means better quality
  More nodes only means more structure, not automatically better biological utility.

## Quick viva answer

`A pangenome graph stores multiple related sequences in one graph structure, where nodes are sequence segments and alternate paths represent variation. In our project, minigraph builds healthy, unhealthy, and combined graphs for APP, PSEN1, and PSEN2, and PanAligner aligns held-out sequences onto those graphs.`
