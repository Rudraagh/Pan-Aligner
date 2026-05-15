# 04. GFA and GAF

## What they are

### GFA

`GFA` stands for `Graphical Fragment Assembly`.

It is a text format used to represent sequence graphs.

### GAF

`GAF` stands for `Graph Alignment Format`.

It is a text format used to store alignments of a query sequence against a graph.

## Why they are needed

The project needs:

- one format to store the graph
- another format to store alignments onto that graph

That is exactly what `GFA` and `GAF` provide.

## Simple intuition

- `GFA` is the map
- `GAF` is the travel report showing how a query moved through the map

## GFA in simple words

A GFA file mostly contains lines like:

- `S`: segment or node
- `L`: link or edge

### Very simple example

```text
S   s1   ACGT
S   s2   TTAA
L   s1   +   s2   +   0M
```

Meaning:

- node `s1` has sequence `ACGT`
- node `s2` has sequence `TTAA`
- there is a directed link from `s1` to `s2`

## How our project reads GFA

The script `scripts/visualize.py`:

- reads `S` lines
- creates graph nodes
- stores sequence length
- reads `L` lines
- creates directed edges

The theory layer also parses GFA and keeps actual segment sequences for educational examples.

## GAF in simple words

A GAF line tells how a query sequence aligned to a graph path.

It includes fields such as:

- query name
- query length
- query start and end
- strand
- graph path
- path length
- path start and end
- residue matches
- alignment block length
- mapping quality
- optional tags

## How our project parses GAF

The script `scripts/parse_gaf.py` defines an `AlignmentResult` object and computes:

- `coverage`
- `identity`
- `path_span`
- `alignment_score`
- `normalized_score`
- `traversed_nodes`

## Important derived measures

### Coverage

```text
alignment_block_length / query_length
```

This tells how much of the query was aligned.

### Identity

```text
residue_matches / alignment_block_length
```

This tells how similar the aligned part is.

### Normalized score

```text
alignment_score / query_length
```

This makes scores easier to compare across queries of different lengths.

## How best alignment is chosen

The function `best_alignment()` selects the best GAF record using:

1. normalized score
2. identity
3. coverage
4. mapping quality

That means the project does not blindly pick the first alignment.

## Path field intuition

In GAF, the path can look like a walk through graph nodes.

The project extracts traversed nodes by splitting the path on direction markers like `<` and `>`.

That helps visualize which graph nodes were involved in the alignment.

## How our project uses GFA and GAF together

### Step 1

`minigraph` produces `.gfa` files.

### Step 2

`PanAligner` takes:

- one `GFA` graph
- one query FASTA

and outputs a `GAF` file.

### Step 3

The project parses that GAF to produce alignment summaries and evaluation metrics.

## Real example from repo

Output directories include:

- `graphs/train/*.gfa`
- `outputs/alignments/*.gaf`
- `outputs/evaluation/alignment_records.json`

So the actual pipeline is:

```text
FASTA -> GFA -> GAF -> JSON metrics
```

## ASCII picture

```text
Reference + samples
        |
        v
   minigraph
        |
        v
      GFA
        |
        v
   PanAligner + query
        |
        v
      GAF
        |
        v
   parser + evaluation
```

## Common viva questions

### What is GFA?

A graph representation format that stores sequence segments and links between them.

### What is GAF?

A graph alignment format that stores how a query sequence aligns to a graph path.

### Which tool creates GFA in your project?

`minigraph`

### Which tool creates GAF in your project?

`PanAligner`

### What is parsed from GAF in your project?

Identity, coverage, path span, mapping quality, score, and traversed graph nodes.

### Why do you need both formats?

Because one stores the graph itself and the other stores alignments onto that graph.

## Common mistakes and confusions

- Saying GAF stores the graph
  No. GFA stores the graph.

- Saying GFA stores alignment scores
  No. GAF stores alignment results.

- Confusing graph path with query sequence
  The query is the read/sequence being aligned; the path is the route through the graph.

- Assuming one GAF file contains only one alignment
  It can contain multiple alignment records.

## Fast viva answer

`GFA stores the pangenome graph, while GAF stores a query-to-graph alignment. In our project, minigraph creates GFA files, PanAligner produces GAF files, and our parser extracts identity, coverage, score, and traversed nodes for evaluation.`
