# 03. Cyclic Graphs and SCC

## What it is

### Cyclic graph

A cyclic graph is a graph that contains a loop. That means if you start at some node and follow directed edges, you can eventually return to the same node.

### SCC

`SCC` means `Strongly Connected Component`.

An SCC is a group of nodes where every node can reach every other node by directed paths.

## Why it is needed

The PanAligner paper studies co-linear chaining on pangenome graphs, including cyclic graphs. Cycles are important because they make ordering harder.

In a DAG, ordering is easier because there are no loops.

In a cyclic graph:

- topological ordering is not directly possible
- a path can circle back
- chaining decisions become more complex

That is why SCC analysis is important.

## Simple intuition

Think of a cyclic graph like a road roundabout.

- if the road network has no loops, you always move forward
- if there is a roundabout, you can come back around

An SCC is like a strongly connected neighborhood where every road can lead back to every other road.

## Tiny example

```text
A -> B -> C
     ^    |
     |    v
     E <- D
```

Here `B, C, D, E` can form a strongly connected region if each is reachable from the others.

## Technical explanation

In directed graphs:

- a cycle is any directed closed walk
- an SCC is a maximal strongly connected subgraph

If a node is alone and does not loop back, it can still form an SCC of size `1`.

That is exactly why the APP graph in this repo has:

- `115` SCCs
- largest SCC size `1`

This means every node forms its own SCC, so the graph is effectively acyclic in that example.

## How our project uses SCC

The script `scripts/theory/scc_demo.py` implements `Tarjan's algorithm`.

### Purpose

- detect SCCs in the educational graph
- visualize SCC structure
- help explain cyclic graph decomposition in viva

### Important detail

The theory layer may add synthetic back edges to create an educational cyclic scenario, because the real locus graphs in the current dataset are often DAG-like.

That is not cheating. It is a teaching strategy.

## Tarjan's algorithm in simple words

Tarjan's algorithm does a DFS and keeps track of:

- discovery index
- lowlink value
- a stack of active nodes

When a node is found to be the root of an SCC, all nodes up to that point are popped from the stack as one component.

## Why Tarjan's algorithm is good

- efficient
- standard
- directly gives SCC groups

## ASCII intuition for lowlink

```text
If DFS goes:
A -> B -> C -> D

and D can reach B again,
then B, C, D may belong to the same SCC.
```

The `lowlink` value helps detect this back reachability.

## How the repo demonstrates cyclicity

In `scripts/theory/common.py`:

1. a representative subgraph is extracted
2. if it is acyclic, synthetic back edges may be added
3. SCC logic is then demonstrated on that educational copy

This is very important to explain honestly:

`The educational cyclic graph is derived from the real graph, but synthetic back edges may be added only for teaching the paper's cyclic case.`

## Real repo observation

From `outputs/theory/graph_analysis.txt`:

- original APP combined graph: acyclic
- representative subgraph: acyclic
- educational graph: synthetic back edges added

This means:

- the practical graph is real
- the cyclic teaching case is educational

## How SCC helps later stages

Once SCCs are known, we can:

- identify cyclic regions
- reason about where loops happen
- simplify cyclic handling
- build DAG-like approximations for downstream processing

## DAG relation

If we compress each SCC into one meta-node, the SCC-compressed graph is always a DAG.

This is a major theoretical idea:

- cycles exist inside SCCs
- the overall relationship between SCCs becomes acyclic

## How our project explains SCC visually

The file `outputs/theory/scc_graph.png` colors nodes by SCC.

This helps in viva because you can say:

- same color means same SCC
- different colors mean different SCCs
- a single-node SCC means no cycle around that node

## Common viva questions

### What is a cycle?

A directed loop where you can start from a node and return to it by following edge directions.

### What is an SCC?

A maximal set of nodes in which every node can reach every other node.

### Why are SCCs important?

Because they identify cyclic regions and help convert a hard cyclic problem into a more manageable form.

### What algorithm did you use for SCC detection?

Tarjan's algorithm in the theory demo.

### Does every SCC mean a large cycle?

No. A single node can also be an SCC.

### Did your real APP graph contain cycles?

The current APP graph statistics show it is DAG-like, so the theory layer adds synthetic back edges only in the educational copy to demonstrate the cyclic case from the paper.

## Common mistakes and confusions

- Confusing a cycle with an SCC
  A cycle is a loop; an SCC is a strongly connected region.

- Thinking SCC size `1` means no SCC exists
  It still exists. It is just a trivial SCC.

- Claiming synthetic edges were added to the real graph
  They were added only to the educational theory graph.

- Saying SCC detection aligns sequences
  It does not. It is a graph analysis step.

## Fast viva answer

`A cyclic graph contains loops, and an SCC is a set of nodes that are mutually reachable. In our project, SCCs are detected with Tarjan's algorithm in the educational theory layer so we can explain how cyclic regions are identified before DAG-style processing and chaining.`
