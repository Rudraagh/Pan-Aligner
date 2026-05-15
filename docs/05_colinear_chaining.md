# 05. Co-linear Chaining

## What it is

Co-linear chaining is the process of selecting an ordered sequence of good local matches, called anchors, so that together they describe a strong overall alignment.

## Why it is needed

One query does not usually match a graph perfectly in one single chunk.

Instead, alignment is often built from smaller good matches. Chaining helps combine them in a sensible order.

Without chaining:

- local matches may be scattered
- wrong anchors may be chosen
- the final alignment may be inconsistent

## Simple intuition

Imagine matching a sentence to a book.

You find small matching phrases:

- phrase 1 matches page 10
- phrase 2 matches page 11
- phrase 3 matches page 12

If they appear in the same order in both the query and the target, they form a good chain.

That is co-linear chaining.

## In graph language

Each anchor is a local match between:

- part of the query
- part of a graph node or graph path

The best chain should:

- follow increasing query order
- follow valid graph reachability
- maximize reward
- penalize unreasonable gaps

## Anchor representation in our project

In `scripts/theory/anchors_demo.py`, an anchor is represented as:

```text
(vertex, [x..y], [c..d], weight)
```

Where:

- `vertex`: graph node
- `[x..y]`: interval on the graph node
- `[c..d]`: interval on the query
- `weight`: how strong that anchor is

## Why anchors matter

Anchors are the building blocks of chaining.

Good chaining requires:

- useful anchors
- correct precedence between anchors
- a scoring rule

## Precedence in simple words

Anchor `A` should come before anchor `B` if:

- `A` appears earlier in the query
- `A` can reach `B` in the graph

This prevents impossible chains.

## Gap cost intuition

Even if two anchors are in order, the gap between them matters.

If the jump between them is too large or inconsistent, the chain should be penalized.

In simple terms:

- closer, more consistent anchors are better
- huge jumps are less believable

## Simplified DP used in this repo

The educational demo uses:

```text
dp[j] = best score of a chain ending at anchor j
```

Then for each possible earlier anchor `i`, it checks:

```text
candidate = dp[i] + weight(j) - gap_cost(i, j)
```

If that candidate is better, `dp[j]` is updated.

## Why this is dynamic programming

Because the best chain ending at one anchor depends on best chains ending at earlier anchors.

This reuses solved subproblems instead of recomputing everything.

## Query gap vs graph gap

Our demo separates:

- `query_gap`
- `graph_gap`

### Query gap

Distance between anchors on the query sequence.

### Graph gap

Distance between anchors in the graph.

### Total gap

```text
total_gap = query_gap + graph_gap
```

## Simplified graph gap in our project

In `scripts/theory/chaining_demo.py`, graph distance is approximated using shortest-path distance between anchor nodes, plus small prefix/suffix adjustments.

This is not the full PanAligner internals.

It is an educational version that is easy to explain.

## Important distinction

### Real PanAligner

- optimized internal data structures
- real production alignment engine
- more advanced than the Python demo

### Our educational chaining demo

- small anchor set
- simplified precedence
- simplified gap cost
- simple DP recurrence
- used for understanding and viva explanation

## Real repo result

From `outputs/theory/chaining_results.txt`:

- best chain: `A6 -> A7 -> A8`
- best chain score: `12.00`

### What that means

Among the demo anchors, the DP found that anchors `A6`, `A7`, and `A8` formed the best consistent chain after considering weights and gap penalties.

## ASCII chaining example

```text
Query:   [A1]---[A2]--[A3]

Graph:   n1 -> n2 -> n3

If:
A1 maps to n1
A2 maps to n2
A3 maps to n3

and order is preserved,
then A1 -> A2 -> A3 is a valid co-linear chain.
```

## Why cyclic graphs make chaining harder

In a DAG, order is easier because you cannot loop back.

In cyclic graphs:

- same vertex may appear again
- reachability becomes trickier
- ordering is less obvious

That is why the paper pays special attention to cyclic graphs.

## Iterative convergence demo

The repo also includes an iterative score propagation demo.

Why?

Because some paper ideas are easier to explain as repeated score updates until values stabilize.

This helps show:

- how scores can improve over iterations
- why convergence matters conceptually

## Common viva questions

### What is co-linear chaining?

It is the process of selecting an ordered, compatible sequence of anchors that together form the best alignment explanation.

### What is an anchor?

A local match between a query interval and a graph location, with an associated weight.

### Why is chaining needed?

Because alignments are often built from multiple local matches, not one perfect global match.

### What is the DP idea?

For each anchor, compute the best score of any valid chain ending at that anchor.

### What is a gap cost?

A penalty for inconsistency or distance between consecutive anchors.

### Did you implement full PanAligner chaining?

No. We implemented a simplified educational DP for explanation, while real alignment is still done by the actual PanAligner binary.

## Common mistakes and confusions

- Thinking chaining means simple sorting
  It is not. It is optimization under constraints.

- Ignoring graph reachability
  Query order alone is not enough.

- Thinking the demo DP equals the full paper implementation
  It is only a simplified teaching version.

- Forgetting gap penalties
  They are central to choosing realistic chains.

## Fast viva answer

`Co-linear chaining selects an ordered set of compatible anchors that preserve query order and graph reachability while maximizing score and penalizing gaps. In our project, we explain this with a simplified dynamic programming demo, while the real sequence-to-graph alignment is still performed by the PanAligner binary.`
