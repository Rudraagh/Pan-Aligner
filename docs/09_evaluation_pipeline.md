# 09. Evaluation Pipeline

## What it is

The evaluation pipeline measures how well held-out test sequences align to train-derived combined pangenome graphs using the real `PanAligner`.

## Why it is needed

If we only build graphs and never test them, we cannot show whether the workflow works in practice.

Evaluation answers:

- do held-out sequences align?
- how strong are the alignments?
- are the graphs useful for unseen data?

## Simple intuition

This is like building a map from training routes and then checking whether new routes can still be explained well by that map.

## Evaluation design in this project

### Training side

Train sequences are used to build the graphs.

### Testing side

Held-out sequences are aligned later as queries.

This avoids the trivial case where the query was already used to construct the graph.

## Workflow

```text
Preprocess
   ->
Train/test split
   ->
Build graphs from training sequences only
   ->
Take test sequences one by one
   ->
Align each query to the corresponding combined train graph
   ->
Parse GAF
   ->
Store best alignment record
   ->
Compute overall and per-gene metrics
```

## Evaluation script

The main evaluation logic is in:

`scripts/paper_evaluation.py`

## What it does internally

For each gene:

1. load the combined train graph
2. iterate over healthy and unhealthy test sequences
3. run `PanAligner`
4. parse the resulting GAF
5. choose the best alignment
6. record summary values

## Metrics used

### Alignment found

Whether at least one valid best alignment exists.

### Alignment score

Uses the parsed `AS` tag if available, otherwise residue matches.

### Normalized score

Alignment score divided by query length.

### Identity

Residue matches divided by alignment block length.

### Coverage

Alignment block length divided by query length.

### Mapping quality

MAPQ value from the GAF line.

### Path span

Span on the graph path.

### Matched nodes

Number of traversed graph nodes in the alignment path.

## Overall metrics in this repo

From `outputs/evaluation/alignment_metrics.json`:

- query count: `46`
- aligned query count: `46`
- alignment rate: `1.0`
- mean identity: `0.9995539276605313`
- mean coverage: `1.0000659121670492`
- mean MAPQ: `60.0`
- mean alignment score: `128845.0`

## Per-gene metrics in this repo

| Gene | Query count | Aligned | Mean identity | Mean coverage | Mean MAPQ |
|---|---:|---:|---:|---:|---:|
| APP | 10 | 10 | 0.999756 | 1.000048 | 60 |
| PSEN1 | 32 | 32 | 0.999437 | 1.000117 | 60 |
| PSEN2 | 4 | 4 | 0.999983 | 0.999705 | 60 |

## How to interpret these numbers

### Alignment rate = 1.0

All held-out test queries aligned successfully.

### Identity near 1.0

The aligned sequences are extremely similar to the matched graph paths.

### Coverage near 1.0

Almost the whole query is being aligned.

### MAPQ = 60

This suggests high-confidence alignments in the current output.

## Important viva point

These results do not mean:

- the graphs are universally perfect
- the aligner solves every graph problem
- the theory demo is identical to the production tool

They mean:

- for this dataset and current workflow, held-out queries align very well to the train-derived combined graphs

## Why combined graph is used

The combined graph contains more variation than healthy-only or unhealthy-only alone. That makes it a more complete target for held-out alignment.

## Why no classification metric is the main focus now

The current project framing is paper reproduction and alignment behavior, not disease prediction.

So the key evaluation is:

- alignment quality

not:

- clinical classification accuracy

## Outputs generated

- `outputs/evaluation/alignment_metrics.json`
- `outputs/evaluation/alignment_records.json`
- `outputs/evaluation/evaluation_report.txt`
- `outputs/evaluation/alignments/.../combined.gaf`

## Example record meaning

Each stored record includes:

- gene
- sequence id
- bucket
- query FASTA path
- GAF path
- whether alignment was found
- score
- identity
- coverage
- mapping quality
- path span
- matched nodes

So the evaluation is both summary-level and per-sequence.

## Strengths of this evaluation

- held-out testing
- real PanAligner usage
- per-gene breakdown
- multiple quality metrics
- saved raw alignment outputs

## Limitations

- current dataset is limited to three loci
- many real graphs are mostly DAG-like
- alignment success alone does not fully capture all algorithmic behaviors in the paper

## Common viva questions

### Why use held-out sequences?

To test alignment on unseen data and avoid evaluating on sequences already used to build the graph.

### Why evaluate on combined graphs?

Because combined graphs capture both healthy and unhealthy variation and are more comprehensive targets.

### What does identity near 1 mean?

It means the aligned portion of the query matches the graph path very closely.

### What does coverage near 1 mean?

It means almost the full query was aligned.

### What is alignment rate?

The fraction of queries that produced a valid alignment.

### Why is evaluation important in a reproduction project?

Because reproduction should demonstrate both conceptual understanding and practical behavior.

## Common mistakes and confusions

- Confusing alignment score with identity
  Score is an overall alignment measure; identity is fraction matched.

- Thinking coverage > 1 is impossible
  In practice, due to alignment definitions and rounding, reported values can be near or slightly above 1 in derived summaries.

- Saying the evaluation proves disease prediction
  It does not. It proves strong graph alignment behavior in this setup.

## Fast viva answer

`Our evaluation aligns held-out test sequences to train-derived combined pangenome graphs using the real PanAligner, parses the best GAF alignments, and reports alignment rate, identity, coverage, MAPQ, and score. In the current outputs, all 46 held-out queries aligned successfully with identity and coverage very close to 1.`
