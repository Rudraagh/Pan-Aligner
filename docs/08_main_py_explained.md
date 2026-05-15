# 08. main.py Explained

## What it is

`main.py` is the single master entry point for the current paper-reproduction project.

If someone asks:

`Which file controls the project?`

the answer is:

`main.py`

## Why it is needed

Instead of manually running many scripts one by one, `main.py` organizes the whole workflow and offers clear modes.

## Supported modes

`main.py` supports three main modes:

```text
python main.py --full-pipeline
python main.py --theory-only
python main.py --evaluate
```

## What each mode means

### `--full-pipeline`

Runs:

- preprocessing
- train/test split
- graph building
- graph visualization
- theory suite
- PanAligner held-out evaluation
- report writing

### `--theory-only`

Runs only the educational theory demos:

- graph analysis
- SCC
- DAG conversion
- path cover
- anchors
- precedence
- chaining

### `--evaluate`

Runs the held-out PanAligner alignment evaluation on the current train/test setup.

## Main functions in `main.py`

### `default_fasta_inputs()`

Finds all root FASTA files matching:

```text
*_combined.fasta
```

These are the default starting input files.

### `build_training_artifacts(...)`

This is one of the most important functions.

It does:

1. preprocess FASTA files
2. split into train/test
3. build train graphs
4. visualize healthy, unhealthy, and combined train graphs

### `choose_theory_graph(graph_manifest)`

Chooses which graph to use for the theory suite.

Preference:

- use `APP` combined graph if available
- otherwise use the first gene alphabetically

### `run_theory_suite(graph_path)`

Calls all educational demo scripts and returns their results.

### `run_alignment_evaluation(...)`

Calls `evaluate_panaligner_workflow(...)` from `paper_evaluation.py`.

### `write_theoretical_report(...)`

Creates a human-readable summary of:

- SCC logic
- DAG conversion
- path cover
- anchors
- precedence
- chaining

### `write_final_summary(...)`

Writes the overall project summary report.

## Flow of `--full-pipeline`

```text
Read FASTA
   ->
Preprocess
   ->
Split into train/test
   ->
Build train graphs
   ->
Visualize train graphs
   ->
Run theory suite on chosen graph
   ->
Run held-out alignment evaluation
   ->
Write reports
   ->
Print JSON summary
```

## Flow of `--theory-only`

```text
Select theory graph
   ->
Run theory demos
   ->
Write theory report
   ->
Write final summary
```

## Flow of `--evaluate`

```text
Preprocess and split
   ->
Build train graphs
   ->
Run held-out PanAligner evaluation
   ->
Write final summary
```

## Important parameters

`main.py` also accepts:

- `--threads`
- `--split-seed`
- `--test-fraction`
- `--minigraph-bin`
- `--panaligner-bin`
- `--theory-graph`
- `--input-fastas`

## Why this structure is good

### Clean orchestration

The orchestration is centralized, so the project is easier to:

- run
- explain
- demo

### Clear user-facing modes

You can quickly say:

- full pipeline for everything
- theory only for viva explanation
- evaluate for held-out alignment testing

## Real vs educational separation inside `main.py`

`main.py` makes the distinction very clear:

- real pipeline uses `minigraph` and `PanAligner`
- theory suite uses the Python demo modules

This design is excellent for viva because it avoids confusion.

## Example viva explanation

`main.py` is the controller file. It first prepares the data and builds training graphs, then optionally runs the theory suite and held-out alignment evaluation depending on the selected mode. This makes the entire project reproducible from one entry point.

## Common viva questions

### Why did you create a single `main.py`?

To make the project reproducible, easier to run, and easier to present.

### Which mode is used to run everything?

`--full-pipeline`

### Which mode is best for viva demonstrations?

`--theory-only` is useful for concept explanation, while `--full-pipeline` shows the full project capability.

### Does `main.py` directly implement PanAligner internals?

No. It orchestrates scripts and external tools.

### Why choose APP as the default theory graph?

Because the code prefers APP when available, likely to keep demo behavior consistent and familiar.

## Common mistakes and confusions

- Thinking `main.py` contains all algorithms directly
  It mainly coordinates them.

- Forgetting that `--evaluate` still depends on prepared graph/data workflow
  Evaluation needs train/test and train graphs.

- Saying theory-only mode runs PanAligner
  It does not.

## Fast viva answer

`main.py` is the master driver of the project. It provides full-pipeline, theory-only, and evaluate modes, and coordinates preprocessing, train/test split, graph building, theory demos, real PanAligner alignment evaluation, and report generation from one place.`
