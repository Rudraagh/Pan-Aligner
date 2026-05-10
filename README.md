# Implementation and Reproduction of the PanAligner Paper

This repository is organized as a paper-reproduction project centered on the real [PanAligner](https://github.com/at-cg/PanAligner) aligner and [minigraph](https://github.com/lh3/minigraph) graph construction workflow.

The final project focus is:

- pangenome graph construction
- cyclic graph handling
- co-linear chaining concepts
- graph alignment with PanAligner
- theoretical reproductions of paper ideas
- held-out alignment evaluation

The repository is intentionally presented as a PanAligner paper implementation project, not as a prediction or machine-learning system.

## Core scope

The main workflow now focuses only on paper-relevant components:

1. FASTA preprocessing
2. stratified train/test split
3. pangenome graph construction with minigraph
4. graph visualization and topology analysis
5. SCC detection
6. DFS-based back-edge removal
7. cyclic graph to DAG approximation
8. path cover generation
9. anchor representation
10. precedence relation
11. simplified co-linear chaining DP
12. gap-cost analysis
13. reachability analysis
14. iterative chaining and convergence demonstration
15. real PanAligner execution
16. GAF parsing
17. held-out alignment evaluation
18. theory and paper-reproduction reports

## Important distinction

The project has two clearly separated layers:

- `Real PanAligner implementation`
  This uses the actual PanAligner binary for sequence-to-graph alignment.

- `Simplified theoretical reproduction`
  This adds educational Python modules that demonstrate the main algorithmic ideas from the PanAligner paper for viva explanation, visualization, and reporting.

The theory modules do not replace, modify, or compete with the real PanAligner internals.

## Project structure

```text
FINAL YEAR PROJECT PHASE/
├── main.py
├── app_combined.fasta
├── psen1_combined.fasta
├── psen2_combined.fasta
├── PanAligner/
├── minigraph/
├── data/
│   ├── healthy/
│   ├── unhealthy/
│   ├── raw/
│   └── metadata/
├── graphs/
├── outputs/
│   ├── graphs/
│   ├── alignments/
│   ├── evaluation/
│   ├── theory/
│   └── reports/
├── scripts/
│   ├── common.py
│   ├── preprocess.py
│   ├── split_dataset.py
│   ├── build_graph.py
│   ├── align.py
│   ├── parse_gaf.py
│   ├── visualize.py
│   ├── paper_evaluation.py
│   └── theory/
│       ├── common.py
│       ├── graph_analysis_demo.py
│       ├── scc_demo.py
│       ├── dag_converter_demo.py
│       ├── path_cover_demo.py
│       ├── anchors_demo.py
│       ├── precedence_demo.py
│       ├── chaining_demo.py
│       └── visualization_demo.py
├── requirements.txt
├── run_pipeline.sh
└── README.md
```

Note:

- Some older exploratory scripts may still remain in the repository for archival reasons.
- They are not part of the main paper-reproduction workflow.
- `main.py` is the authoritative entry point.

## Setup

This project is intended for Linux or WSL.

### Automated setup

```bash
chmod +x scripts/setup_tools.sh run_pipeline.sh
./scripts/setup_tools.sh
```

### Manual setup

```bash
sudo apt-get update
sudo apt-get install -y build-essential zlib1g-dev python3 python3-pip

git clone https://github.com/at-cg/PanAligner
git clone https://github.com/lh3/minigraph.git

cd PanAligner && make -j"$(nproc)"
cd ../minigraph && make -j"$(nproc)"
cd ..

python3 -m pip install -r requirements.txt
```

## Build verification

PanAligner:

```bash
./PanAligner/PanAligner -cx lr PanAligner/test/MT.gfa PanAligner/test/MT-orangA.fa > out.gaf
head out.gaf
```

minigraph:

```bash
./minigraph/minigraph --version
```

## Single entry point

`main.py` is the master driver for the project.

Supported modes:

```bash
python3 main.py --full-pipeline
python3 main.py --theory-only
python3 main.py --evaluate
```

Meaning:

- `--full-pipeline`: preprocess, split, build graphs, run theory modules, execute PanAligner on held-out queries, generate evaluation outputs, and write final reports
- `--theory-only`: run only the educational reproductions of paper concepts
- `--evaluate`: run the held-out PanAligner alignment evaluation without rerunning the theory suite

## Workflow overview

### 1. FASTA preprocessing

`scripts/preprocess.py`:

- parses the gene FASTA files
- validates DNA sequences
- normalizes sequence formatting
- separates reference and labeled sequences
- writes processed FASTA files and `preprocess_manifest.json`

### 2. Train/test split

`scripts/split_dataset.py`:

- creates a reproducible stratified split
- preserves the current project’s existing train/test structure
- writes:
  - `data/metadata/train_manifest.json`
  - `data/metadata/test_manifest.json`
  - `data/metadata/split_manifest.json`
  - `data/metadata/split_summary.txt`

The split is retained because it is useful for held-out PanAligner evaluation even though classification is no longer the project focus.

### 3. Pangenome graph construction

`scripts/build_graph.py`:

- builds healthy-only, unhealthy-only, and combined graphs per gene
- uses the locus reference as the graph base
- uses minigraph with parameters suitable for this dataset

Although the directory names still reflect the original labeled FASTA organization, the paper-reproduction workflow focuses on the graph construction and alignment behavior rather than label prediction.

### 4. Graph analysis and visualization

`scripts/visualize.py`:

- parses GFA structure
- computes graph statistics
- renders graph visualizations

Generated outputs go under `outputs/graphs/`.

### 5. Educational theory modules

The folder `scripts/theory/` provides simplified conceptual reproductions of major PanAligner paper ideas:

- `graph_analysis_demo.py`
  representative graph extraction, reachability context, cyclic demo graph preparation
- `scc_demo.py`
  Tarjan SCC detection
- `dag_converter_demo.py`
  DFS back-edge detection and DAG approximation
- `path_cover_demo.py`
  simplified greedy path cover on the DAG approximation
- `anchors_demo.py`
  simplified anchor representation `(vertex, [x..y], [c..d], weight)`
- `precedence_demo.py`
  simplified precedence relation
- `chaining_demo.py`
  simplified co-linear chaining DP, gap costs, chain reconstruction, convergence demo
- `visualization_demo.py`
  wrapper to generate all theory visuals

These modules produce educational outputs under `outputs/theory/`.

### 6. Real PanAligner execution

`scripts/align.py` wraps the real PanAligner binary:

```bash
python3 scripts/align.py \
  --graph graphs/app.combined.gfa \
  --query data/test/healthy/app/APP_H_3.fa \
  --panaligner-bin ./PanAligner/PanAligner \
  --output-gaf outputs/alignments/app_query.gaf
```

This executes the actual PanAligner command shape:

```bash
./PanAligner/PanAligner -t 4 -cx lr graph.gfa query.fa > out.gaf
```

### 7. GAF parsing

`scripts/parse_gaf.py` extracts alignment information such as:

- alignment score
- query span
- path span
- identity
- coverage
- mapping quality
- traversed graph nodes

### 8. Held-out alignment evaluation

`scripts/paper_evaluation.py` evaluates the reproduction pipeline by:

- aligning held-out test sequences against the train-derived combined graphs
- collecting best-alignment summaries
- reporting:
  - alignment rate
  - mean identity
  - mean coverage
  - mean MAPQ
  - mean alignment score
  - per-gene summaries

This evaluation is paper-aligned because it focuses on sequence-to-graph alignment behavior rather than downstream classification.

## Output directories

The main workflow produces:

- `outputs/graphs/`
  graph visualizations and graph statistics
- `outputs/alignments/`
  GAF files and alignment overlays
- `outputs/evaluation/`
  held-out PanAligner evaluation summaries
- `outputs/theory/`
  educational theory outputs and plots
- `outputs/reports/`
  final written reports

## Key generated files

Examples of final artifacts:

- `outputs/theory/scc_analysis.txt`
- `outputs/theory/scc_graph.png`
- `outputs/theory/back_edges.txt`
- `outputs/theory/dag_graph.png`
- `outputs/theory/path_cover.txt`
- `outputs/theory/path_cover_graph.png`
- `outputs/theory/anchors.csv`
- `outputs/theory/precedence_analysis.txt`
- `outputs/theory/chaining_results.txt`
- `outputs/theory/chaining_graph.png`
- `outputs/theory/iteration_log.txt`
- `outputs/theory/convergence_plot.png`
- `outputs/theory/gap_cost_analysis.txt`
- `outputs/evaluation/alignment_metrics.json`
- `outputs/evaluation/alignment_records.json`
- `outputs/evaluation/evaluation_report.txt`
- `outputs/reports/theoretical_reproduction_report.txt`
- `outputs/reports/final_project_summary.txt`

## Sample usage

### Run everything

```bash
./run_pipeline.sh
```

### Run only the theory layer

```bash
python3 main.py --theory-only
```

### Run only the held-out alignment evaluation

```bash
python3 main.py --evaluate
```

## Viva / report explanation

### What the project reproduces

The project reproduces the PanAligner paper at two levels:

- practical workflow reproduction using the real PanAligner binary
- conceptual reproduction using simplified educational graph algorithms

### What co-linear chaining means here

Co-linear chaining selects an ordered set of anchors that is consistent in query order and graph reachability, while balancing anchor weights and gap costs. The theory layer demonstrates this idea using a simplified dynamic programming implementation.

### How cyclic graphs are handled

The theory modules explicitly demonstrate:

- SCC identification
- back-edge removal
- conversion of a cyclic graph into a DAG approximation
- path cover generation on the DAG approximation

This mirrors the paper’s graph-preprocessing logic conceptually, while keeping the real PanAligner alignment engine untouched.

### Why the theory modules are separate

PanAligner’s internal implementation is optimized and more sophisticated than a classroom-scale Python script. The theory modules are therefore intentionally educational reproductions rather than replacements.

## Final framing

The repository should be understood as:

**Implementation and Reproduction of the PanAligner Paper**

not as:

- a machine learning project
- a disease prediction framework
- an AI pipeline
