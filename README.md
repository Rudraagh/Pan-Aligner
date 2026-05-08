# Alzheimer’s Disease Prediction with PanAligner and Pangenome Graphs

This project is a research-oriented prototype that builds Alzheimer’s Disease (AD) sequence classification on top of the real [PanAligner](https://github.com/at-cg/PanAligner) graph aligner and a real [minigraph](https://github.com/lh3/minigraph) GFA construction workflow. It uses three AD-related gene FASTA files (`APP`, `PSEN1`, `PSEN2`), preprocesses healthy and disease-associated sequences, builds pangenome graphs, aligns new queries with PanAligner, analyzes graph topology, and predicts `HEALTHY` vs `UNHEALTHY`.

## What this repository contains

- Real PanAligner integration using the official command shape:
  - `./PanAligner -cx lr graph.gfa query.fa > out.gaf`
- Real minigraph-based graph construction:
  - `minigraph -cxggs reference.fa sample1.fa sample2.fa ... > graph.gfa`
- Python orchestration for preprocessing, graph building, alignment, parsing, scoring, and visualization
- A Linux/WSL shell runner for reproducible end-to-end execution

## Project structure

```text
FINAL YEAR PROJECT PHASE/
├── app_combined.fasta
├── psen1_combined.fasta
├── psen2_combined.fasta
├── PanAligner/
├── minigraph/                       # clone separately with setup_tools.sh or git clone
├── data/
│   ├── healthy/
│   ├── unhealthy/
│   ├── queries/
│   ├── raw/
│   └── metadata/
├── graphs/
├── outputs/
│   ├── graphs/
│   └── alignments/
├── scripts/
│   ├── common.py
│   ├── preprocess.py
│   ├── build_graph.py
│   ├── align.py
│   ├── parse_gaf.py
│   ├── predictor.py
│   ├── visualize.py
│   ├── pipeline.py
│   └── setup_tools.sh
├── requirements.txt
├── run_pipeline.sh
└── README.md
```

## Data assumptions

Your FASTA headers already contain class labels:

- `APP_H_*`, `PSEN1_H_*`, `PSEN2_H_*` -> healthy
- `APP_U_*`, `PSEN1_U_*`, `PSEN2_U_*` -> unhealthy
- `ref|...` -> locus reference sequence

The preprocessing step uses these headers directly and does not invent labels.

## Setup

This pipeline is intended for Linux or WSL, which matches your requirement that the project be runnable on Linux/WSL.

### Option 1: automated setup

```bash
chmod +x scripts/setup_tools.sh run_pipeline.sh
./scripts/setup_tools.sh
```

### Option 2: manual setup

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

## Pipeline overview

### 1. FASTA preprocessing

`scripts/preprocess.py`:

- parses all gene FASTA files with Biopython
- validates DNA alphabet (`A/C/G/T/N`)
- normalizes uppercase sequences
- removes duplicates within each gene/class bucket
- writes split FASTA files to:
  - `data/healthy/<gene>/`
  - `data/unhealthy/<gene>/`
  - `data/raw/`
- writes `data/metadata/preprocess_manifest.json`

Run:

```bash
python3 scripts/preprocess.py
```

### 2. Pangenome graph construction

`scripts/build_graph.py` builds three graphs per gene:

- `graphs/<gene>.healthy.gfa`
- `graphs/<gene>.unhealthy.gfa`
- `graphs/<gene>.combined.gfa`

It uses the locus reference as the base and then adds each labeled sample FASTA as a separate minigraph input, which is important because minigraph graph construction works best when each input assembly/sample is its own file.

Run:

```bash
python3 scripts/build_graph.py --threads 8 --minigraph-bin ./minigraph/minigraph
```

### 3. Graph analysis and visualization

`scripts/visualize.py`:

- parses `S` and `L` records from GFA
- computes:
  - node count
  - edge count
  - weakly connected components
  - strongly connected components
  - self loops
  - cycle detection
- saves graph statistics JSON and PNG figures under `outputs/graphs/`

Run:

```bash
python3 scripts/visualize.py
```

### 4. Sequence alignment with PanAligner

`scripts/align.py` wraps the real PanAligner binary:

```bash
python3 scripts/align.py \
  --graph graphs/app.combined.gfa \
  --query data/queries/query.fa \
  --panaligner-bin ./PanAligner/PanAligner \
  --output-gaf outputs/alignments/app_query.gaf
```

This ultimately executes the real PanAligner-style command:

```bash
./PanAligner/PanAligner -t 4 -cx lr graphs/app.combined.gfa data/queries/query.fa > outputs/alignments/app_query.gaf
```

### 5. GAF parsing

`scripts/parse_gaf.py` converts GAF output into structured Python objects.

Extracted fields include:

- alignment score (`AS` if present, otherwise residue matches)
- query span
- path span
- mapping quality
- traversed graph nodes
- identity
- coverage
- normalized score

### 6. Healthy vs unhealthy prediction

`scripts/predictor.py` uses a graph-based evidence comparison:

1. Infer the most likely gene by aligning the query to each gene’s combined graph.
2. Align the same query to that gene’s healthy graph and unhealthy graph.
3. Compute a composite class score from:
   - identity
   - coverage
   - mapping quality
   - normalized alignment score
4. Predict the class with the stronger graph alignment evidence.

This keeps the disease layer on top of the real graph aligner instead of replacing it with a synthetic classifier.

Run:

```bash
python3 scripts/predictor.py \
  --query-fasta data/queries/query.fa \
  --panaligner-bin ./PanAligner/PanAligner \
  --threads 8
```

If you already generated GAFs externally, you can score them directly without rerunning alignment:

```bash
python3 scripts/predictor.py \
  --query-fasta data/queries/query.fa \
  --gene APP \
  --healthy-gaf outputs/alignments/app_vs_healthy.gaf \
  --unhealthy-gaf outputs/alignments/app_vs_unhealthy.gaf \
  --combined-gaf outputs/alignments/app_vs_combined.gaf
```

Output:

- `outputs/alignments/prediction.json`
- class-specific GAF files
- combined-graph alignment PNG

### 7. End-to-end run

```bash
./run_pipeline.sh data/queries/query.fa
```

Optional explicit gene:

```bash
./run_pipeline.sh data/queries/query.fa APP
```

## Why there are class-specific graphs and combined graphs

The project produces both:

- combined per-gene graphs for topology analysis and visualization
- healthy-only and unhealthy-only graphs for prediction

This is deliberate. In practice, PanAligner gives alignment evidence on graph paths and segments, but the easiest reproducible way to convert that into disease classification without modifying PanAligner internals is to compare the same query against class-specific pangenome graphs built from the same locus. That preserves a real graph-alignment workflow while making the classification step explainable and reproducible.

## Alzheimer’s disease prediction logic

The disease layer is intentionally simple and defendable for a final-year prototype:

- each gene has a healthy graph and an unhealthy graph
- the same query is aligned to both
- the better-supported class wins

The predictor returns:

- `prediction`
- `confidence`
- `selected_gene`
- `matched_regions`
- `combined_graph_nodes`
- a plain-language explanation

## Viva / report explanation

### How PanAligner works

PanAligner is a sequence-to-graph aligner for cyclic and acyclic pangenome graphs. It follows a seed-chain-extend strategy:

1. find seed matches between the query and graph
2. chain anchors across graph structure
3. perform extension/alignment along promising graph paths

The project uses PanAligner as-is and does not reimplement its chaining or alignment core.

### What co-linear chaining means here

Co-linear chaining selects a high-quality ordered chain of anchors so that a read can be mapped through a graph in a biologically consistent order. In pangenome graphs, this is harder than linear genomes because variants create branches and some graphs can contain cycles. PanAligner’s contribution is that it supports exact co-linear chaining generalized to cyclic graphs.

### How cyclic pangenome graphs are handled

The GFA/rGFA graph is constructed with minigraph, and PanAligner then aligns reads directly to that graph. Because PanAligner supports cyclic graphs, it can reason over loops and revisited graph structure rather than requiring the graph to be flattened into a single linear reference.

### How the disease predictor is layered on top

The predictor is separate from the aligner:

- minigraph builds healthy/unhealthy graph references
- PanAligner computes graph alignments
- Python parses GAF and compares graph alignment quality
- the final label is assigned from alignment evidence

That separation is useful in a project review because you can clearly explain which part is established research software and which part is your disease-classification contribution.

## Expected outputs

- `data/metadata/preprocess_manifest.json`
- `data/metadata/graph_manifest.json`
- `graphs/*.gfa`
- `outputs/graphs/*.png`
- `outputs/graphs/*.stats.json`
- `outputs/alignments/*.gaf`
- `outputs/alignments/prediction.json`
- `outputs/alignments/*.png`

## Current environment note

PanAligner was cloned into this workspace, but the current Windows PowerShell session does not have a working GNU `make` toolchain for the repository’s default Linux build, and direct WSL invocation from this session returned an access-denied error. The provided scripts and README therefore target Linux/WSL explicitly, which is also the platform you requested for reproducible execution.

## Sample walkthrough

```bash
# 1. Build tools
./scripts/setup_tools.sh

# 2. Create a query FASTA
cat > data/queries/query.fa <<'EOF'
>query_1
ACGTACGTACGTACGTACGT
EOF

# 3. Run the full pipeline
./run_pipeline.sh data/queries/query.fa

# 4. Inspect the prediction
cat outputs/alignments/prediction.json
```

## Suggested next enhancement

If you want to improve biological interpretability after the baseline pipeline is running, the next best extension is to annotate known AD-associated variant coordinates per gene and project the PanAligner path hits onto those loci so the explanation can mention specific disease-associated regions instead of only graph-support scores.
