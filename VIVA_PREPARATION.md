# Viva Preparation: Pangenome Graph Project

## 1. Project overview

### The 30-second description

This repository has two deliberately separate stories.  The retained baseline is a reproducible PanAligner/minigraph-oriented graph/alignment and downstream ML workflow.  The **proposed prototype** is an alignment-free, minimizer-resolution pangenome-graph workflow: FASTA assemblies become minimizer walks, the walks become a graph with one explicit path per assembly, path-derived structural consistency annotates the graph, anchors are ranked, and a bounded PG-SCUnK-style controller can rebuild the graph with neighbouring parameter settings.  The final proposed-only run covered APP (36), PSEN1 (124), and PSEN2 (14) assemblies on the CPU.  Its reported scores are structural coherence metrics, not biological accuracy or superiority over PanAligner.

### Motivation and vocabulary

A reference genome is one representative sequence.  It can omit insertions, alternative alleles, structural variation, and population-specific sequence carried by other haplotypes.  A **pangenome** represents sequence diversity across many assemblies.  A pangenome graph represents shared sequence as shared graph elements and alternatives as branches; an input assembly can be represented by a path through it.

Here, a *haplotype-resolved assembly* means an input FASTA is retained as a separately named observed assembly/path, with a label such as `HEALTHY` or `UNHEALTHY`.  This does **not** mean the code performs phasing: it consumes assemblies already provided in that form.  Graph construction is difficult because repeated sequence, variation, orientation, errors, scale, and choice of graph resolution all affect the result.

The prototype uses **k-mers** (length-`k` substrings) and **minimizers** (a representative low-hash k-mer selected from each sliding window).  Minimizers reduce the number of tokens compared with all k-mers, while preserving an ordered sketch.  The graph nodes are minimizer identities, not base-level variation graph segments.  This is attractive because it avoids pairwise alignment in this code path, but it also sacrifices base-level resolution.

Paths preserve which ordered minimizer walk came from each assembly.  Anchors are selected graph nodes intended to be informative under the implementation's score.  Quality matters because a graph can be connected yet have weakly supported or unstable structure.  Parameters (`k`, window size, support threshold, anchor budget and spacing) change graph granularity and anchor behaviour; feedback is useful because a fixed choice is not automatically good for every input set.

### Actual architecture and status

```text
AssemblyInput FASTA
  -> compute_minimizers() [deterministic CPU routine]
  -> build_proposed_pangenome_graph()
  -> HaplotypePath walks + nodes + observed transitions
  -> calculate_graph_path_consistency()
  -> build_graph_anchor_candidates()/select_informative_anchors()
  -> run_pg_scunk_feedback_loop()
  -> graph/report bundle + final validation
```

Implemented: CPU minimizers, minimizer-node graph building, explicit paths, private-node rescue, graph-derived scores, deterministic anchor selection, a bounded local feedback heuristic, artifacts, unit tests, and a completed three-gene proposed evaluation.  Simplified/prototype: all of those algorithms are small project-local implementations, not a full biological graph assembler or a reproduction of PanAligner internals.  GPU: minimizer extraction has a CUDA implementation with CPU/GPU parity and a minimizer-only speed-up measured on one laptop GPU (Section 2, Contribution A); graph construction and every later stage remain CPU-only.  Proposed/future: GPU validation on HPC hardware, faster host-side record handling, richer orientation/base-level graph representation, and biological benchmark evidence.  Not implemented: GPU graph construction, exact PanAligner algorithm reproduction, and proof of biological superiority.

## 2. Novelty analysis and safe claims

### Contribution A — GPU minimizer extraction

**Idea and difference.** A conventional CPU sketcher computes all windows serially; the architecture exposes a `backend` choice and detects CuPy.  **Code:** `scripts/proposed/minimizer_engine.py:compute_minimizers`.  **What is actually implemented:** `_compute_gpu_minimizers` runs two CuPy `RawKernel`s (SHA-1 64-bit prefix per k-mer; sliding-window minimum, earliest position on ties), compacts consecutive duplicates on the device and copies only the selected minimizers back.  `backend="gpu"` raises if CUDA is unusable; `auto` falls back to CPU and records why.  Graph construction, consistency, anchors and feedback always run on the CPU, and the graph records the backend per stage.  The final proposed evaluation in Section 10 was run with `--backend cpu`.

**Measured evidence (one laptop GPU only):** NVIDIA GeForce RTX 3070 Ti Laptop GPU, CuPy 14.2.0, CUDA 13.2 runtime / 13.3 driver, commit `f2d3530`, results in `outputs/gpu_validation/laptop/`.  GPU output equals CPU output on 400 seeded random cases (IUPAC symbols, k up to the kernel limit 55) and all 12 benchmark cases; the CUDA-gated unit tests, including graph-level parity, pass.  Minimizer-only end-to-end speed-up is 2.8–12.4×, and 2.8× / 3.2× / 3.1× on the APP / PSEN1 / PSEN2 references at k=11, w=6.  For inputs of 100 kb or more, device time is only about 3–7% of the GPU path; the rest is Python record construction, which caps the speed-up.  The GPU memory pool held at most 27 MB.  Not measured: process RSS, HPC GPUs, whole-pipeline time.

**Safe viva claim:** “Minimizer extraction has a CUDA implementation whose output matched the CPU reference exactly on a laptop RTX 3070 Ti, with a 3–12× minimizer-only speed-up; graph construction is still CPU-only and the GPU path has not been validated on HPC hardware.”  Do not claim pipeline or graph-construction speed-up, HPC GPU validation, or that the final proposed evaluation used the GPU.  Examiner attack: *“Why only ~3× on your genes if the kernels are so fast?”* Answer: “Device time is a few percent of the GPU path; building Python records dominates.  Reducing host-side object construction is the next step, and the phase breakdown in `results.json` shows it.”

### Contribution B — path-preserving, consistency-aware construction and anchoring

**Idea.** Rather than discard low-support minimizers, the constructor retains every observed minimizer walk and marks nodes below `min_node_support` as rescued.  It records one `HaplotypePath` per assembly and computes context-derived scores before anchor ranking.  **Code:** `graph_constructor.py`, `graph_consistency.py`, `anchor_selection.py`.  **Math:** node/path/edge formulas are in Section 6; anchor score is `0.35U + 0.35D + 0.30C` by default.

**Implemented:** minimizer identity nodes, support sets, ordered paths, observed directed edges, rescue flags, and graph-derived anchor features.  **Limitations:** no sequence alignment, no reverse-complement canonicalization, no explicit orientation inference, no base-level paths, and “uniqueness” means low path coverage—not biological uniqueness.  **Safe claim:** “The proposed code preserves observed minimizer paths and uses their structural context in scoring; it is an alignment-free prototype.”

Examiner attack: *“Is keeping private nodes noise preservation?”* Answer: “It can preserve errors as well as real divergence.  The implementation intentionally protects observed divergent walks from silent deletion, then marks them as rescued/private and exposes their count.  It is a recall-oriented policy, not proof that every retained node is biological variation.”

### Contribution C — PG-SCUnK-style quality feedback

`quality_feedback.py` implements a project-local, deterministic heuristic controller.  It diagnoses the weakest structural signal, proposes a small unseen parameter neighbourhood, rebuilds the entire graph for each candidate, retains best-so-far only after a meaningful improvement (`> 1e-4`), and stops under bounded conditions.  The name must be framed as **PG-SCUnK-style/project-local**, not an external canonical implementation unless an external specification is separately demonstrated.

**Safe claim:** “It demonstrates closed-loop parameter retuning against a stated structural objective.”  It does not learn parameters, guarantee a global optimum, validate a biological objective, or prove the objective is correct.  Examiner attack: *“Can it optimize a bad objective?”* Answer: “Yes.  The controller is only as valid as the quality function; bounded search and validation protect software behaviour, not biological truth.”

### Not novel / do not claim as ours

- k-mers, minimizers, graph paths, GFA, graph traversal, hashes, and CPU dictionaries/sets are established ideas.
- minigraph and PanAligner are external baseline tools/concepts, not proposed algorithms implemented here.
- The legacy alignment, feature extraction, ML, GUI, and educational theory modules are retained project components, not evidence that the proposed graph algorithm is PanAligner.
- High structural-consistency values are not accuracy, sensitivity, precision, or superiority.

## 3. Codebase map and execution traces

| File/module | Actual role and important interfaces | Viva focus |
|---|---|---|
| `proposed_methodology.py` | CLI; `--final-proposed-evaluation`, `--pg-scunk-stress-demo`, normal workflow options. | `argparse`, `sys.path` setup, mutually ordered branches. |
| `scripts/proposed/models.py` | Dataclasses: `AssemblyInput`, `MinimizerRecord`, `ProposedGraphNode`, `ProposedGraphEdge`, `HaplotypePath`, `AnchorCandidate`, `ProposedGraph`. | Typed in-memory contract between stages. |
| `minimizer_engine.py` | SHA-1-prefix stable hashing; CPU reference and CUDA minimizer extraction. | CPU/GPU parity, `int32` kernel arguments, duplicate suppression. |
| `workflow.py` | Reads manifests/FASTA; runs per-gene proposed workflow and writes a normal summary. | Input validation and lightweight default cap of 8 per bucket. |
| `graph_constructor.py` | Builds nodes/paths/edges, rescue metadata, quality, JSON/GFA reports, integrity validation. | This is the core proposed graph constructor. |
| `graph_consistency.py` | Mutates graph annotations with path-derived node/edge/path scores and emits CSV/JSON. | Circularity and score ranges. |
| `anchor_selection.py` | Builds/scorers candidates and greedily applies spacing. | `U,D,C`, deterministic sort key. |
| `quality_feedback.py` | Candidate representation, diagnosis, bounded feedback loop, reports. | Best versus final candidate. |
| `stress_experiment.py` | Builds controlled four-haplotype demo and invokes the real constructor/controller. | Demonstration, not biology. |
| `final_evaluation.py` | Runs all available APP/PSEN1/PSEN2 train inputs without cap, validates and separates baseline/proposed reports. | CPU forced; no resume logic. |
| `scripts/build_graph.py`, `align.py`, `evaluate.py`, `paper_evaluation.py` | Legacy/baseline minigraph/PanAligner construction/alignment/evaluation. | Separate alignment-oriented evidence. |
| `extract_features.py`, `train_ml.py`, `predictor.py`, `predict_ml.py` | Legacy graph-alignment feature and prediction layer. | Not part of proposed graph quality calculation. |
| `gui_app.py`, `main.py` | Existing desktop/primary project workflow. | Do not conflate GUI with proposed implementation. |
| `tests/test_proposed_*.py` | Unit/smoke coverage for each proposed layer. | Correctness checks, not biological validation. |

**Final-run path:** `proposed_methodology.main()` -> `run_final_proposed_evaluation()` -> `load_assembly_inputs(train_manifest, None)` -> for APP/PSEN1/PSEN2: `run_pg_scunk_feedback_loop(..., backend="cpu")` -> `build_proposed_pangenome_graph()` -> `write_proposed_graph_bundle()`/`visualize_graph()`/`write_pg_scunk_reports()` -> `validate_final_evaluation()` -> CSV/JSON/README/comparison.

**Normal-run path:** CLI -> `run_proposed_methodology()` -> manifest input/filter/cap -> same feedback loop -> per-gene artifacts including `*.pg_scunk_feedback.json`.  It uses `feedback.pop("graph")`; final evaluation instead explicitly selects `feedback["best_graph"]`.

## 4. Minimizer engine

`_stable_hash(kmer)` computes `int(sha1(kmer.encode("ascii")).hexdigest()[:16], 16)`: a stable, deterministic 64-bit-prefix integer, not Python's randomized `hash()`.  `_compute_cpu_minimizers(sequence, k, window_size, assembly_id)` validates positive parameters, returns `[]` if the sequence is shorter than `k`, materializes `(kmer, hash, position)` tuples, then slides windows of `window_size` k-mers.  Each window chooses the lexicographically deterministic minimum by `(hash, position, kmer)`.  Consecutive identical `(hash, position)` selections are suppressed.  If windows produce no minimizer (for example when the number of k-mers is less than the window size), it emits the single global minimum k-mer.

There is **no reverse-complement canonicalization** in this source: it hashes the forward substring exactly as sliced.  Thus reverse-complement-equivalent sequences are not normalized to the same token.  Repeated minimizers at different positions are allowed in records; only immediate repeated node IDs are compressed when a path is built.

Worked conceptual example (hash values illustrative, not actual SHA-1): sequence `ACGTA`, `k=3` gives k-mers at positions 0/1/2: `ACG(8), CGT(3), GTA(6)`. With window 2: window `[ACG,CGT]` selects `CGT@1`; window `[CGT,GTA]` selects the same `CGT@1`, so duplicate suppression emits only `CGT@1`.  With a different minimum in the second window, a second record would be emitted.

Time is `O(L)` k-mer generation plus an `O(w)` `min` for each of approximately `L-k-w+2` windows: roughly `O(Lw)` as written, not an optimized deque implementation.  It materializes all k-mers and records, so memory is `O(L)` per assembly before graph aggregation.  A real GPU implementation would need device-side encoding, hashing/window reduction, duplicate compaction, memory bounds, and CPU/GPU equivalence tests.

## 5. Graph construction and models

`AssemblyInput` carries ID, gene, label and raw sequence.  `MinimizerRecord` carries a token, hash, minimizer position, window start and assembly provenance.  `ProposedGraphNode` has node ID (`n1`, `n2`, ...), minimizer `sequence`, support, mean position, provenance, rescue flag and consistency fields.  `ProposedGraphEdge` is a directed `(source,target)` transition with supporting assemblies.  `HaplotypePath` stores the ordered node IDs/positions for one assembly plus derived fields.  `ProposedGraph` aggregates them.

`build_proposed_pangenome_graph()` requires nonempty input and exactly one gene.  For each assembly it calls `compute_minimizers`.  It groups records by minimizer string while counting each minimizer once per assembly (`seen_in_assembly`); node support therefore counts distinct assemblies, not record multiplicity.  Nodes are created in sorted minimizer-string order.  A node is `is_rescued` when `support < min_node_support`; it is not removed.

For each assembly, records are converted to node IDs.  Only an immediately repeated node is omitted from the path; a repeated minimizer separated by another node remains a repeated visit.  Adjacent path nodes create directed edge keys and each edge accumulates a set of assembly IDs.  Therefore edges are observed path transitions, not inferred alignment edges.  Metadata records shared (`support > 1`), private (`support == 1`), rescued and discarded counts; discarded observed minimizer/node counts are explicitly zero.

Toy graph: haplotype H1 minimizers `[A,B,C]`, H2 `[A,D,C]`.  The graph has nodes A/B/C/D, edges A->B->C and A->D->C, and two paths.  If D occurs only in H2 and `min_node_support=2`, it remains as a rescued private node.  That protects the H2 walk, but cannot distinguish a genuine variant from a sequencing/assembly artefact.

GFA writes `S` lines with sequence, `RC` support and `CS` consistency; `L` lines are always `+` to `+` with overlap `0M`; `P` lines list `node+` segments and use `HP:Z` for the input label.  Hence it is forward-only orientation and uses `0M`, not a computed base overlap.  JSON serializes all dataclasses; large JSON is a full internal representation, not an external standard.

`validate_proposed_graph()` checks duplicate assembly paths, equal node/position list lengths, sorted positions, node existence and each path transition’s edge existence.  It does not prove biological validity, alignment correctness, coordinate accuracy, or absence of cycles.

## 6. Haplotype consistency: exact implemented formulas

`calculate_graph_path_consistency(graph)` is explicitly based only on observed `HaplotypePath` traversals.  Every score is intended to lie in `[0,1]`.

For node `v`:

```text
C(v) = .25 path_coverage + .25 predecessor_consistency
     + .25 successor_consistency + .25 local_path_consistency
```

`path_coverage` is the fraction of all paths containing `v`.  Predecessor/successor consistency are the dominant fraction among predecessor/successor contexts, using `START`/`END` at boundaries.  Local consistency averages, over node occurrences, the fraction of expected adjacent transitions that exist in the graph; singleton path nodes receive 1.0 because there is no expected transition to violate.  A private but perfectly locally connected node can therefore have a respectable context score despite low coverage.

For edge `(u,v)`:

```text
C(u,v) = .50 path_support_fraction + .50 source_transition_agreement
```

`path_support_fraction = edge_support / number_of_paths`; `source_transition_agreement = edge_support / number_of_paths_containing_source`.  The latter asks what share of source-node paths take that outgoing transition.

For a path `p`:

```text
C(p) = .40 valid_transition_fraction
     + .30 mean_node_consistency + .30 mean_edge_consistency
```

No-transition paths get transition fraction 1.0; an empty list of valid edges defaults to mean edge consistency 1.0.  It records rescued and low-score (`< .50`) node counts.

For the graph:

```text
C(G) = .40 mean_node_consistency
     + .30 mean_edge_consistency + .30 mean_path_consistency
```

Reports are `node_consistency.csv`, `edge_consistency.csv`, `path_consistency.csv`, and `graph_consistency.json`.  The legacy pairwise haplotype scorer is separate from this proposed construction flow; it is not the mechanism that produces nodes/edges or selects anchors.

**Critical circularity statement:** the edges were constructed from the observed paths, then path/edge validity is measured against those edges.  High values therefore show internal path/graph coherence under the construction, not independent evidence of biological correctness, variant truth, or alignment accuracy.  The source itself states this limitation in `quality_feedback.py` reports and `final_evaluation.py` README.

## 7. Anchor selection

For every graph node the code creates one `AnchorCandidate`.  Let `P(v)` be node path coverage, `occ(v)` its path occurrence count, and `maxocc` the largest such count:

```text
U(a) = 1 - P(v)                         (path selectivity)
D(a) = occ(v) / maxocc                  (relative path-occurrence depth)
C(a) = C(v)                             (Step-2 node consistency)
Score(a) = alpha U + beta D + gamma C
```

Defaults are `alpha=.35`, `beta=.35`, `gamma=.30`; `validate_anchor_weights` rejects negative values and values that do not sum to 1 within `1e-9`.  The values are user-chosen defaults in source, not learned from a training procedure.  Example: `U=.8,D=.6,C=.9` under `.4,.3,.3` gives `.4*.8 + .3*.6 + .3*.9 = .77`.

Candidates are ranked deterministically by `(-final_score, -consistency_score, -depth, node_id, anchor_id)`.  Greedy selection rejects a candidate if the absolute difference in **mean minimizer position** from any accepted anchor is less than `min_spacing`, stops at `max_anchors`, then returns selected anchors ordered by `(mean_position,node_id,anchor_id)`.  Ties therefore do not rely on dictionary iteration.

Complexity is `O(V + total_path_visits + V log V + VA)` where `V` is node/candidate count and `A` is accepted anchors (spacing scan); with small anchor budgets, the last factor is bounded.  `U` is not biological specificity, `D` depends on this dataset/path representation, and `C` inherits the structural circularity.  Anchor coverage is the fraction of input assembly IDs appearing in any selected anchor’s `haplotypes`, so it can be low even when graph consistency is high.  Because `U` favours private nodes that cover one haplotype each, coverage mostly tracked the `max_anchors` cap (APP 12/36 haplotypes).  Since commit `22076bc` anchor coverage is **reported only**: it is not part of `Q`, not diagnosed and not an acceptance condition, because no aligner consumes the selected anchors.

## 8. PG-SCUnK-style controller

`run_pg_scunk_feedback_loop()` starts from `_adaptive_default_parameters`: `k` is clamped to 3..11 from half the shortest sequence; window is clamped to 2..6; minimum support is 2 for at least three assemblies; anchor count is clamped to 4..12 based on assembly count; spacing is 4.0 for shortest sequence >=40 else 2.0.  Overrides may replace these values.

`canonical_configuration()` makes an immutable ordered tuple over `k, window_size, min_node_support, max_anchors, min_anchor_spacing, anchor_alpha, anchor_beta, anchor_gamma`; float values are rounded to 12 decimals.  `visited` prevents recomputing the same configuration.  Each candidate actually calls `build_proposed_pangenome_graph`, so it rebuilds the graph rather than re-labeling an old score.

Diagnostic signals are: node consistency, node support ratio, edge consistency, path consistency, and `1 - min(1, branching_node_ratio)`.  The smallest signal (with fixed priority tie-breaking) chooses an action label.  Every action proposes the same neighbourhood: `k - 2`, `window - 1`, or both (floored at 3 and 2).  Only `k` and `window_size` are searched: `min_node_support` only labels rescued nodes (no node is removed, so it never changes `Q`), and anchor settings are not scored.  The other fields stay in the configuration key and are still passed to the builder.  This is a small deterministic neighbourhood—not Bayesian optimisation, gradient descent, or exhaustive search.

```text
pending = [initial]; visited = {}; best = -infinity
while pending and round < max_rounds:
  evaluate each unvisited candidate by rebuilding graph
  record quality and acceptance
  replace best only if Q > best + tolerance
  break if acceptance passes
  stop on no records, stagnation, or no unseen candidates
  diagnose best graph and make unseen neighbours
return protected best graph, plus final candidate record and trajectory
```

Acceptance: `_quality_passes` requires overall quality >= `.72`, node consistency >= `.45`, and path consistency >= `.45` by default.  (Before commit `22076bc` it also required anchor coverage >= `.70`; that gate was removed with the anchor term.)  Stopping reasons are `quality_threshold_reached`, `max_rounds_reached`, `no_unseen_configurations`, or `stagnation_patience_exceeded`; default maximum rounds is 3 and default stagnation patience is 2.

Best-so-far protection matters because the final evaluated candidate can be worse.  `best_graph` is returned separately from `final_graph`; final evaluation serializes the best graph.  A duplicate candidate is skipped; an immediate threshold success ends after its first accepted evaluation; worse candidates remain in the audit trail but do not replace best.

### Stress-demo trace

`outputs/proposed_step5_pgscunk_demo/experiment_metadata.json` labels the four short controlled-substitution haplotypes **DEMONSTRATION ONLY — not a biological benchmark**.  Regenerated at commit `f84cb42` with `python main.py --proposed --pg-scunk-stress-demo --backend cpu` (current structure-only formula).  It began at `Q=0.70938`, target `.99`, `k=6, window=5`.  Seven unique configurations were rebuilt in three rounds, varying only `k` and `window`.  The best was candidate `c004`, round 2: `k=4, window=4`, `Q=0.76872`; round 3 tried `k=3` and/or `window=3` and all were worse (the final candidate `c007` had `Q=0.66806`).  It stopped at `no_unseen_configurations`, not success.  This demonstrates rebuilding, diagnostics and best-so-far retention; it does not validate a gene or disease claim.  It is currently the only run in which the controller actually retunes anything (see Section 10).  (With the old formula, the same demo had reached `Q=0.80341` at `c006`.)

## 9. Exact quality function and validation

The source in `graph_constructor.py` calculates:

```text
Q = ( .25 mean_node_consistency + .20 mean_edge_consistency
    + .15 mean_path_consistency + .15 mean_node_support_ratio
    + .10 (1 - min(1, branching_node_ratio)) ) / .85
```

`mean_node_support_ratio` averages `node.support / assembly_count`; `branching_node_ratio` is nodes with more than one outgoing edge divided by nodes.  This differs from graph consistency `C(G)`, which uses only node/edge/path consistency `.40/.30/.30`.  Until commit `22076bc`, `Q` also had a `.15 anchor_coverage_ratio` term; it was removed and the remaining original weights (sum `.85`) were rescaled to sum to 1, keeping their proportions.  Anchor coverage is still reported in `quality.json`.  Two inputs are constant by construction (path transition fraction and node local path consistency are always 1, so every path score is at least `.40`), which keeps `Q` high; `tests/test_metric_definitions.py` pins this.  `validation.json` is a software/data-integrity check: paths are valid, score fields are in range, PG-SCUnK keys are unique, selected anchor positions obey spacing, best tracking is consistent, and CSV/JSON report counts/quality agree.  Validation is not biological correctness, statistical significance, or benchmark accuracy.

## 10. Final real-data evaluation

Source: `outputs/proposed_final_evaluation/proposed_final_summary.json` (regenerated at commit `22076bc` with `python proposed_methodology.py --final-proposed-evaluation --backend cpu`, 194 s); each per-gene `validation.json` passed.  The final evaluator uses all listed train-manifest inputs and the requested minimizer backend (CPU in this run).

| Gene | Assemblies | Nodes | Edges | rescued/private | anchors | node / edge / path C | graph C | anchor cov. (reported) | Q | rounds/configs | accepted | backend |
|---|---:|---:|---:|---:|---:|---|---:|---:|---:|---|---|---|
| APP | 36 | 72,646 | 79,053 | 116 / 116 | 12 | .98248 / .99867 / .98755 | .98886 | .33333 | .98399 | 1 / 1 | yes | cpu |
| PSEN1 | 124 | 21,397 | 23,118 | 417 / 417 | 12 | .98380 / .98488 / .98879 | .98562 | .08871 | .97999 | 1 / 1 | yes | cpu |
| PSEN2 | 14 | 15,635 | 16,089 | 46 / 46 | 12 | .98929 / .98891 / .99217 | .99004 | 1.00000 | .98648 | 1 / 1 | yes | cpu |

All three genes pass on their initial configuration, so **the controller performs no retuning on the real inputs**; its retuning has only been exercised on the synthetic stress demo.  The `.72` threshold is far below these scores, which are high largely because of the constant-by-construction terms (Section 9) and because the inputs are closely related locus sequences (the held-out haplotypes in `outputs/truth_eval/REPORT.md` differ from the reference by only one or two small variants).  The graphs themselves are unchanged from the earlier run (same node and edge counts).  In the earlier run (old formula), APP (`Q = .88639`) and PSEN1 (`Q = .84871`) were rejected only by the anchor-coverage gate, which mainly reflected the anchor cap; PSEN2 (`Q = .98851`) passed.  Differences between the three genes are small and carry no biological meaning.

The run covered 174 assemblies.  It validates serialisation and structural properties, not a held-out biological ground truth.  The baseline alignment metrics in `outputs/evaluation/` are intentionally separate and not numerically interchangeable with this table.

## 10b. Machine-learning layer (full detail: `docs/17_ml_layer.md`)

**Model A, alignment reliability.**
- **Data and labels:** labels come only from simulation truth. PBSIM2 reads from the 46 held-out haplotypes are aligned by PanAligner and minigraph to their own gene's graph and to the other two genes' graphs as decoys. That gives 113,295 correct, 14,073 incorrect (13,977 off-target, 96 within-locus) and 259,920 unaligned rows.
- **Features:** only things visible after alignment.
- **Evaluation:** folds over 20 kb windows of true position with neighbouring windows removed, plus a leave-one-gene-out test; isotonic calibration on a held-back part of each training fold; a MAPQ-only model and a plain MAPQ threshold as baselines.
- **Results:** PR-AUC for the incorrect class is 0.979 (gradient boosting) and 0.986 (logistic) vs 0.564 (MAPQ only) under window folds, and 0.929 and 0.956 when a gene is held out. The first figure of 0.995 used a haplotype split that was too weak, because every haplotype is the same reference plus one variant.
- **Compared with MAPQ 60:** at the same 89.6% coverage, error among accepted is 0.84% (gradient boosting, window folds) vs 4.24%.
- **Focused second dataset (`outputs/ml_alignment_v2/`):** 75% accuracy reads with 536 within-locus errors, 96% of them minigraph's. Gradient boosting's within-locus PR-AUC is 0.965 (window) and 0.946 (gene held out). Under strict evaluation with tie-aware comparison, the model is reliably better than MAPQ 10–40 (at 94.6% accepted, 0.16–0.17% vs 0.30% error), indistinguishable at 50, and unresolved at 60 (0.17% vs 0.13%; intervals include zero).
- **Abstention (v2, strict evaluation):** coverage targets open the band in every fold and abstained alignments are enriched for errors (7–76% vs 0.83%), but error targets of 0.2% or lower are not reached on held-out data and over half of rows get p(incorrect) = 0 while about 0.25% of them are wrong, so the three-way output is not dependable.
- **Abstention (first dataset):** with the 1% / 90% targets the "insufficient evidence" band never opened. With 0.1% / 99% it covers about 2% of alignments, and held-out error among accepted is 0.12–0.17%, above the 0.1% target; that run used only the weak haplotype split.

**Model B, variant classification experiment: a negative result.** Full account: `docs/17_ml_layer.md`. Evidence: `outputs/model_b_training/20261030_training_01/`.

- **What it investigated:** whether simple sequence-derived features recover ClinVar classification direction (P/LP against B/LB) in APP, PSEN1 and PSEN2 better than simple baselines.
- **Dataset:** 693 variants from ClinVar at 1 star or more: 162 P/LP, 531 B/LB, 186 missense (153 P / 33 B). `VCV001231596` is excluded because its REF equals its ALT. The freeze still has 694 rows; the final data is `20261030_validation_02`, and `20261030_validation_01` is superseded.
- **Effective sample:** 13 clusters hold pathogenic variants, and only 12 among missense.
- **Why it is not an Alzheimer's classifier:** labels classify variants in three genes, not people. There are no patients, phenotypes or outcomes, and the disease assertions are not Alzheimer-specific. Benign does not mean protective.
- **Features:** variant type, length change, transition flag, strand-neutral substitution, trinucleotide context, CpG context, GC content and homopolymer run. For missense: amino-acid change, BLOSUM62, change in hydropathy, volume and charge, and proline, glycine and cysteine flags. Gene, position and consequence are baselines only.
- **Evaluation:** clusters (1 kb chains plus same residue) never cross roles; 5 folds, 5 seeds; settings and thresholds chosen on validation; 95% intervals from resampling clusters.

| Experiment | Balanced accuracy |
|---|---:|
| Missense model | 0.734 |
| Best baseline, chosen on validation | 0.848 |
| Missense model without GC | 0.622 |
| Amino-acid features alone | 0.582 |
| All-variant model | 0.702 |
| All-variant model without GC | 0.578 |
| Consequence-only baseline | 0.969 |

Primary difference −0.114. Headline seed 0: −0.177 (95% interval −0.311 to −0.021). Seeds 1–4 run from −0.139 to −0.078, with intervals including zero. So: no demonstrated gain.

**Why the baselines won.** They exploit how the dataset is composed. PSEN1 missense is 125 P to 4 B, so gene alone scores 0.848. Benign records are mostly synonymous, intronic or UTR, so consequence scores 0.969. The sequence model was not given either. Neither baseline is a biological predictor of an individual variant, and gene may reflect curation and reporting differences.

**What GC analysis revealed.** Removing GC lowers the mean from 0.734 to 0.622. `gc_200` has the largest coefficient (−0.72, negative in all 25 folds), 88% of its variance lies between clusters, and its direction reverses by gene (AUC 0.05, 0.63, 0.78). It appears to carry regional or gene-associated signal, mostly through APP. On seed 0 the paired difference is only 0.023 (−0.108 to 0.176). This does not prove sequence holds no biological information.

**Why it is still useful.** It stops an inflated claim. The protocol could have shown a gain, and it did not. It shows that leakage-aware design matters, that strong-looking features can be proxies for region, and that any future model must beat gene and consequence. It also shows how little independent data these three genes offer.

**Possible improvements, each a new declared analysis:** a later frozen ClinVar release with at least 20 P/LP and 20 B/LB in new clusters; more genes; protein-domain or conservation features, watching for circularity; REVEL, CADD or AlphaMissense as comparisons only; control for PSEN1/PSEN2 paralogy.

**Likely viva questions:**

- *"Did Model B fail?"* "The experiment worked and its answer was negative: these features didn't beat simple baselines under a fixed protocol. That's a finding about this representation and dataset."
- *"Can it predict Alzheimer's or someone's risk?"* "No. It only tries to recover a database annotation of variants in three genes."
- *"Why not just report the gene baseline's 0.848 as your model?"* "It reflects PSEN1 having 125 pathogenic and 4 benign missense variants. It says nothing about a given variant's biology, and may reflect how genes were studied and reported."
- *"Why use balanced accuracy?"* "Classes are imbalanced, 162 to 531 overall and 153 to 33 in missense, so plain accuracy would reward predicting the majority."
- *"Why group by cluster?"* "Nearby variants and those on one residue are not independent. Random splits let a model recognise a region. That's why 12 pathogenic missense clusters, not 153 variants, sets the real sample size."
- *"Does losing accuracy without GC prove GC matters biologically?"* "No. It points the other way: GC mostly marks region and gene. It's an inference from the pattern, not a proof, and the seed 0 paired interval includes zero."
- *"So sequence has no information about pathogenicity?"* "We can't say that. These features are simple and the data are few."
- *"Could you improve it by tuning now?"* "Not as confirmation. We've read the test folds, so any further run on them is exploratory. Independent evaluation needs new variants or a later frozen release."
- *"Holding out APP gave 0.660 against 0.500. Doesn't that show generalisation?"* "No. The baseline can't exceed 0.5 there by construction, the model includes GC, and three genes can't show generalisation. Other held-out genes gave about 0.47–0.52."
- *"PSEN2 amino-acid features reached 0.736. Isn't that a finding?"* "No. It's an after-the-fact observation on 14 variants."
- *"Is 693 your sample size?"* "Only as a row count. Thirteen clusters hold pathogenic variants, and 64% of records rest on one submitter."
- *"What happened to the original 220-sequence dataset?"* "Its label source is undocumented, though an audit found 216 of 220 match ClinVar records. We replaced it with documented data and kept the old study as history."

**Safe claim:** "Model A detects off-target alignments on simulated, locus-scale reads far better than MAPQ alone. Under cluster-separated evaluation, Model B's sequence-derived features did not recover ClinVar classification direction better than gene or consequence alone, in APP, PSEN1 and PSEN2. It carries no claim about diagnosis, individual risk or Alzheimer's disease."

Examiner attack: *"Your PR-AUC of about 0.98 looks too good."* Answer: "Our first figure of 0.995 was optimistic because of a weak split, which we replaced. Almost every incorrect alignment is an off-target read aligned to a single-locus graph, which low identity and coverage give away. It says nothing about whole-genome misplacement; only 96 within-locus errors exist, too few for strong claims."

## 11. Testing evidence

| Test file | Coverage / important cases | It catches | It does not prove |
|---|---|---|---|
| `test_proposed_minimizer_engine.py` | deterministic records, parameter validation, short inputs, backend metadata; CUDA-gated parity (incl. 120 seeded random cases), `int32` kernel arguments, phase timings. | unstable/minimum/empty-input regressions; CPU/GPU divergence on a CUDA machine. | speed-up, or parity on GPUs other than the one tested; sequencing correctness. |
| `test_proposed_graph_constructor.py` | node/edge/path construction, ordered paths, private-node rescue, validation and reports. | missing path edges or accidental low-support deletion. | biological graph accuracy. |
| `test_proposed_graph_consistency.py` | node/edge/path/graph formulas and reports. | formula/report regressions. | independence from construction circularity. |
| `test_proposed_anchor_selection.py` | weight validation/formula, ranking effects, deterministic spacing, reports. | incorrect score/ranking/spacings. | biological anchor specificity. |
| `test_proposed_haplotype_scoring.py` | legacy haplotype scorer behaviour. | legacy scorer regressions. | proposed graph construction quality. |
| `test_proposed_quality_feedback.py` | candidate generation, visited keys, improvement/best retention, stopping conditions. | duplicate/replacement/control-flow errors. | optimal tuning. |
| `test_proposed_pg_scunk_step4.py` | Step-4 loop/report semantics. | broken controller integration. | external PG-SCUnK equivalence. |
| `test_proposed_pg_scunk_demo.py` | controlled stress demo/smoke evidence. | demo metadata/output regressions. | real-world benchmark performance. |

Tests demonstrate deterministic software properties on selected synthetic/small cases.  They do not constitute a full published benchmark, biological validation or memory scalability study.  CPU/GPU parity has been checked on one laptop GPU only (CUDA tests plus `scripts/gpu_minimizer_bench.py`).

## 12. Baseline, PanAligner, minigraph, and legacy ML

The retained baseline uses project scripts such as `build_graph.py`, `align.py`, `evaluate.py`, and `paper_evaluation.py` to orchestrate external minigraph/PanAligner binaries and parse GAF/alignment-oriented outputs.  `paper_evaluation.py` contains local minigraph-first logic and falls back/selects PanAligner under its stated conditions.  Feature extraction and ML (`extract_features.py`, `train_ml.py`, `predict_ml.py`) consume PanAligner-derived alignment features for downstream classification; they are not part of proposed `Q`.

The proposed constructor is not “just minigraph”; it is a separate Python minimizer-walk prototype.  Conversely, it is not a replacement claim for PanAligner, and it is not an exact reproduction of PanAligner internal algorithms.  Baseline files are preserved so a fair future comparison can use identical datasets, splits, hardware, time/memory measurement, fixed tool versions, stated parameters, and appropriate comparable outcomes (for example truth-backed variant/path/alignment measures), rather than comparing structural Q directly to alignment accuracy.

## 13. Limitations and future work

| Limitation | Impact/current-review wording | Realistic future work |
|---|---|---|
| GPU limited to minimizers, one device | Minimizer-only 2.8–12.4× on a laptop GPU; graph construction is CPU; final evaluation ran on CPU. | Faster host-side record handling, HPC GPU validation, GPU-backed graph stages if justified by profiling. |
| Forward-only/minimizer graph | Does not model reverse orientation or base-level overlaps. | Canonical/strand-aware tokens, oriented edges, sequence segments and overlaps. |
| `O(Lw)` window scan / materialization | Memory/time may not scale to chromosome-level assemblies. | Monotonic-deque minimizers, streaming/chunking, profiling. |
| Rescue retains all low support | Preserves divergence and potential errors alike. | Read/assembly quality evidence and support calibration. |
| Structural circularity | Scores are not independent biological truth. | Truth sets, alignment/variant/path benchmarks and orthogonal validation. |
| Heuristic controller | Small local search can miss a better configuration. | Larger constrained search, held-out tuning, multi-objective and robustness studies. |
| Fixed anchor weights/budget | Dataset-sensitive coverage; APP/PSEN1 show this. | Adaptive budgeting/learned weights with held-out evaluation. |
| No PanAligner reproduction claim | Cannot claim internal equivalence. | Reproduce only where specifications/tests permit and cite exact versions. |
| No full benchmark suite | No superiority/statistical claim. | External data, baselines, runtime/memory, uncertainty and repeats. |

For the current review, these limitations are not fatal if the project is defended as an implemented, tested structural prototype with clearly delimited claims.

## 14. Future work

Future work must be described as future work, not quietly folded into the implementation claim: (1) validating the GPU minimizer path on HPC GPUs and reducing its Python record-construction overhead; (2) canonical, orientation-aware paths and base/segment-level graph overlaps; (3) independent biological/path/variant validation; (4) adaptive anchor budgets and weights selected on held-out data; (5) stronger multi-objective parameter search; (6) external benchmark datasets and repeated runtime/memory experiments; and (7) a fair, versioned comparison against minigraph/PanAligner.  The expected payoff is stronger biological and systems evidence; the required discipline is to measure it rather than infer it from the current structural score.

## 15. Code walkthroughs

1. **Minimizers:** input `sequence,k,window,id` -> validate -> build k-mer tuple list -> select minimum tuple per window -> suppress same `(hash,position)` -> `MinimizerRecord` list.  Edge cases: short sequence and window longer than k-mer list.  Complexity approximately `O(Lw)`/`O(L)` memory.
2. **Graph:** input homogeneous-gene assemblies -> per-assembly minimizers -> string-grouped nodes -> directed transitions from paths -> consistency -> anchors -> quality -> integrity check.  Edge case: empty/mixed-gene input raises `ValueError`.
3. **Paths:** records map through `node_lookup`; immediate equal consecutive node IDs are skipped; positions remain ordered record positions.  Repeated nonconsecutive tokens remain meaningful revisits.
4. **Consistency:** build node context and observed-edge maps from paths -> annotate nodes -> edges -> paths -> graph summary.  Edge case: no path transitions defaults as described in Section 6.
5. **Anchors:** count path node occurrences -> derive U/D/C -> weighted score -> deterministic greedy spacing.  Edge case: no candidates gives zero-valued summary.
6. **Candidate generation:** diagnose lowest signal -> emit 2–3 parameter alternatives -> canonicalise -> remove visited/local duplicates.  No random sampling.
7. **Feedback:** evaluate all pending current-round candidates; protect best; possibly stop; generate next round around best.  `max_rounds` bounds total rounds, not total candidates.
8. **Quality:** compute all six components in constructor, then `_quality_passes` separately applies stricter gates.  Q and acceptance are different.
9. **Final evaluation:** uses uncapped manifest data, serial genes, CPU, writes per-gene artifacts first and top-level summaries after all genes succeed.  It has no resume/skip-completed logic.

## 16. Likely “explain this code” prompts

| Short source pattern | Explanation / possible modification answer |
|---|---|
| `@dataclass(frozen=True)` | Immutable value-like input/record; mutation would raise.  Use frozen only where identity should not change. |
| `default_factory=list` | Gives each instance its own list; avoids shared mutable defaults. |
| `sha1(... )[:16]` | Deterministic 64-bit-prefix hash integer; collision risk remains theoretically possible. |
| `min(window, key=lambda...)` | Selects deterministic minimizer with explicit tie-breaks. |
| `last_signature` | Suppresses only consecutive duplicate selections across overlapping windows. |
| `backend.lower()` | Normalizes user input before validation. |
| `_load_usable_cupy()` | Proves a CUDA device can allocate before the GPU path is used; the kernels themselves are in `_compute_gpu_minimizers`. |
| `defaultdict(list)` | Creates grouping list on first key access. |
| `seen_in_assembly` | Prevents repeated occurrences inflating node support within one assembly. |
| `sorted(grouped_records.items())` | Stable node-ID assignment across runs. |
| `support < min_node_support` | Classification as rescued, not filtering. |
| `node_lookup[minimizer]` | Maps token identity to created graph node. |
| `if not ordered_nodes or ...` | Avoids indexing empty path and adjacent self duplication. |
| `zip(nodes, nodes[1:])` | Produces adjacent directed transitions. |
| `set[assembly_id]` support | Counts distinct paths, not number of observations. |
| `validate_proposed_graph` | Structural/integrity validation only. |
| `Counter(values)` | Obtains context mode for dominant-fraction score. |
| `START` / `END` | Makes boundary context explicit in consistency calculation. |
| `.50*... + .50*...` | Weighted convex score; interpretation needs each component. |
| `< LOW_CONSISTENCY_THRESHOLD` | Flags diagnostic counts, not deletion. |
| `math.isclose(sum,1)` | Prevents invalid anchor weight mixtures despite floating point. |
| `-final_score` sort key | Python ascending sort thereby ranks larger scores first. |
| `abs(pos-pos)<spacing` | Greedy coexistence rule based on mean positions. |
| `round(float,12)` key | Stabilizes float configuration identity. |
| `visited` set | Prevents repeated rebuilds/infinite cycling through same candidates. |
| `Q > best + tolerance` | Rejects negligible numerical/noise changes as meaningful improvements. |
| `final_graph` vs `best_graph` | Audit’s last candidate may be worse; returned result protects best. |
| `csv.DictWriter` | Enforces an explicit, reproducible report column order. |
| `P\t...node+` | GFA path serialisation is forward orientation in this prototype. |
| `--backend cpu` final run | The Section 10 results were produced on the CPU; GPU evidence is separate, in `outputs/gpu_validation/laptop/`. |

## 17. Viva question bank

Use the answer pattern below: **short** is the opening sentence; **deeper** is the defensible expansion; **file** tells you where to point; **follow-up** anticipates pressure; **safe wording** avoids overclaiming.

| Area / question | Short + deeper answer | File / math / follow-up / safe wording |
|---|---|---|
| A. What is a pangenome? | A representation of sequence diversity across multiple genomes/assemblies. A graph can share common material and branch at alternatives. | Concept; follow-up “is every graph biologically correct?” Safe: “It represents inputs at chosen resolution.” |
| B. Why not one reference? | One reference may miss alternative sequence and variation. This project retains one path per supplied assembly. | `models.py`; follow-up phasing. Safe: “Inputs are already assembly-resolved.” |
| C. Is this a variation graph? | It is a directed minimizer-node graph, not a base-level, fully oriented variation graph. | `graph_constructor.py`; safe: “prototype graph representation.” |
| D. Why minimizers? | They reduce sketch density versus all k-mers while keeping order. They trade resolution for computational compactness. | `minimizer_engine.py`; `O(Lw)` current implementation. |
| E. Why not all k-mers? | All k-mers create more tokens/edges and memory pressure; minimizers are a sampling heuristic. | Follow-up information loss; safe: “we quantify structural, not exact base accuracy.” |
| E. What does k do? | Larger k is more specific but can fragment/noise sensitivity; smaller k shares more and may be ambiguous. | adaptive `k` in `quality_feedback.py`. |
| E. What does window do? | It controls minimizer sampling density; larger windows generally select fewer representatives. | current selection loop; follow-up why no deque. |
| E. Is GPU implemented? | For minimizer extraction only: CUDA kernels with exact CPU parity on a laptop RTX 3070 Ti and a 2.8–12.4× minimizer-only speed-up. Graph construction is CPU, the final evaluation ran on CPU, and HPC GPUs are untested. | `minimizer_engine.py`, `outputs/gpu_validation/laptop/`; exact wording. |
| F. How are nodes created? | One node per distinct minimizer string after per-assembly support grouping. | `grouped_records`, `node_lookup`; support counts assemblies. |
| F. What if a minimizer repeats? | Records may repeat; support counts an assembly once, immediate path duplicates are compressed, later revisits remain. | graph path loop; follow-up cycles. |
| F. What is an edge? | A directed adjacent transition observed in one or more input minimizer paths. | `edge_supports`; safe: “observed path-derived edge.” |
| F. Why `0M` GFA overlap? | The prototype does not compute base overlaps; `0M` is a simplified link representation. | `write_proposed_gfa`; limitation. |
| G. Why preserve private nodes? | To avoid silently erasing an observed divergent assembly walk. They are marked rescued/private, not certified biological variants. | rescue metadata. |
| G. What about reverse complements? | They are not canonicalized/orientation-aware in current source. | forward substring/hash; future work. |
| H. Is consistency circular? | Partly: edge/path structure is derived from the same observed paths being scored. | formulas Section 6; safe: “internal structural coherence, not accuracy.” |
| H. What does .99 mean? | Approximately .99 on this bounded structural formula, not 99% biological accuracy. | `graph_consistency.py`; follow-up validation. |
| I. Why U/D/C anchors? | U favours low path coverage, D relative recurrence, C local structural coherence. Weighted scoring makes trade-offs explicit. | `anchor_selection.py`; formula Section 7. |
| I. Why `.35/.35/.30`? | Fixed design defaults; they are not learned. The weights sum to one and tests show each changes ranking. | constants/tests; future tune held-out. |
| I. Why should I trust anchor score? | It is explainable and deterministic, not independently biologically validated. | report CSV; safe limitation. |
| I. What if no anchors? | Candidate summary handles empty lists; coverage becomes zero and acceptance fails. | `anchor_summary`, quality gates. |
| J. What is PG-SCUnK here? | A project-local bounded feedback controller, not a claimed canonical external implementation. | `quality_feedback.py`; safe phrase. |
| J. What prevents infinite loops? | `visited`, finite candidate generation, `max_rounds`, and stagnation stopping. | `canonical_configuration`, loop. |
| J. Why visited set? | It prevents duplicate expensive graph rebuilds and cycles caused by parameter clamping. | `visited`; follow-up canonical floats. |
| J. Can it optimize a bad objective? | Yes; it optimizes its Q, so objective validity must be tested independently. | Q formula; future external validation. |
| J. Best vs final candidate? | The final tried candidate may be worse; the returned graph is protected best-so-far. | `best_graph`, `final_graph`; stress `c004/c007`. |
| K. Is this global optimization? | No, it is deterministic local neighbourhood search. | candidate generator; safe: “bounded heuristic.” |
| L. Complexity? | Current minimizers are roughly `O(sum Lw)` and store `O(sum L)` records; graph/report costs depend on distinct nodes, edges and path visits. | loops in engine/constructor. |
| M. Why dataclasses? | They make stage contracts explicit and serializable with `asdict`. | `models.py`; follow-up mutability. |
| N. Why dictionaries/sets? | They aggregate nodes/edges and guarantee each assembly supports an item once. | `defaultdict`, sets. |
| O. What do tests prove? | Deterministic implementation behaviour and report invariants on selected cases. | test table; not biological performance. |
| P. Why preserve baseline? | It keeps established alignment-oriented outputs distinct for future fair comparison. | `baseline_proposed_comparison.md`. |
| Q. Why was anchor coverage removed from Q? | No aligner consumes the selected anchors, and because `U` favours private nodes, coverage mostly tracked the `max_anchors` cap (APP 12/36, PSEN1 11/124). It previously rejected APP and PSEN1 for that reason alone. | commit `22076bc`; `test_metric_definitions.py`. |
| Q. Why does the controller not retune on real data? | All three genes exceed the `.72` threshold on the initial configuration (Q ≈ .98), so it stops after one evaluation. | final summary; stress demo is the only retuning run. |
| R. Why not simply PanAligner? | PanAligner is retained for alignment; this prototype explores an alignment-free, path/quality-control layer. | baseline scripts; no replacement claim. |
| S. Is minigraph your algorithm? | No. It is an external baseline tool used in legacy orchestration. | `paper_evaluation.py`, `align.py`. |
| T. How was the GPU path checked? | Seeded random and benchmark parity against the CPU reference, CUDA-event phase timings, memory-pool bytes, and a run manifest with commit and environment. | `scripts/gpu_minimizer_bench.py`; one laptop GPU only. |
| U. What exactly is novel? | The project’s proposed combination: explicit minimizer paths, rescue, graph-derived scoring, consistency-aware anchors and bounded retuning. | Safe: “proposed integration/prototype,” not invention of primitives. |
| V. Biggest limitation? | Structural scores reuse observed construction paths and do not establish biological accuracy. | Section 6/13. |
| W. Next experiment? | Truth-backed external benchmark with same inputs/baselines, orientation-aware graph and runtime/memory comparison. | future-work table. |
| X. What would invalidate hypothesis? | If orthogonal benchmark shows no improvement/robustness, poor scaling, or score has no relation to external truth. | Honest scientific answer. |

## 18. Quick revision sheet

**30 seconds:** “I retained the PanAligner/minigraph baseline and added a separate CPU alignment-free prototype.  It turns assembly minimizers into explicit per-assembly graph paths, preserves low-support nodes as rescued, measures internal structural coherence, ranks anchors, and uses a bounded feedback controller to try alternative parameters.  The final 174-assembly run completed and validated, but its scores are not biological accuracy.  Minimizer extraction also has a CUDA path that matched the CPU exactly on a laptop GPU with a 3–12× minimizer-only speed-up; the graph itself is built on the CPU.”

**1 minute:** Add the three results: PSEN2 accepted with Q `.98851`/coverage `1.0`; APP Q `.88639` but coverage `.333`; PSEN1 Q `.84871`, coverage `.10484`, 7 configurations and no acceptance.  Explain acceptance requires Q plus coverage/node/path gates.

**3 minutes:** Explain minimizers -> node/edge/path construction -> `C(v), C(e), C(p), C(G)` -> anchors `U,D,C` -> Q -> controller.  State rescue and circularity.  Close with fair future benchmark/GPU/orientation work.

**Five formulas:** node `(.25,.25,.25,.25)`; edge `(.5,.5)`; path `(.4,.3,.3)`; graph `(.4,.3,.3)`; anchor `.35U+.35D+.30C`; Q weights `.25,.20,.15,.15,.15,.10`.

**Must know parameters:** default final initial `k=11`, window `6`, support `2`, anchors `12`, spacing `4.0`, weights `.35/.35/.30`, threshold `.72` (no anchor-coverage gate since `22076bc`), searched parameters `k` and `window` only, tolerance `1e-4`, max rounds `3`.

**20 must-know questions:** What is a minimizer? Why paths? What is rescued? What exactly does the GPU do, and what not? How is node support counted? What is an observed edge? Is consistency circular? Q vs acceptance? What is U/D/C? Why fixed weights? What does visited do? Why best-so-far? Why was anchor coverage removed from Q? Why does the controller not retune on real data? Why is min_node_support not searched? What validates JSON? What do tests not prove? Is this PanAligner? Is minigraph proposed? What future evidence is required? Where do Model A's labels come from? Why is MAPQ a weak baseline here? Why is Model B's result negative, and why is it not an Alzheimer's classifier? Why was the legacy ML path retired?

## 19. CLAIMS I MUST NOT MAKE

| Unsafe claim | Safe, accurate replacement |
|---|---|
| “Our GPU implementation accelerates the pipeline.” | “On one laptop GPU, minimizer extraction is 2.8–12.4× faster with exact CPU parity; graph construction and later stages are CPU-only and HPC GPUs are untested.” |
| “A .99 consistency score means 99% accuracy.” | “It is a bounded graph-derived structural coherence score, not an accuracy percentage.” |
| “We reproduced PanAligner.” | “We retained PanAligner/minigraph baseline components and implemented a separate proposed alignment-free prototype.” |
| “Our graph is biologically correct because paths validate.” | “Path validation checks internal graph/path integrity; biological correctness needs independent evidence.” |
| “Private nodes are real variants.” | “Private low-support nodes are retained observed structure and may include real divergence or artefacts.” |
| “PG-SCUnK learns optimal parameters.” | “The project-local controller deterministically explores a small heuristic neighbourhood.” |
| “All three genes are accepted, so the graphs are high quality.” | “All three pass a structural-coherence threshold on the first configuration; the score is high partly by construction and is not biological accuracy.” |
| “The feedback loop improved the real-data graphs.” | “On APP/PSEN1/PSEN2 the controller accepted the initial configuration; retuning was only demonstrated on a synthetic stress demo.” |
| “PSEN1 is biologically worse.” | “The per-gene differences are small structural-score differences; no biological conclusion follows.” |
| “Our model replaces MAPQ.” | “On simulated locus-scale reads it beats MAPQ when few alignments are discarded, but MAPQ 60 is about equal at its own operating point.” |
| “Our model detects misaligned reads with 0.98 PR-AUC.” | “On simulated locus-scale reads, mostly off-target alignments; within-locus misplacements are too few (96) to judge, and whole-genome alignment is untested.” |
| “Model B classifies disease-associated variants.” | “Model B tested whether sequence-derived features recover ClinVar classification direction in three genes. It did not beat simple baselines (0.734 vs 0.848 on missense), and it is a negative research finding, not a usable classifier.” |
| “GC content is a biological signal of pathogenicity.” | “GC appears to carry regional or gene-associated signal here; that doesn’t prove sequence lacks biological information.” |
| “Model B generalises to held-out genes.” | “Three genes cannot show generalisation, and existing test folds have been read, so further runs on them are exploratory.” |
| “The ML predicts Alzheimer’s risk.” | “Nothing in this project is a clinical or diagnostic prediction.” |
| “The proposed method beats PanAligner.” | “No direct biological performance comparison is established by these non-interchangeable metrics.” |
| “Forward GFA means strand handling is solved.” | “The current GFA emits forward `+` paths/links; orientation-aware modelling is future work.” |

## 20. Review checklist before the viva

- Open `outputs/proposed_final_evaluation/proposed_final_summary.csv` and be able to state the three rows.
- Be ready to navigate from `proposed_methodology.py` to `final_evaluation.py`, then `quality_feedback.py` and `graph_constructor.py`.
- Say “structural coherence” rather than “accuracy” unless discussing actual baseline alignment metrics separately.
- State the CPU/GPU limitation before an examiner has to find it.
- If unsure, say what the source proves, what it does not prove, and the experiment needed to establish the missing claim.
