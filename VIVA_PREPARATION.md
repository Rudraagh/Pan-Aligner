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

Implemented: CPU minimizers, minimizer-node graph building, explicit paths, private-node rescue, graph-derived scores, deterministic anchor selection, a bounded local feedback heuristic, artifacts, unit tests, and a completed three-gene proposed evaluation.  Simplified/prototype: all of those algorithms are small project-local implementations, not a full biological graph assembler or a reproduction of PanAligner internals.  Proposed/future: genuine GPU kernels, richer orientation/base-level graph representation, and biological benchmark evidence.  Not implemented: GPU-accelerated minimizer computation in the current kernel, exact PanAligner algorithm reproduction, and proof of biological superiority.

## 2. Novelty analysis and safe claims

### Contribution A — GPU-parallel minimizer concept

**Idea and difference.** A conventional CPU sketcher computes all windows serially; the architecture exposes a `backend` choice and detects CuPy.  **Code:** `scripts/proposed/minimizer_engine.py:compute_minimizers`.  **What is actually implemented:** `_compute_cpu_minimizers` always performs the work.  Even if CuPy imports, the code only reports `backend_used = "gpu"` and says it uses the shared deterministic minimizer logic; it does not transfer arrays or launch a GPU kernel.  The final evaluator explicitly calls `backend="cpu"`.

**Safe viva claim:** “GPU acceleration is a stated architecture/future contribution; the validated implementation is deterministic CPU minimizer extraction with GPU readiness/detection only.”  Do not claim speed-up, CUDA execution, GPU parity, or parallel kernels.  Examiner attack: *“Where is the GPU kernel?”* Strong answer: “There is none in the current implementation; CuPy detection is not computational acceleration.  A real version would batch encoded bases/k-mers in device arrays, compute hashes/window minima on device, compact duplicates, transfer records, and validate bit-for-bit or semantically against CPU.”

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
| `minimizer_engine.py` | SHA-1-prefix stable hashing and CPU minimizer extraction. | GPU caveat and duplicate suppression. |
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

Complexity is `O(V + total_path_visits + V log V + VA)` where `V` is node/candidate count and `A` is accepted anchors (spacing scan); with small anchor budgets, the last factor is bounded.  `U` is not biological specificity, `D` depends on this dataset/path representation, and `C` inherits the structural circularity.  Anchor coverage is the fraction of input assembly IDs appearing in any selected anchor’s `haplotypes`, so it can be low even when graph consistency is high.

## 8. PG-SCUnK-style controller

`run_pg_scunk_feedback_loop()` starts from `_adaptive_default_parameters`: `k` is clamped to 3..11 from half the shortest sequence; window is clamped to 2..6; minimum support is 2 for at least three assemblies; anchor count is clamped to 4..12 based on assembly count; spacing is 4.0 for shortest sequence >=40 else 2.0.  Overrides may replace these values.

`canonical_configuration()` makes an immutable ordered tuple over `k, window_size, min_node_support, max_anchors, min_anchor_spacing, anchor_alpha, anchor_beta, anchor_gamma`; float values are rounded to 12 decimals.  `visited` prevents recomputing the same configuration.  Each candidate actually calls `build_proposed_pangenome_graph`, so it rebuilds the graph rather than re-labeling an old score.

Diagnostic signals are: anchor coverage, node consistency, node support ratio, edge consistency, path consistency, and `1 - min(1, branching_node_ratio)`.  The smallest signal (with fixed priority tie-breaking) chooses an action.  Anchor weakness proposes lower spacing, larger anchor budget, or a normalized alpha boost.  Path weakness proposes lower `k` and/or window.  Other weaknesses propose lower `k`, lower window, or lower support.  This is a small deterministic neighbourhood—not Bayesian optimisation, gradient descent, or exhaustive search.

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

Acceptance is stricter than the overall threshold: `_quality_passes` requires overall quality >= `.72`, anchor coverage >= `.70`, node consistency >= `.45`, and path consistency >= `.45` by default.  Thus a high quality score alone does not mean accepted.  Stopping reasons are `quality_threshold_reached`, `max_rounds_reached`, `no_unseen_configurations`, or `stagnation_patience_exceeded`; default maximum rounds is 3 and default stagnation patience is 2.

Best-so-far protection matters because the final evaluated candidate can be worse.  `best_graph` is returned separately from `final_graph`; final evaluation serializes the best graph.  A duplicate candidate is skipped; an immediate threshold success ends after its first accepted evaluation; worse candidates remain in the audit trail but do not replace best.

### Stress-demo trace

`outputs/proposed_step5_pgscunk_demo/experiment_metadata.json` labels the four short controlled-substitution haplotypes **DEMONSTRATION ONLY — not a biological benchmark**.  It began at `Q=0.67796875`, target `.99`, `k=6, window=5, support=2`.  Seven unique configurations were rebuilt in three rounds.  The best was candidate `c006`, round 3: `k=4, window=4, support=2`, `Q=0.80341146`; the final candidate `c007` was worse (`Q=.76959245`).  It stopped at `max_rounds_reached`, not success.  This demonstrates rebuilding, diagnostics and best-so-far retention; it does not validate a gene or disease claim.

## 9. Exact quality function and validation

The source in `graph_constructor.py` calculates:

```text
Q = .25 mean_node_consistency + .20 mean_edge_consistency
  + .15 mean_path_consistency + .15 mean_node_support_ratio
  + .15 anchor_coverage_ratio
  + .10 (1 - min(1, branching_node_ratio))
```

`mean_node_support_ratio` averages `node.support / assembly_count`; `branching_node_ratio` is nodes with more than one outgoing edge divided by nodes.  This differs from graph consistency `C(G)`, which uses only node/edge/path consistency `.40/.30/.30`.  Anchor coverage is separate.  `validation.json` is a software/data-integrity check: paths are valid, score fields are in range, PG-SCUnK keys are unique, selected anchor positions obey spacing, best tracking is consistent, and CSV/JSON report counts/quality agree.  Validation is not biological correctness, statistical significance, or benchmark accuracy.

## 10. Final real-data evaluation

Source: `outputs/proposed_final_evaluation/proposed_final_summary.json`; each per-gene `validation.json` passed.  The final evaluator uses all listed train-manifest inputs and forces CPU.

| Gene | Assemblies | Nodes | Edges | rescued/private | anchors | node / edge / path C | graph C | anchor cov. | Q | rounds/configs | accepted | backend |
|---|---:|---:|---:|---:|---:|---|---:|---:|---:|---|---|---|
| APP | 36 | 72,646 | 79,053 | 116 / 116 | 12 | .98248 / .99867 / .98755 | .98886 | .33333 | .88639 | 2 / 4 | no | cpu |
| PSEN1 | 124 | 21,397 | 23,118 | 417 / 417 | 16 | .98380 / .98488 / .98879 | .98562 | .10484 | .84871 | 3 / 7 | no | cpu |
| PSEN2 | 14 | 15,635 | 16,089 | 46 / 46 | 12 | .98929 / .98891 / .99217 | .99004 | 1.00000 | .98851 | 1 / 1 | yes | cpu |

PSEN2 is the strongest result by Q, graph consistency, path consistency and acceptance: its anchor coverage meets `.70`.  PSEN1 is weakest by Q because its selected anchors cover only about 10.5% of paths, despite high structural consistency; it increased its anchor limit from 12 to 16 in the best configuration, yet did not satisfy coverage and stopped at the 3-round limit.  APP is notable/surprising: very high graph consistency and nearly perfect edge consistency coexist with only one-third anchor coverage, so it is not accepted.  This is exactly why “Q is high” must not be conflated with controller acceptance.

The run covered 174 assemblies.  It validates serialisation and structural properties, not a held-out biological ground truth.  The baseline alignment metrics in `outputs/evaluation/` are intentionally separate and not numerically interchangeable with this table.

## 11. Testing evidence

| Test file | Coverage / important cases | It catches | It does not prove |
|---|---|---|---|
| `test_proposed_minimizer_engine.py` | deterministic records, parameter validation, short inputs, backend metadata. | unstable/minimum/empty-input regressions. | GPU acceleration or sequencing correctness. |
| `test_proposed_graph_constructor.py` | node/edge/path construction, ordered paths, private-node rescue, validation and reports. | missing path edges or accidental low-support deletion. | biological graph accuracy. |
| `test_proposed_graph_consistency.py` | node/edge/path/graph formulas and reports. | formula/report regressions. | independence from construction circularity. |
| `test_proposed_anchor_selection.py` | weight validation/formula, ranking effects, deterministic spacing, reports. | incorrect score/ranking/spacings. | biological anchor specificity. |
| `test_proposed_haplotype_scoring.py` | legacy haplotype scorer behaviour. | legacy scorer regressions. | proposed graph construction quality. |
| `test_proposed_quality_feedback.py` | candidate generation, visited keys, improvement/best retention, stopping conditions. | duplicate/replacement/control-flow errors. | optimal tuning. |
| `test_proposed_pg_scunk_step4.py` | Step-4 loop/report semantics. | broken controller integration. | external PG-SCUnK equivalence. |
| `test_proposed_pg_scunk_demo.py` | controlled stress demo/smoke evidence. | demo metadata/output regressions. | real-world benchmark performance. |

Tests demonstrate deterministic software properties on selected synthetic/small cases.  They do not constitute a full published benchmark, biological validation, memory scalability study, or CPU/GPU parity study.

## 12. Baseline, PanAligner, minigraph, and legacy ML

The retained baseline uses project scripts such as `build_graph.py`, `align.py`, `evaluate.py`, and `paper_evaluation.py` to orchestrate external minigraph/PanAligner binaries and parse GAF/alignment-oriented outputs.  `paper_evaluation.py` contains local minigraph-first logic and falls back/selects PanAligner under its stated conditions.  Feature extraction and ML (`extract_features.py`, `train_ml.py`, `predict_ml.py`) consume PanAligner-derived alignment features for downstream classification; they are not part of proposed `Q`.

The proposed constructor is not “just minigraph”; it is a separate Python minimizer-walk prototype.  Conversely, it is not a replacement claim for PanAligner, and it is not an exact reproduction of PanAligner internal algorithms.  Baseline files are preserved so a fair future comparison can use identical datasets, splits, hardware, time/memory measurement, fixed tool versions, stated parameters, and appropriate comparable outcomes (for example truth-backed variant/path/alignment measures), rather than comparing structural Q directly to alignment accuracy.

## 13. Limitations and future work

| Limitation | Impact/current-review wording | Realistic future work |
|---|---|---|
| CPU routine only | No measured GPU speedup; final backend is CPU. | Device kernels, batched encoding/reductions, duplicate compaction and parity tests. |
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

Future work must be described as future work, not quietly folded into the implementation claim: (1) a genuine GPU minimizer kernel with CPU/GPU parity and timing; (2) canonical, orientation-aware paths and base/segment-level graph overlaps; (3) independent biological/path/variant validation; (4) adaptive anchor budgets and weights selected on held-out data; (5) stronger multi-objective parameter search; (6) external benchmark datasets and repeated runtime/memory experiments; and (7) a fair, versioned comparison against minigraph/PanAligner.  The expected payoff is stronger biological and systems evidence; the required discipline is to measure it rather than infer it from the current structural score.

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
| `try: import cupy` | Capability detection, not a computation kernel. |
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
| `backend="cpu"` final call | Final evidence explicitly avoids claiming GPU execution. |

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
| E. Is GPU implemented? | No; execution is CPU. CuPy is detected only and final evaluation explicitly requests CPU. | `minimizer_engine.py`, `final_evaluation.py`; safe exact wording. |
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
| J. Best vs final candidate? | The final tried candidate may be worse; the returned graph is protected best-so-far. | `best_graph`, `final_graph`; stress `c006/c007`. |
| K. Is this global optimization? | No, it is deterministic local neighbourhood search. | candidate generator; safe: “bounded heuristic.” |
| L. Complexity? | Current minimizers are roughly `O(sum Lw)` and store `O(sum L)` records; graph/report costs depend on distinct nodes, edges and path visits. | loops in engine/constructor. |
| M. Why dataclasses? | They make stage contracts explicit and serializable with `asdict`. | `models.py`; follow-up mutability. |
| N. Why dictionaries/sets? | They aggregate nodes/edges and guarantee each assembly supports an item once. | `defaultdict`, sets. |
| O. What do tests prove? | Deterministic implementation behaviour and report invariants on selected cases. | test table; not biological performance. |
| P. Why preserve baseline? | It keeps established alignment-oriented outputs distinct for future fair comparison. | `baseline_proposed_comparison.md`. |
| Q. Why APP not accepted? | Anchor coverage is .333, below `.70`, despite Q .886 and high structural consistency. | APP summary; exact answer. |
| Q. Why PSEN1 only ~10% anchors? | The selected 16 anchors cover only .10484 of assembly IDs under current score/spacing; controller ran 7 configs but hit max rounds. | PSEN1 summary; no speculative biology. |
| R. Why not simply PanAligner? | PanAligner is retained for alignment; this prototype explores an alignment-free, path/quality-control layer. | baseline scripts; no replacement claim. |
| S. Is minigraph your algorithm? | No. It is an external baseline tool used in legacy orchestration. | `paper_evaluation.py`, `align.py`. |
| T. How would GPU be done correctly? | Parallel device k-mer encoding/hashing/window minima and compaction, plus CPU/GPU parity, timing and memory experiments. | future work; safe no current claim. |
| U. What exactly is novel? | The project’s proposed combination: explicit minimizer paths, rescue, graph-derived scoring, consistency-aware anchors and bounded retuning. | Safe: “proposed integration/prototype,” not invention of primitives. |
| V. Biggest limitation? | Structural scores reuse observed construction paths and do not establish biological accuracy. | Section 6/13. |
| W. Next experiment? | Truth-backed external benchmark with same inputs/baselines, orientation-aware graph and runtime/memory comparison. | future-work table. |
| X. What would invalidate hypothesis? | If orthogonal benchmark shows no improvement/robustness, poor scaling, or score has no relation to external truth. | Honest scientific answer. |

## 18. Quick revision sheet

**30 seconds:** “I retained the PanAligner/minigraph baseline and added a separate CPU alignment-free prototype.  It turns assembly minimizers into explicit per-assembly graph paths, preserves low-support nodes as rescued, measures internal structural coherence, ranks anchors, and uses a bounded feedback controller to try alternative parameters.  The final 174-assembly run completed and validated, but its scores are not biological accuracy and GPU execution is future work.”

**1 minute:** Add the three results: PSEN2 accepted with Q `.98851`/coverage `1.0`; APP Q `.88639` but coverage `.333`; PSEN1 Q `.84871`, coverage `.10484`, 7 configurations and no acceptance.  Explain acceptance requires Q plus coverage/node/path gates.

**3 minutes:** Explain minimizers -> node/edge/path construction -> `C(v), C(e), C(p), C(G)` -> anchors `U,D,C` -> Q -> controller.  State rescue and circularity.  Close with fair future benchmark/GPU/orientation work.

**Five formulas:** node `(.25,.25,.25,.25)`; edge `(.5,.5)`; path `(.4,.3,.3)`; graph `(.4,.3,.3)`; anchor `.35U+.35D+.30C`; Q weights `.25,.20,.15,.15,.15,.10`.

**Must know parameters:** default final initial `k=11`, window `6`, support `2`, anchors `12`, spacing `4.0`, weights `.35/.35/.30`, threshold `.72`, anchor-coverage gate `.70`, tolerance `1e-4`, max rounds `3`.

**20 must-know questions:** What is a minimizer? Why paths? What is rescued? Why not GPU claim? How is node support counted? What is an observed edge? Is consistency circular? Q vs acceptance? What is U/D/C? Why fixed weights? What does visited do? Why best-so-far? Why APP fails? Why PSEN1 fails? Why PSEN2 passes? What validates JSON? What do tests not prove? Is this PanAligner? Is minigraph proposed? What future evidence is required?

## 19. CLAIMS I MUST NOT MAKE

| Unsafe claim | Safe, accurate replacement |
|---|---|
| “Our GPU implementation accelerates minimizers.” | “GPU is an architectural/future goal; the validated minimizer computation runs on CPU.” |
| “A .99 consistency score means 99% accuracy.” | “It is a bounded graph-derived structural coherence score, not an accuracy percentage.” |
| “We reproduced PanAligner.” | “We retained PanAligner/minigraph baseline components and implemented a separate proposed alignment-free prototype.” |
| “Our graph is biologically correct because paths validate.” | “Path validation checks internal graph/path integrity; biological correctness needs independent evidence.” |
| “Private nodes are real variants.” | “Private low-support nodes are retained observed structure and may include real divergence or artefacts.” |
| “PG-SCUnK learns optimal parameters.” | “The project-local controller deterministically explores a small heuristic neighbourhood.” |
| “APP succeeded because Q is high.” | “APP’s Q is high but it was not accepted because anchor coverage was below the required gate.” |
| “PSEN1 is biologically worse.” | “Under this anchor-selection metric PSEN1 had low selected-anchor coverage; no biological conclusion follows.” |
| “The proposed method beats PanAligner.” | “No direct biological performance comparison is established by these non-interchangeable metrics.” |
| “Forward GFA means strand handling is solved.” | “The current GFA emits forward `+` paths/links; orientation-aware modelling is future work.” |

## 20. Review checklist before the viva

- Open `outputs/proposed_final_evaluation/proposed_final_summary.csv` and be able to state the three rows.
- Be ready to navigate from `proposed_methodology.py` to `final_evaluation.py`, then `quality_feedback.py` and `graph_constructor.py`.
- Say “structural coherence” rather than “accuracy” unless discussing actual baseline alignment metrics separately.
- State the CPU/GPU limitation before an examiner has to find it.
- If unsure, say what the source proves, what it does not prove, and the experiment needed to establish the missing claim.
