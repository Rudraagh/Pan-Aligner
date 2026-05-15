# 10. Viva Questions

This file is designed for confident viva preparation.

For each question, you get:

- short answer
- detailed answer
- follow-up questions
- easy way to remember

There are more than `100` questions here.

## Project Basics

### 1. What is the title of your project?

Short answer: `Implementation and Reproduction of PanAligner for Co-linear Chaining on Cyclic Pangenome Graphs.`

Detailed answer: The project reproduces the PanAligner paper by combining real pangenome graph construction and real graph alignment with a separate educational theory layer that explains SCCs, DAG conversion, path cover, precedence, and simplified co-linear chaining.

Follow-up questions: Why did you choose this title? What parts are implementation and what parts are reproduction?

Easy way to remember: `real tool + theory explanation`

### 2. What is the main objective of your project?

Short answer: To understand, implement, and explain the PanAligner workflow and paper concepts clearly.

Detailed answer: The objective is to build pangenome graphs, align held-out sequences to them using the real PanAligner, and separately create simple demos for the paper's difficult graph algorithms so the project is both practical and teachable.

Follow-up questions: Is your goal prediction or alignment? Why is explanation important?

Easy way to remember: `run it, understand it, explain it`

### 3. What paper is your project based on?

Short answer: `Co-linear chaining on pangenome graphs.`

Detailed answer: The project is based on the paper about co-linear chaining on pangenome graphs, especially the challenge of handling cyclic graph structure while finding a valid chain of anchors.

Follow-up questions: What is the central idea of the paper? Why is chaining difficult on graphs?

Easy way to remember: `paper about chaining on graph paths`

### 4. What does your project do in one sentence?

Short answer: It builds pangenome graphs and aligns held-out sequences to them using PanAligner, while also teaching the paper's algorithms with simplified demos.

Detailed answer: The project preprocesses genomic FASTA data, creates train/test splits, builds train graphs with minigraph, aligns test sequences using PanAligner, evaluates the alignments, and separately reproduces the core theory in a beginner-friendly Python layer.

Follow-up questions: Why two layers? What is the benefit of held-out evaluation?

Easy way to remember: `graph build + graph align + graph explain`

### 5. Is your project a machine learning project?

Short answer: No, not in its current main form.

Detailed answer: Although some older ML-related scripts remain in the repository from earlier exploration, the current project is framed as a PanAligner paper-reproduction and graph-alignment project. The authoritative flow is `main.py`, which focuses on preprocessing, graph construction, theory demos, PanAligner alignment, and evaluation.

Follow-up questions: Why are ML files still present? Which file is the authoritative entry point?

Easy way to remember: `current scope is graphs, not ML`

### 6. What is the difference between implementation and reproduction here?

Short answer: Implementation means running the practical pipeline, and reproduction means recreating the paper concepts and workflow.

Detailed answer: The implementation side includes preprocessing, graph building, PanAligner integration, and held-out evaluation. The reproduction side includes simplified educational demos of SCCs, DAG conversion, path cover, precedence, and chaining so the paper can be understood and defended clearly.

Follow-up questions: Did you reproduce exact internals? Did you use the official tool?

Easy way to remember: `implementation runs, reproduction explains`

### 7. Which file is the main entry point?

Short answer: `main.py`

Detailed answer: `main.py` is the single master driver of the project. It controls preprocessing, splitting, graph building, theory execution, evaluation, and report generation through clearly defined command-line modes.

Follow-up questions: What are the modes? Why use a single entry point?

Easy way to remember: `main.py controls everything`

### 8. What are the three main modes in `main.py`?

Short answer: `--full-pipeline`, `--theory-only`, and `--evaluate`.

Detailed answer: Full pipeline runs almost everything end to end, theory-only runs just the educational graph-algorithm demos, and evaluate runs held-out PanAligner alignment evaluation using the current graph setup.

Follow-up questions: Which mode is best for viva demonstration? Which mode uses PanAligner?

Easy way to remember: `everything, theory, evaluation`

### 9. Which genes are used in your project?

Short answer: `APP`, `PSEN1`, and `PSEN2`.

Detailed answer: The repository processes three genes commonly associated with Alzheimer's-related studies in this dataset: APP, PSEN1, and PSEN2. Each gene has a reference sequence and healthy and unhealthy sample sequences.

Follow-up questions: Are the graphs built per gene? Why three genes?

Easy way to remember: `APP-PSEN1-PSEN2`

### 10. How many reference sequences are there?

Short answer: `3`

Detailed answer: There is exactly one reference sequence per gene, so one for APP, one for PSEN1, and one for PSEN2. The preprocessing script enforces that each gene must have exactly one reference sequence.

Follow-up questions: Why exactly one reference? What happens if more than one is found?

Easy way to remember: `one reference per gene`

## Pangenome Graph Basics

### 11. What is a pangenome graph?

Short answer: A graph representation of many related genome sequences together.

Detailed answer: Instead of storing only one linear reference, a pangenome graph represents shared regions and variant branches using nodes and edges, so multiple possible sequence paths can coexist in the same structure.

Follow-up questions: Why is that better than one reference? What do nodes represent?

Easy way to remember: `many genomes, one graph`

### 12. Why use a graph instead of a linear genome reference?

Short answer: Because a graph can represent variation more naturally.

Detailed answer: A linear reference forces every query to match one path, but real biological sequences contain variation. A graph can store multiple alternative routes, which makes alignment more flexible and better suited to diverse populations or variant-rich regions.

Follow-up questions: What kind of variation? How does this help alignment?

Easy way to remember: `graph handles alternate paths`

### 13. What do graph nodes represent in your project?

Short answer: Sequence segments.

Detailed answer: In the GFA-based graphs used here, each node corresponds to a segment of sequence. Different sequences can share nodes if they share sequence content, and branch where variation occurs.

Follow-up questions: Do nodes store full genomes? How are edges used?

Easy way to remember: `node = sequence piece`

### 14. What do graph edges represent?

Short answer: Valid transitions between sequence segments.

Detailed answer: Edges represent adjacency or allowed continuation from one graph segment to another. If a query traverses a path through the graph, it follows these edges between nodes.

Follow-up questions: Are edges directed? What happens at branches?

Easy way to remember: `edge = can go next`

### 15. What is the practical meaning of a path in a pangenome graph?

Short answer: One possible sequence through the graph.

Detailed answer: A path is a valid ordered traversal through nodes and edges. It can represent the reference sequence, a sample-specific variant path, or some shared combination of graph segments.

Follow-up questions: Can multiple samples share part of a path? Why is that useful?

Easy way to remember: `path = one possible genome route`

### 16. What tool builds the graph in your project?

Short answer: `minigraph`

Detailed answer: The project uses the real `minigraph` tool to construct pangenome graphs from a reference FASTA and multiple sample FASTA sequences. This produces `.gfa` graph files.

Follow-up questions: What does PanAligner do then? Why not use PanAligner for construction?

Easy way to remember: `minigraph makes graphs`

### 17. What graph types are built for each gene?

Short answer: Healthy, unhealthy, and combined graphs.

Detailed answer: For each gene, the build step creates three graph views: one graph using healthy samples, one using unhealthy samples, and one combined graph using both. The combined graph is especially useful for held-out evaluation.

Follow-up questions: Why build all three? Which one is used for evaluation?

Easy way to remember: `H, U, and H+U`

### 18. Why is the combined graph important?

Short answer: It captures the broadest variation.

Detailed answer: The combined graph contains information from both healthy and unhealthy samples, so it provides a more inclusive target when aligning held-out sequences.

Follow-up questions: Could healthy-only miss some variants? Why not evaluate on all three?

Easy way to remember: `combined = most complete`

### 19. What is the current APP graph size in your repo?

Short answer: `115 nodes and 154 edges` for one APP combined graph stats file.

Detailed answer: The graph statistics file `outputs/graphs/app.combined.stats.json` shows the APP combined graph has 115 nodes and 154 edges, one weakly connected component, and no cycles detected in that specific graph snapshot.

Follow-up questions: Does no cycle mean the paper is irrelevant? Why still teach cyclic handling?

Easy way to remember: `APP 115/154`

### 20. Are all real graphs in your repo strongly cyclic?

Short answer: No.

Detailed answer: Many of the current locus graphs are mostly DAG-like. That is why the theory layer creates an educational cyclic copy when needed, so the cyclic case from the paper can still be explained clearly.

Follow-up questions: Did you alter the real graph? Why synthetic back edges?

Easy way to remember: `real graphs mostly DAG-like, theory adds cycles`

## FASTA, Preprocessing, and Split

### 21. What is FASTA?

Short answer: A common text format for storing biological sequences.

Detailed answer: FASTA stores sequence entries with a header line starting with `>` followed by sequence lines. In this project, the root combined FASTA files contain reference and labeled sequences for the three genes.

Follow-up questions: What information is in the headers? Why preprocess them?

Easy way to remember: `header + sequence`

### 22. Why do you preprocess the FASTA files?

Short answer: To clean, validate, and organize them for downstream use.

Detailed answer: Preprocessing separates reference, healthy, and unhealthy sequences, validates DNA, normalizes sequence format, writes one-sample FASTA files, writes combined class FASTA files, and saves metadata for later stages.

Follow-up questions: What if a sequence is invalid? What if duplicates exist?

Easy way to remember: `clean and organize before graph work`

### 23. What does `preprocess.py` validate?

Short answer: DNA content and expected structure.

Detailed answer: It normalizes sequences, validates DNA characters, checks the parsed metadata from headers, enforces one reference per gene, and can deduplicate repeated sequences by gene and label.

Follow-up questions: Why one reference? What happens to invalid sequences?

Easy way to remember: `validate, separate, save`

### 24. What does deduplication mean in your preprocessing?

Short answer: Removing repeated identical sequences within the same gene and label group.

Detailed answer: The script tracks seen sequences using a key based on gene and label. If the same sequence appears again for the same group and deduplication is enabled, it is skipped.

Follow-up questions: Why deduplicate? Could duplicates bias graph construction?

Easy way to remember: `same gene + same label + same sequence = skip`

### 25. How many healthy and unhealthy sequences were processed?

Short answer: `110 healthy` and `110 unhealthy`.

Detailed answer: According to `preprocess_manifest.json`, the current project processed 110 healthy sequences, 110 unhealthy sequences, and 3 reference sequences, with no duplicates and no invalid sequences skipped in that snapshot.

Follow-up questions: How many genes? Were any sequences rejected?

Easy way to remember: `110 H, 110 U, 3 ref`

### 26. Why do you perform a train/test split?

Short answer: To test graph alignment on unseen sequences.

Detailed answer: If the same sequences used to build the graph are also used for evaluation, the test is less meaningful. The train/test split creates held-out sequences that act as unseen queries for PanAligner.

Follow-up questions: Why is this scientifically better? Is it classification?

Easy way to remember: `build on train, test on unseen`

### 27. Is the train/test split random?

Short answer: Yes, but reproducibly random.

Detailed answer: The split is created using a random seed, which is `42` by default, so the same split can be reproduced consistently across runs.

Follow-up questions: Why reproducibility? What is the test fraction?

Easy way to remember: `random but repeatable`

### 28. What is the default test fraction?

Short answer: `0.20`

Detailed answer: By default, 20 percent of sequences from each class are held out for testing. The rest are used for training and graph construction.

Follow-up questions: How is class balance maintained? What if a class is too small?

Easy way to remember: `80 train / 20 test`

### 29. How does the split preserve class structure?

Short answer: Healthy and unhealthy records are split separately for each gene.

Detailed answer: The script iterates through healthy and unhealthy buckets independently, shuffles indices with the seeded RNG, then assigns train and test records within each label category per gene.

Follow-up questions: Why not split all records together? Does this keep balance?

Easy way to remember: `split each class separately`

### 30. What is the overall current split count?

Short answer: `174 train sequences and 46 test sequences` excluding references.

Detailed answer: The split summary shows 87 healthy train, 87 unhealthy train, 23 healthy test, and 23 unhealthy test. That totals 174 train and 46 test sample sequences, plus 3 references.

Follow-up questions: How many queries were evaluated? Why 46?

Easy way to remember: `174 build, 46 test`

## GFA and GAF

### 31. What is GFA?

Short answer: A graph file format.

Detailed answer: GFA stands for Graphical Fragment Assembly and stores sequence graph information such as segments and links. In this project, minigraph outputs the pangenome graphs in GFA format.

Follow-up questions: What do `S` and `L` lines mean? How is it parsed?

Easy way to remember: `GFA = graph`

### 32. What is GAF?

Short answer: A graph alignment file format.

Detailed answer: GAF stands for Graph Alignment Format. It stores how a query sequence aligns to a graph path, including query span, path span, matches, block length, mapping quality, and optional tags.

Follow-up questions: Which tool outputs it? What metrics do you derive from it?

Easy way to remember: `GAF = alignment`

### 33. Which tool outputs GFA in your project?

Short answer: `minigraph`

Detailed answer: The graph-building stage calls minigraph on the reference and sample FASTA files and writes the resulting graph in GFA format.

Follow-up questions: What about PanAligner? Does it also output GFA?

Easy way to remember: `minigraph -> GFA`

### 34. Which tool outputs GAF in your project?

Short answer: `PanAligner`

Detailed answer: PanAligner aligns query FASTA sequences to a graph GFA and prints its alignments as GAF records, which the project saves to output files.

Follow-up questions: How do you parse the best alignment? Where are GAF files stored?

Easy way to remember: `PanAligner -> GAF`

### 35. What are `S` lines in GFA?

Short answer: Segment or node lines.

Detailed answer: In GFA, `S` lines define graph segments. In the project parser, they become nodes in a directed graph and can store sequence content or at least length.

Follow-up questions: What about `L` lines? What attributes are stored on nodes?

Easy way to remember: `S = segment`

### 36. What are `L` lines in GFA?

Short answer: Link or edge lines.

Detailed answer: `L` lines define directed links between segments. The project parser converts them into directed graph edges, optionally storing orientation metadata in the visualization layer.

Follow-up questions: Why directed? Does orientation matter?

Easy way to remember: `L = link`

### 37. What is coverage in your GAF parser?

Short answer: Aligned block length divided by query length.

Detailed answer: Coverage measures how much of the query was aligned. If nearly the whole query aligns, coverage approaches 1.

Follow-up questions: Can coverage be close to 1 with low identity? What does that mean?

Easy way to remember: `how much query got covered`

### 38. What is identity in your GAF parser?

Short answer: Residue matches divided by alignment block length.

Detailed answer: Identity measures how similar the aligned region is by comparing the number of matching residues to the total aligned block length.

Follow-up questions: Why is identity different from coverage? Which one measures similarity?

Easy way to remember: `how correct the aligned part is`

### 39. What is normalized score?

Short answer: Alignment score divided by query length.

Detailed answer: Since longer queries can naturally have larger raw scores, normalized score helps compare alignments more fairly across different query lengths.

Follow-up questions: Why not use raw score alone? How do you choose best alignment?

Easy way to remember: `score per query length`

### 40. How do you choose the best alignment from multiple GAF records?

Short answer: By normalized score first, then identity, coverage, and mapping quality.

Detailed answer: The parser's `best_alignment()` function ranks alignments primarily by normalized score, and breaks ties using identity, coverage, and MAPQ.

Follow-up questions: Why this order? Why not just use MAPQ?

Easy way to remember: `score first, quality next`

## Cycles, SCC, and DAG

### 41. What is a cycle in a graph?

Short answer: A path that returns to the starting node.

Detailed answer: In a directed graph, a cycle means you can follow edge directions and eventually come back to a previously visited node, often the starting node.

Follow-up questions: Why do cycles make analysis harder? Can a DAG have cycles?

Easy way to remember: `loop in directed path`

### 42. What is an SCC?

Short answer: A strongly connected component.

Detailed answer: An SCC is a maximal set of nodes in which every node can reach every other node through directed paths.

Follow-up questions: Can one node be an SCC? Why do SCCs matter?

Easy way to remember: `mutually reachable group`

### 43. Which SCC algorithm did you use?

Short answer: Tarjan's algorithm.

Detailed answer: The theory layer uses Tarjan's algorithm to identify SCCs by tracking DFS index values, lowlink values, and a stack of active nodes.

Follow-up questions: Why Tarjan? Could Kosaraju also work?

Easy way to remember: `DFS + lowlink + stack`

### 44. Why are SCCs important in this project?

Short answer: Because they identify cyclic regions.

Detailed answer: SCCs help isolate where cycles exist in the graph. This is useful because the paper's cyclic chaining problem becomes easier to reason about once strongly connected regions are identified.

Follow-up questions: How do SCCs relate to DAG conversion? Why not ignore cycles?

Easy way to remember: `SCC tells where loops live`

### 45. What is a DAG?

Short answer: A directed acyclic graph.

Detailed answer: A DAG is a directed graph with no cycles. DAGs are much easier for ordered dynamic programming because they admit topological ordering.

Follow-up questions: Why is topological order useful? Why convert cyclic structure?

Easy way to remember: `directed and no loops`

### 46. Why convert a cyclic graph into a DAG approximation?

Short answer: To make downstream ordering and chaining easier to explain and compute.

Detailed answer: Cycles break simple ordering. By removing back edges in an educational approximation, the project produces a DAG where path cover and chaining ideas become easier to demonstrate.

Follow-up questions: Is this done on the real graph? Does this equal PanAligner internals?

Easy way to remember: `remove loops to simplify order`

### 47. How do you find back edges in the theory demo?

Short answer: With DFS color states.

Detailed answer: The demo marks nodes white, gray, and black during depth-first search. An edge to a gray node is treated as a back edge, which indicates a cycle in the DFS tree.

Follow-up questions: What do gray nodes mean? Why are back edges important?

Easy way to remember: `edge to gray = back edge`

### 48. Does your project remove back edges from the real production graph?

Short answer: No, only from the educational copy.

Detailed answer: The back-edge removal and DAG conversion happen in the theory layer, not in the real PanAligner alignment graph used for practical evaluation.

Follow-up questions: Why keep them separate? Is this an honest reproduction?

Easy way to remember: `theory only, not production graph`

### 49. Why add synthetic back edges in the theory layer?

Short answer: To demonstrate the cyclic case when real graphs are mostly DAG-like.

Detailed answer: The current real graphs for these loci often do not show strong cyclic structure. To still explain the paper's cyclic algorithms, the project adds synthetic back edges to a representative educational graph only.

Follow-up questions: Which file does this? Did you modify real data?

Easy way to remember: `teach cycles without changing real graph`

### 50. What does it mean if the largest SCC size is 1?

Short answer: There are no nontrivial strongly connected regions.

Detailed answer: If every SCC has size 1 and there are no self-loops, the graph is effectively acyclic at the node level.

Follow-up questions: Does that mean SCC analysis is useless? Why still include it?

Easy way to remember: `singletons only`

## Path Cover, Anchors, and Precedence

### 51. What is a path cover?

Short answer: A set of paths that covers all nodes in a DAG.

Detailed answer: In a DAG, a path cover is a collection of directed paths such that every node belongs to at least one path. It helps structure the graph for ordered reasoning.

Follow-up questions: Why is it useful? Is your method exact?

Easy way to remember: `paths that cover all nodes`

### 52. What path cover method do you use?

Short answer: A simplified greedy path cover.

Detailed answer: The theory demo uses a greedy method on the DAG approximation. It repeatedly starts from an uncovered source-like node and extends forward through uncovered successors.

Follow-up questions: Is this the same as the full paper method? Why use greedy?

Easy way to remember: `simple greedy cover for teaching`

### 53. Is your path cover implementation the exact PanAligner preprocessing algorithm?

Short answer: No.

Detailed answer: It is an educational analogue designed to illustrate the concept. The real tool's internals are more advanced and optimized.

Follow-up questions: Why still include it? What is the benefit in viva?

Easy way to remember: `conceptual, not exact`

### 54. What is an anchor?

Short answer: A local match between the query and the graph.

Detailed answer: In the educational model, an anchor is represented by the graph vertex, graph interval, query interval, and a weight indicating match usefulness.

Follow-up questions: Why are anchors important? What does weight mean?

Easy way to remember: `small local match`

### 55. How is an anchor represented in your theory demo?

Short answer: `(vertex, [x..y], [c..d], weight)`

Detailed answer: `vertex` is the graph node, `[x..y]` is the interval on that node, `[c..d]` is the interval on the query, and `weight` is the anchor score.

Follow-up questions: Which script creates anchors? Are there noise anchors?

Easy way to remember: `where on graph, where on query, how good`

### 56. Why do you add distractor anchors?

Short answer: To show that chaining must choose wisely.

Detailed answer: The theory demo includes a few synthetic noise anchors so the DP does not just have one obvious path. This makes the chaining explanation more realistic and educational.

Follow-up questions: Where are they added? Why is this helpful?

Easy way to remember: `bad choices make good demo`

### 57. What is precedence between anchors?

Short answer: A rule saying one anchor can validly come before another.

Detailed answer: In the simplified model, anchor A precedes anchor B if query coordinates increase and graph reachability allows A to come before B.

Follow-up questions: Why need precedence? Is query order alone enough?

Easy way to remember: `valid before-after relation`

### 58. How do you check precedence in different vertices?

Short answer: By graph reachability.

Detailed answer: If the two anchors are on different vertices, anchor A can precede anchor B only if B's vertex is reachable from A's vertex in the graph.

Follow-up questions: How is reachability computed? What if there is no path?

Easy way to remember: `can we get there in graph?`

### 59. How do you treat same-vertex precedence in the simplified model?

Short answer: It is allowed if intervals are ordered, or in cyclic cases when appropriate.

Detailed answer: If two anchors are on the same vertex, the simplified rule checks whether the graph interval ordering makes sense, and in cyclic cases it can allow same-vertex precedence more flexibly.

Follow-up questions: Why special handling? What if the graph is cyclic?

Easy way to remember: `same vertex needs special rule`

### 60. What file computes precedence?

Short answer: `scripts/theory/precedence_demo.py`

Detailed answer: That script builds a precedence graph over the generated anchors and saves both a visualization and a text report.

Follow-up questions: What output file does it create? How is it used in chaining?

Easy way to remember: `precedence_demo builds anchor order graph`

## Co-linear Chaining

### 61. What is co-linear chaining?

Short answer: Choosing the best ordered chain of compatible anchors.

Detailed answer: It is the process of selecting local matches that preserve query order and graph consistency while maximizing score and penalizing unreasonable gaps.

Follow-up questions: Why not take all anchors? What makes anchors compatible?

Easy way to remember: `best ordered anchor sequence`

### 62. Why is chaining difficult on graphs?

Short answer: Because graph order is not as simple as linear order.

Detailed answer: In a linear sequence, order is obvious. In graphs, especially cyclic ones, multiple paths may exist and a vertex may be revisited, so compatibility requires graph reachability and careful ordering logic.

Follow-up questions: Why are cycles especially hard? What role do SCCs play?

Easy way to remember: `many paths, maybe loops`

### 63. What DP recurrence do you use in the simplified demo?

Short answer: `dp[j] = best chain score ending at anchor j`

Detailed answer: For each anchor `j`, the demo starts with the anchor's own weight, then checks all earlier anchors `i` that can precede `j`. If `dp[i] + weight(j) - gap(i,j)` is better, it updates `dp[j]`.

Follow-up questions: Why dynamic programming? What is stored besides scores?

Easy way to remember: `best ending here`

### 64. What is a gap cost?

Short answer: A penalty between two consecutive anchors.

Detailed answer: The gap cost reflects how far apart two anchors are in the query and graph. Larger or less consistent gaps reduce the desirability of chaining those anchors together.

Follow-up questions: What are query gap and graph gap? Why penalize gaps?

Easy way to remember: `distance penalty`

### 65. What is query gap?

Short answer: Distance between anchor intervals on the query.

Detailed answer: It measures how many query positions lie between the end of one anchor and the start of the next.

Follow-up questions: What if anchors overlap? Is negative gap allowed?

Easy way to remember: `space on the read`

### 66. What is graph gap?

Short answer: Distance between anchor positions in the graph.

Detailed answer: In the educational demo, it is approximated using shortest-path distance between anchor vertices plus small prefix and suffix adjustments based on local positions.

Follow-up questions: Is this the full PanAligner graph gap? Why approximate?

Easy way to remember: `space on the graph`

### 67. What is the best chain found in your current theory output?

Short answer: `A6 -> A7 -> A8`

Detailed answer: According to `outputs/theory/chaining_results.txt`, the simplified DP selected anchors A6, A7, and A8 as the best consistent chain for the current educational anchor set, with best score `12.00`.

Follow-up questions: Why not include A1 to A5? What role do gap costs play?

Easy way to remember: `best chain starts at A6`

### 68. What is the current best chain score?

Short answer: `12.00`

Detailed answer: The best chain score in the current educational output is 12.00 after accounting for anchor weights and gap penalties.

Follow-up questions: Does score depend only on weights? How can score decrease?

Easy way to remember: `best score twelve`

### 69. Why not say your simplified DP is the full paper algorithm?

Short answer: Because it is not.

Detailed answer: The educational DP is intentionally simplified for clarity. The real PanAligner uses more advanced internal logic, data structures, and optimizations than the demo.

Follow-up questions: Then why include the DP? Is it still useful?

Easy way to remember: `teach concept, not exact internals`

### 70. What is iterative convergence in your chaining demo?

Short answer: Repeated score updates until values stop changing.

Detailed answer: The demo also runs a score-propagation style process across the precedence graph to illustrate how chaining values can stabilize after repeated updates.

Follow-up questions: Why include it? Does PanAligner use exactly this?

Easy way to remember: `update until stable`

## PanAligner and Practical Alignment

### 71. What is PanAligner?

Short answer: A sequence-to-pangenome-graph aligner.

Detailed answer: PanAligner is the real tool used in the project to align held-out query FASTA sequences to graphs stored in GFA format. It outputs alignments in GAF format.

Follow-up questions: Who builds the graph then? What format does it output?

Easy way to remember: `aligner, not graph builder`

### 72. What command pattern do you use to run PanAligner?

Short answer: `PanAligner -t 4 -cx lr graph.gfa query.fa > out.gaf`

Detailed answer: The wrapper script constructs a command that passes thread count, long-read preset, the graph file, and the query file, and saves stdout as a GAF file.

Follow-up questions: What does `-cx lr` mean in practice? Why save to GAF?

Easy way to remember: `graph + query -> gaf`

### 73. Which file wraps PanAligner execution?

Short answer: `scripts/align.py`

Detailed answer: The `run_panaligner(...)` function in `align.py` runs the external binary and immediately parses the generated GAF file for later use.

Follow-up questions: Does it choose best alignment? Which parser is used?

Easy way to remember: `align.py runs the aligner`

### 74. What is the role of `paper_evaluation.py`?

Short answer: It performs held-out alignment evaluation.

Detailed answer: It reads the test manifest and train graph manifest, runs PanAligner for each held-out query, parses the best alignment, stores records, and computes overall and per-gene metrics.

Follow-up questions: What outputs are produced? Which metrics are reported?

Easy way to remember: `evaluate held-out queries`

### 75. Why align to train-derived graphs?

Short answer: To avoid evaluating on sequences already used to build the graph.

Detailed answer: Graphs are built from the training set only, and evaluation queries come from the test set. This makes the evaluation more meaningful and closer to an unseen-query scenario.

Follow-up questions: Is this classification? What would happen if we used all sequences in the graph?

Easy way to remember: `unseen queries on train graph`

### 76. What is the current overall alignment rate?

Short answer: `1.0`

Detailed answer: The current evaluation output shows that all 46 held-out test queries aligned successfully to their corresponding combined train graphs.

Follow-up questions: Is this surprising? What does it say about the dataset?

Easy way to remember: `46 out of 46`

### 77. What is the current mean identity?

Short answer: Approximately `0.99955`

Detailed answer: The mean identity across aligned held-out queries is extremely close to 1, showing very high sequence similarity between the queries and the chosen graph paths.

Follow-up questions: Why so high? Does identity alone prove everything?

Easy way to remember: `almost perfect identity`

### 78. What is the current mean MAPQ?

Short answer: `60`

Detailed answer: The current output reports mean mapping quality 60, which generally indicates very confident alignments in this evaluation set.

Follow-up questions: Does MAPQ depend on the aligner? Why also inspect identity and coverage?

Easy way to remember: `max-looking confidence`

### 79. Which graph is used for each evaluation query?

Short answer: The corresponding gene's combined train graph.

Detailed answer: APP test sequences align to the APP combined train graph, PSEN1 test sequences to the PSEN1 combined train graph, and similarly for PSEN2.

Follow-up questions: Why not align every query to every graph? Could that be another experiment?

Easy way to remember: `same gene, combined train graph`

### 80. What practical outputs are produced by evaluation?

Short answer: GAF files, JSON metrics, JSON records, and a text report.

Detailed answer: The evaluation creates per-query alignment folders, `alignment_metrics.json`, `alignment_records.json`, and `evaluation_report.txt`, and also example alignment visualizations.

Follow-up questions: Where are these stored? Why store both summary and detailed records?

Easy way to remember: `raw alignments + summaries`

## Theory Layer vs Real Layer

### 81. What is the single most important distinction in your project?

Short answer: Real PanAligner alignment versus simplified educational demos.

Detailed answer: The real pipeline uses official tools to build graphs and align queries, while the theory layer is a separate Python educational reproduction of the paper's main ideas for explanation and viva preparation.

Follow-up questions: Why is this distinction important? Did you modify PanAligner?

Easy way to remember: `real tool, demo theory`

### 82. Did you reimplement the full PanAligner internals?

Short answer: No.

Detailed answer: The project uses the real PanAligner binary for alignment. The Python theory scripts only provide simplified conceptual reproductions to make the paper understandable.

Follow-up questions: Then what exactly did you implement? Why still call it reproduction?

Easy way to remember: `used real aligner, explained theory ourselves`

### 83. Why is the theory layer still valuable if it is simplified?

Short answer: Because it makes the paper understandable.

Detailed answer: The full aligner internals are optimized and harder to present directly in viva. The theory layer isolates and visualizes the central concepts so they can be taught clearly without pretending to replace the real tool.

Follow-up questions: Which concepts are covered? How does this help defense?

Easy way to remember: `clarity over complexity`

### 84. Did the theory layer affect the real PanAligner results?

Short answer: No.

Detailed answer: The theory scripts are separate and educational. They do not alter the actual train graphs used for evaluation or the PanAligner binary's behavior.

Follow-up questions: Then how are they connected? Which graph do they read?

Easy way to remember: `separate for explanation only`

### 85. Why do you say your project is honest?

Short answer: Because it clearly separates what is real and what is educational.

Detailed answer: The project does not falsely claim the simplified DP or SCC demos are the actual PanAligner internals. It explicitly states that the real binary performs alignment and the theory layer explains the ideas conceptually.

Follow-up questions: Why is this important in academic work? How would you defend it?

Easy way to remember: `clear scope, no false claims`

## Graph Analysis and Visualization

### 86. Why visualize the graphs?

Short answer: To understand structure and communicate it clearly.

Detailed answer: Visualization helps reveal branching, connectivity, and highlighted alignment paths. It also makes concepts like SCCs, path covers, and traversed nodes much easier to explain in viva.

Follow-up questions: Which script handles it? What library is used?

Easy way to remember: `see structure, explain better`

### 87. Which library is used for graph modeling in the Python demos?

Short answer: `networkx`

Detailed answer: `networkx` is used in the visualization and theory modules to represent directed graphs, compute SCCs and reachability, and create graph plots.

Follow-up questions: Why choose networkx? Is it used in PanAligner internals?

Easy way to remember: `networkx for readable graph demos`

### 88. Which library is used for plotting?

Short answer: `matplotlib`

Detailed answer: The project uses matplotlib to save graph and chaining plots such as graph topology images, path cover images, best-chain plots, and convergence plots.

Follow-up questions: Why plot to PNG? Where are outputs stored?

Easy way to remember: `matplotlib makes the figures`

### 89. What is highlighted in alignment visualization?

Short answer: Nodes traversed by the alignment.

Detailed answer: The project parses traversed graph nodes from the best GAF path and highlights those nodes when drawing the graph, giving an intuitive view of where the query aligned.

Follow-up questions: How are traversed nodes extracted? Why is this useful?

Easy way to remember: `highlight the route`

### 90. What is reachability in your theory context?

Short answer: Which nodes can be reached from each node.

Detailed answer: The theory layer precomputes descendants for each node so that precedence logic can quickly check whether one anchor's vertex can reach another anchor's vertex.

Follow-up questions: Where is reachability used? Why not recompute every time?

Easy way to remember: `who can reach whom`

## Debugging and Practical Questions

### 91. If a query does not align, where would you debug first?

Short answer: Input files, graph choice, and GAF output.

Detailed answer: I would first verify that the query FASTA is valid, confirm the correct graph path and PanAligner binary are being used, inspect the produced GAF file, and then check whether the query belongs to the expected gene and whether the graph was built correctly.

Follow-up questions: Would you inspect preprocessing or graph construction next? What if GAF is empty?

Easy way to remember: `query, graph, aligner, output`

### 92. If GAF parsing fails, what might be the reason?

Short answer: Malformed GAF structure or unexpected fields.

Detailed answer: The parser expects at least 12 tab-separated fields. If the line is incomplete, empty, corrupted, or in an unexpected format, parsing can fail.

Follow-up questions: Which function checks this? What exception is raised?

Easy way to remember: `need the required 12 fields`

### 93. If a graph is missing, which manifest would you inspect?

Short answer: `train_graph_manifest.json` or `graph_manifest.json`

Detailed answer: These manifest files store the resolved graph paths for each gene. They are the best place to verify whether the graphs were built and where they were saved.

Follow-up questions: Which script creates them? Are there train and non-train versions?

Easy way to remember: `manifest tells path truth`

### 94. If train/test counts look wrong, which file would you inspect?

Short answer: `split_manifest.json` and `split_summary.txt`

Detailed answer: These files record the split counts per gene and class, making them the first place to verify whether the split was performed correctly.

Follow-up questions: Which script writes them? What is the default seed?

Easy way to remember: `check split metadata first`

### 95. If graph visualization is empty, what might be the reason?

Short answer: The parsed graph may have zero nodes or the file may not contain expected GFA records.

Detailed answer: An empty or malformed GFA, wrong file path, or parsing issue with segment/link lines could lead to an empty networkx graph and therefore no meaningful visualization.

Follow-up questions: Which parser is used? How does the script handle empty graphs?

Easy way to remember: `no parsed nodes, no plot`

### 96. If your theory demo shows no cycles, what would you say in viva?

Short answer: The real graph is mostly DAG-like, so the project adds synthetic back edges only in the educational copy when needed.

Detailed answer: That is an honest and expected outcome for the current dataset. The paper studies cyclic graphs in general, so the educational layer creates a synthetic cyclic example from a representative real subgraph to teach those ideas clearly.

Follow-up questions: Why not force cycles into production evaluation? Is that honest?

Easy way to remember: `real may be acyclic, theory adds cycle`

### 97. Why store outputs in JSON as well as text?

Short answer: JSON is machine-readable and text is easy for humans.

Detailed answer: JSON files help scripts and later processing, while text reports are convenient for quick review, presentation, and viva discussion.

Follow-up questions: Which outputs are JSON? Which are text?

Easy way to remember: `JSON for tools, text for people`

### 98. Why save one-sample FASTA files after preprocessing?

Short answer: They make later stages simpler and more modular.

Detailed answer: Separate FASTA files are convenient for train/test splitting, per-query alignment, and clean organization by gene and class.

Follow-up questions: Could PanAligner align combined FASTA instead? Why use one-sample files?

Easy way to remember: `one query file per sample`

### 99. Why is reproducibility important in this project?

Short answer: Because a paper reproduction project should be repeatable.

Detailed answer: Reproducibility allows the same split, graph build, evaluation behavior, and reports to be generated again, which is important for both scientific credibility and viva defense.

Follow-up questions: How is reproducibility achieved? Which parameters matter?

Easy way to remember: `same inputs, same results`

### 100. What would you improve if you had more time?

Short answer: Larger datasets, richer anchor generation, and closer comparison with PanAligner internals.

Detailed answer: I would extend the evaluation to larger and more cyclic graph datasets, generate anchors from more realistic seed hits instead of only educational examples, and compare theory-layer outputs more directly with selected PanAligner alignment traces.

Follow-up questions: Why these three improvements? What is the current limitation?

Easy way to remember: `bigger data, better anchors, deeper comparison`

## Paper Understanding

### 101. What problem does the paper address?

Short answer: Co-linear chaining on pangenome graphs, especially cyclic ones.

Detailed answer: The paper studies how to select consistent chains of anchors on graph-based genomes where graph topology may include cycles, making the ordering problem more difficult than in linear references or simple DAGs.

Follow-up questions: Why are cycles hard? Why is chaining central?

Easy way to remember: `chaining becomes hard in graph loops`

### 102. Why are cyclic pangenome graphs challenging?

Short answer: Because loops break simple ordering.

Detailed answer: In cyclic graphs, a node or region may be revisited and the graph does not admit a straightforward topological ordering, so anchor precedence and dynamic programming become more subtle.

Follow-up questions: How do SCCs help? Why use a DAG approximation in the demo?

Easy way to remember: `loops break easy order`

### 103. How does your project reflect the paper's ideas?

Short answer: Through both real alignment workflow and educational graph demos.

Detailed answer: The real workflow demonstrates graph construction and graph alignment, while the theory layer reproduces concepts such as SCC detection, cyclic handling, path cover, precedence, and simplified co-linear chaining.

Follow-up questions: Which part is exact and which part is simplified?

Easy way to remember: `workflow + concepts`

### 104. Does your project claim to reproduce the full algorithmic complexity of the paper?

Short answer: No.

Detailed answer: The project reproduces the workflow and core concepts in a practical and educational way, but it does not claim that the simplified Python theory modules exactly match the official implementation's internal complexity or performance.

Follow-up questions: Why is that still valid? What is the benefit?

Easy way to remember: `conceptually faithful, not internally identical`

### 105. What is the educational value of your project?

Short answer: It makes a difficult paper teachable and defendable.

Detailed answer: Many graph-alignment papers are hard to explain directly from optimized source code. This project bridges that gap by combining real tools with clear diagrams, reports, and simplified demos.

Follow-up questions: Which concepts become easier to teach? Why is that important in viva?

Easy way to remember: `turns hard paper into understandable project`

## Short Rapid-fire Questions

### 106. Who builds the graph?

Short answer: `minigraph`

Detailed answer: minigraph constructs the GFA graph from the reference and sample FASTA files.

Follow-up questions: Then who aligns queries? Which format does it output?

Easy way to remember: `builder = minigraph`

### 107. Who aligns the query to the graph?

Short answer: `PanAligner`

Detailed answer: PanAligner is the real query-to-graph aligner in this project and outputs GAF alignments.

Follow-up questions: Which wrapper calls it? What does GAF store?

Easy way to remember: `aligner = PanAligner`

### 108. What format stores the graph?

Short answer: `GFA`

Detailed answer: GFA stores graph segments and links, and is the format read by visualization and alignment stages.

Follow-up questions: What tool creates it? Which lines define nodes?

Easy way to remember: `graph format = GFA`

### 109. What format stores the alignment?

Short answer: `GAF`

Detailed answer: GAF stores the query's alignment to the graph, including path and alignment quality information.

Follow-up questions: How do you compute identity from it? Which script parses it?

Easy way to remember: `alignment format = GAF`

### 110. What is the default split seed?

Short answer: `42`

Detailed answer: The split is reproducibly randomized using seed 42 by default in `main.py` and `split_dataset.py`.

Follow-up questions: Why use a seed? Can it be changed?

Easy way to remember: `answer to reproducibility = 42`

### 111. What is the default test fraction?

Short answer: `0.20`

Detailed answer: Twenty percent of each class per gene is assigned to the test set by default.

Follow-up questions: What is the train fraction then? How is class balance preserved?

Easy way to remember: `20 percent held out`

### 112. What is the current best demo chain?

Short answer: `A6 -> A7 -> A8`

Detailed answer: The current simplified chaining output selects A6, A7, and A8 as the best chain.

Follow-up questions: What is the best score? Why those anchors?

Easy way to remember: `A6-A7-A8`

### 113. What is the current overall alignment rate?

Short answer: `100%`

Detailed answer: All 46 held-out test queries aligned successfully in the current evaluation report.

Follow-up questions: Does that mean perfect biological truth? What are the other metrics?

Easy way to remember: `46 of 46`

### 114. What is the current mean identity?

Short answer: About `0.99955`

Detailed answer: The current overall mean identity is extremely high, meaning the aligned query and chosen graph path match closely.

Follow-up questions: What about coverage and MAPQ? Why all three?

Easy way to remember: `identity almost one`

### 115. What is the most viva-safe summary sentence?

Short answer: We used the real PanAligner for alignment and separate simplified demos for explanation.

Detailed answer: This sentence is the safest because it clearly distinguishes real implementation from educational reproduction and prevents overstating what the demo layer does.

Follow-up questions: Why is this distinction so important? How would you defend it?

Easy way to remember: `real aligner, demo explainer`

## Final Viva Closing Answers

### 116. What is the biggest strength of your project?

Short answer: It combines practical execution with clear conceptual understanding.

Detailed answer: Many projects either run tools without deep understanding or explain theory without practical reproduction. This project does both, which makes it stronger for learning and for viva defense.

Follow-up questions: Which side was harder? Why is this useful academically?

Easy way to remember: `runs and explains`

### 117. What is the biggest limitation of your project?

Short answer: The theory layer is simplified and the current real graphs are often DAG-like.

Detailed answer: While the practical alignment layer is real, the educational demos simplify the theory, and the current dataset does not always naturally show the strongest cyclic behavior, so part of the cyclic demonstration is synthetic.

Follow-up questions: Why is that acceptable? What future work would address it?

Easy way to remember: `simplified theory, limited cyclic data`

### 118. Why should your viva panel consider the project successful?

Short answer: Because it is technically correct, reproducible, and clearly explainable.

Detailed answer: The project builds real pangenome graphs, performs real held-out graph alignment, produces strong evaluation outputs, and explains the underlying paper concepts honestly and clearly, which is exactly what a strong reproduction project should do.

Follow-up questions: What evidence supports success? Which outputs prove it?

Easy way to remember: `correct, reproducible, explainable`

### 119. If asked whether you exactly reproduced PanAligner internals, what should you say?

Short answer: No, but we reproduced the workflow and key concepts faithfully.

Detailed answer: I would say that the practical alignment uses the official PanAligner binary, while the theory layer is a simplified educational reproduction of the paper concepts. So the workflow is real and the explanation is faithful, but the demo is not a literal reimplementation of every internal optimization.

Follow-up questions: Why is this still valuable? What is honest academic framing?

Easy way to remember: `official tool + faithful teaching`

### 120. What is your final 20-second viva summary?

Short answer: This project reproduces the PanAligner paper by building pangenome graphs, aligning held-out sequences with the real PanAligner, and explaining cyclic graph chaining concepts with simple theory demos.

Detailed answer: Our project starts from FASTA preprocessing, creates reproducible train/test splits, builds healthy, unhealthy, and combined pangenome graphs with minigraph, aligns held-out sequences to train-derived graphs using the real PanAligner, parses GAF alignments for evaluation, and separately demonstrates SCCs, DAG conversion, path cover, precedence, and simplified co-linear chaining for clear viva explanation.

Follow-up questions: Which result would you show first? What is the key distinction in the project?

Easy way to remember: `preprocess, build, align, explain`
