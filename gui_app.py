from __future__ import annotations

import json
import os
import queue
import subprocess
import sys
import threading
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk


ROOT = Path(__file__).resolve().parent
MAIN_SCRIPT = ROOT / "main.py"
OUTPUTS_DIR = ROOT / "outputs"
DEFAULT_GRAPH_MANIFEST = ROOT / "data" / "metadata" / "graph_manifest.json"
DEFAULT_CUSTOM_OUTPUT = OUTPUTS_DIR / "alignments" / "custom_query"


def open_path(target: Path) -> None:
    resolved = target.resolve()
    if sys.platform.startswith("win"):
        os.startfile(str(resolved))
        return
    if sys.platform == "darwin":
        subprocess.Popen(["open", str(resolved)])
        return
    subprocess.Popen(["xdg-open", str(resolved)])


class PanAlignerGUI:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("PanAligner Project Console")
        self.root.geometry("1480x930")
        self.root.minsize(1240, 780)

        self.process: subprocess.Popen[str] | None = None
        self.output_queue: queue.Queue[tuple[str, str]] = queue.Queue()
        self.output_buffer: list[str] = []
        self.preview_image: tk.PhotoImage | None = None
        self.last_summary: dict | None = None
        self.image_paths: list[Path] = []
        self.file_paths: list[Path] = []

        self.pipeline_mode = tk.StringVar(value="--custom-query-analysis")
        self.query_source = tk.StringVar(value="sequence")
        self.threads_var = tk.StringVar(value="4")
        self.split_seed_var = tk.StringVar(value="42")
        self.test_fraction_var = tk.StringVar(value="0.20")
        self.minigraph_var = tk.StringVar(value="")
        self.panaligner_var = tk.StringVar(value="")
        self.theory_graph_var = tk.StringVar(value="")
        self.input_fastas_var = tk.StringVar(value="")

        self.query_fasta_var = tk.StringVar(value="")
        self.query_argument_var = tk.StringVar(value="")
        self.query_gene_var = tk.StringVar(value="")
        self.query_output_var = tk.StringVar(value=str(DEFAULT_CUSTOM_OUTPUT))
        self.graph_manifest_var = tk.StringVar(value=str(DEFAULT_GRAPH_MANIFEST))

        self.status_var = tk.StringVar(value="Ready.")

        self._configure_style()
        self._build_layout()
        self._refresh_outputs()
        self.root.after(150, self._poll_output)

    def _configure_style(self) -> None:
        style = ttk.Style()
        if "clam" in style.theme_names():
            style.theme_use("clam")
        self.root.configure(bg="#F5EFE6")
        style.configure("Title.TLabel", font=("Segoe UI Semibold", 18), foreground="#1F3A5F", background="#F5EFE6")
        style.configure("Section.TLabelframe", background="#FFFDF9")
        style.configure("Section.TLabelframe.Label", font=("Segoe UI Semibold", 11), foreground="#264653")
        style.configure("Primary.TButton", font=("Segoe UI Semibold", 10))
        style.configure("Status.TLabel", background="#F5EFE6", foreground="#5A4A42", font=("Segoe UI", 10))

    def _build_layout(self) -> None:
        header = ttk.Frame(self.root, padding=(18, 16))
        header.pack(fill="x")
        ttk.Label(header, text="PanAligner Project Console", style="Title.TLabel").pack(anchor="w")
        ttk.Label(
            header,
            text="Run the pipeline, inspect evaluation scores, read reports, and preview generated graphs in one dashboard.",
            background="#F5EFE6",
            foreground="#5A4A42",
        ).pack(anchor="w", pady=(4, 0))

        body = ttk.Panedwindow(self.root, orient="horizontal")
        body.pack(fill="both", expand=True, padx=16, pady=(0, 16))

        controls = ttk.Frame(body, padding=8)
        results = ttk.Frame(body, padding=8)
        body.add(controls, weight=4)
        body.add(results, weight=7)

        self._build_controls(controls)
        self._build_results(results)

    def _build_controls(self, parent: ttk.Frame) -> None:
        notebook = ttk.Notebook(parent)
        notebook.pack(fill="both", expand=True)

        pipeline_tab = ttk.Frame(notebook, padding=12)
        custom_tab = ttk.Frame(notebook, padding=12)
        notebook.add(pipeline_tab, text="Project Modes")
        notebook.add(custom_tab, text="Custom Query")

        self._build_pipeline_tab(pipeline_tab)
        self._build_custom_tab(custom_tab)

    def _build_pipeline_tab(self, parent: ttk.Frame) -> None:
        mode_frame = ttk.LabelFrame(parent, text="Run Mode", style="Section.TLabelframe", padding=12)
        mode_frame.pack(fill="x", pady=(0, 10))

        for text, value in (
            ("Full Pipeline", "--full-pipeline"),
            ("Theory Only", "--theory-only"),
            ("Evaluate", "--evaluate"),
        ):
            ttk.Radiobutton(mode_frame, text=text, value=value, variable=self.pipeline_mode).pack(anchor="w", pady=2)

        settings = ttk.LabelFrame(parent, text="Settings", style="Section.TLabelframe", padding=12)
        settings.pack(fill="x", pady=(0, 10))

        self._labeled_entry(settings, "Threads", self.threads_var, 0)
        self._labeled_entry(settings, "Split Seed", self.split_seed_var, 1)
        self._labeled_entry(settings, "Test Fraction", self.test_fraction_var, 2)
        self._labeled_file_row(settings, "Minigraph Binary", self.minigraph_var, 3, filetypes=[("Executable", "*")])
        self._labeled_file_row(settings, "PanAligner Binary", self.panaligner_var, 4, filetypes=[("Executable", "*")])
        self._labeled_file_row(settings, "Theory Graph", self.theory_graph_var, 5, filetypes=[("GFA Files", "*.gfa"), ("All Files", "*.*")])

        input_frame = ttk.LabelFrame(parent, text="Optional Input FASTAs", style="Section.TLabelframe", padding=12)
        input_frame.pack(fill="x", pady=(0, 10))
        ttk.Label(input_frame, text="Leave blank to use the default *_combined.fasta files in the project root.").pack(anchor="w")
        ttk.Entry(input_frame, textvariable=self.input_fastas_var).pack(fill="x", pady=(6, 4))
        ttk.Label(input_frame, text="Use semicolons between paths when providing multiple FASTA files.").pack(anchor="w")

        actions = ttk.Frame(parent)
        actions.pack(fill="x", pady=(4, 0))
        ttk.Button(actions, text="Run Selected Mode", style="Primary.TButton", command=self.run_pipeline_mode).pack(side="left")
        ttk.Button(actions, text="Stop Current Run", command=self.stop_run).pack(side="left", padx=8)
        ttk.Button(actions, text="Refresh Outputs", command=self._refresh_outputs).pack(side="left")

    def _build_custom_tab(self, parent: ttk.Frame) -> None:
        source_frame = ttk.LabelFrame(parent, text="Query Source", style="Section.TLabelframe", padding=12)
        source_frame.pack(fill="x", pady=(0, 10))

        radio_row = ttk.Frame(source_frame)
        radio_row.pack(fill="x", pady=(0, 8))
        ttk.Radiobutton(radio_row, text="Inline DNA Sequence", value="sequence", variable=self.query_source, command=self._toggle_query_source).pack(side="left")
        ttk.Radiobutton(radio_row, text="Query FASTA File", value="fasta", variable=self.query_source, command=self._toggle_query_source).pack(side="left", padx=14)

        self.query_sequence_text = tk.Text(source_frame, height=8, wrap="word", font=("Consolas", 10))
        self.query_sequence_text.pack(fill="x")

        fasta_row = ttk.Frame(source_frame)
        fasta_row.pack(fill="x", pady=(10, 0))
        ttk.Label(fasta_row, text="Query FASTA").pack(anchor="w")
        entry_row = ttk.Frame(fasta_row)
        entry_row.pack(fill="x", pady=(4, 0))
        self.query_fasta_entry = ttk.Entry(entry_row, textvariable=self.query_fasta_var)
        self.query_fasta_entry.pack(side="left", fill="x", expand=True)
        ttk.Button(entry_row, text="Browse", command=lambda: self._pick_file(self.query_fasta_var, [("FASTA Files", "*.fa *.fasta"), ("All Files", "*.*")])).pack(side="left", padx=(8, 0))

        details = ttk.LabelFrame(parent, text="Analysis Details", style="Section.TLabelframe", padding=12)
        details.pack(fill="x", pady=(0, 10))
        self._labeled_entry(details, "Threads", self.threads_var, 0)
        self._labeled_entry(details, "Gene Override", self.query_gene_var, 1)
        self._labeled_entry(details, "Custom Argument / Note", self.query_argument_var, 2)
        self._labeled_file_row(details, "Graph Manifest", self.graph_manifest_var, 3, filetypes=[("JSON Files", "*.json"), ("All Files", "*.*")])
        self._labeled_dir_row(details, "Output Directory", self.query_output_var, 4)
        self._labeled_file_row(details, "PanAligner Binary", self.panaligner_var, 5, filetypes=[("Executable", "*")])

        actions = ttk.Frame(parent)
        actions.pack(fill="x", pady=(4, 0))
        ttk.Button(actions, text="Run Custom Query Analysis", style="Primary.TButton", command=self.run_custom_query).pack(side="left")
        ttk.Button(actions, text="Open Output Folder", command=self._open_output_dir).pack(side="left", padx=8)

        self._toggle_query_source()

    def _build_results(self, parent: ttk.Frame) -> None:
        status_frame = ttk.Frame(parent)
        status_frame.pack(fill="x", pady=(0, 10))
        ttk.Label(status_frame, textvariable=self.status_var, style="Status.TLabel").pack(side="left")
        ttk.Button(status_frame, text="Refresh Dashboard", command=self._refresh_outputs).pack(side="right")

        notebook = ttk.Notebook(parent)
        notebook.pack(fill="both", expand=True)

        overview_tab = ttk.Frame(notebook, padding=10)
        metrics_tab = ttk.Frame(notebook, padding=10)
        visuals_tab = ttk.Frame(notebook, padding=10)
        files_tab = ttk.Frame(notebook, padding=10)
        logs_tab = ttk.Frame(notebook, padding=10)

        notebook.add(overview_tab, text="Overview")
        notebook.add(metrics_tab, text="Evaluation Scores")
        notebook.add(visuals_tab, text="Graphs")
        notebook.add(files_tab, text="Reports & Files")
        notebook.add(logs_tab, text="Console")

        self._build_overview_tab(overview_tab)
        self._build_metrics_tab(metrics_tab)
        self._build_visuals_tab(visuals_tab)
        self._build_files_tab(files_tab)
        self._build_logs_tab(logs_tab)

    def _build_overview_tab(self, parent: ttk.Frame) -> None:
        top = ttk.Panedwindow(parent, orient="vertical")
        top.pack(fill="both", expand=True)

        summary_frame = ttk.LabelFrame(top, text="Run Summary", style="Section.TLabelframe", padding=10)
        artifacts_frame = ttk.LabelFrame(top, text="Detected Artifacts", style="Section.TLabelframe", padding=10)
        top.add(summary_frame, weight=3)
        top.add(artifacts_frame, weight=2)

        self.summary_text = tk.Text(summary_frame, wrap="word", font=("Segoe UI", 10), height=12, bg="#FCFBF8")
        self.summary_text.pack(fill="both", expand=True)

        actions = ttk.Frame(summary_frame)
        actions.pack(fill="x", pady=(8, 0))
        ttk.Button(actions, text="Open Primary Output", command=self._open_primary_output).pack(side="left")

        artifact_columns = ("path", "kind")
        self.artifacts_tree = ttk.Treeview(artifacts_frame, columns=artifact_columns, show="headings", height=10)
        self.artifacts_tree.heading("path", text="Artifact")
        self.artifacts_tree.heading("kind", text="Type")
        self.artifacts_tree.column("path", width=520, anchor="w")
        self.artifacts_tree.column("kind", width=120, anchor="center")
        self.artifacts_tree.pack(side="left", fill="both", expand=True)
        self.artifacts_tree.bind("<Double-1>", lambda _event: self._open_selected_artifact())
        artifact_scroll = ttk.Scrollbar(artifacts_frame, orient="vertical", command=self.artifacts_tree.yview)
        artifact_scroll.pack(side="right", fill="y")
        self.artifacts_tree.configure(yscrollcommand=artifact_scroll.set)

    def _build_metrics_tab(self, parent: ttk.Frame) -> None:
        upper = ttk.LabelFrame(parent, text="Overall Evaluation", style="Section.TLabelframe", padding=10)
        upper.pack(fill="x", pady=(0, 10))
        self.metrics_overall_text = tk.Text(upper, wrap="word", font=("Consolas", 10), height=8, bg="#FCFBF8")
        self.metrics_overall_text.pack(fill="x")

        lower = ttk.LabelFrame(parent, text="Per-Gene Scores", style="Section.TLabelframe", padding=10)
        lower.pack(fill="both", expand=True)

        columns = ("gene", "queries", "aligned", "identity", "coverage", "mapq", "score")
        self.metrics_tree = ttk.Treeview(lower, columns=columns, show="headings")
        headings = {
            "gene": "Gene",
            "queries": "Queries",
            "aligned": "Aligned",
            "identity": "Mean Identity",
            "coverage": "Mean Coverage",
            "mapq": "Mean MAPQ",
            "score": "Mean Score",
        }
        widths = {"gene": 90, "queries": 80, "aligned": 80, "identity": 110, "coverage": 110, "mapq": 95, "score": 120}
        for key in columns:
            self.metrics_tree.heading(key, text=headings[key])
            self.metrics_tree.column(key, width=widths[key], anchor="center")
        self.metrics_tree.pack(side="left", fill="both", expand=True)
        tree_scroll = ttk.Scrollbar(lower, orient="vertical", command=self.metrics_tree.yview)
        tree_scroll.pack(side="right", fill="y")
        self.metrics_tree.configure(yscrollcommand=tree_scroll.set)

    def _build_visuals_tab(self, parent: ttk.Frame) -> None:
        paned = ttk.Panedwindow(parent, orient="horizontal")
        paned.pack(fill="both", expand=True)

        left = ttk.LabelFrame(paned, text="Available Graphs", style="Section.TLabelframe", padding=10)
        right = ttk.LabelFrame(paned, text="Preview", style="Section.TLabelframe", padding=10)
        paned.add(left, weight=2)
        paned.add(right, weight=5)

        self.images_listbox = tk.Listbox(left, exportselection=False)
        self.images_listbox.pack(side="left", fill="both", expand=True)
        self.images_listbox.bind("<<ListboxSelect>>", self._on_image_select)
        image_scroll = ttk.Scrollbar(left, orient="vertical", command=self.images_listbox.yview)
        image_scroll.pack(side="right", fill="y")
        self.images_listbox.configure(yscrollcommand=image_scroll.set)

        controls = ttk.Frame(right)
        controls.pack(fill="x", pady=(0, 8))
        ttk.Button(controls, text="Open Image", command=self._open_selected_image).pack(side="left")

        self.image_label = ttk.Label(right, text="Select a generated PNG to preview it here.", anchor="center", justify="center")
        self.image_label.pack(fill="both", expand=True)

    def _build_files_tab(self, parent: ttk.Frame) -> None:
        paned = ttk.Panedwindow(parent, orient="horizontal")
        paned.pack(fill="both", expand=True)

        left = ttk.LabelFrame(paned, text="Generated Files", style="Section.TLabelframe", padding=10)
        right = ttk.LabelFrame(paned, text="File Viewer", style="Section.TLabelframe", padding=10)
        paned.add(left, weight=2)
        paned.add(right, weight=5)

        self.files_listbox = tk.Listbox(left, exportselection=False)
        self.files_listbox.pack(side="left", fill="both", expand=True)
        self.files_listbox.bind("<<ListboxSelect>>", self._on_file_select)
        files_scroll = ttk.Scrollbar(left, orient="vertical", command=self.files_listbox.yview)
        files_scroll.pack(side="right", fill="y")
        self.files_listbox.configure(yscrollcommand=files_scroll.set)

        actions = ttk.Frame(right)
        actions.pack(fill="x", pady=(0, 8))
        ttk.Button(actions, text="Open File", command=self._open_selected_file).pack(side="left")

        self.file_viewer = tk.Text(right, wrap="word", font=("Consolas", 10), bg="#FCFBF8")
        self.file_viewer.pack(fill="both", expand=True)

    def _build_logs_tab(self, parent: ttk.Frame) -> None:
        log_frame = ttk.LabelFrame(parent, text="Live Console Output", style="Section.TLabelframe", padding=10)
        log_frame.pack(fill="both", expand=True)
        self.log_text = tk.Text(log_frame, wrap="word", font=("Consolas", 10), bg="#FAFAFA")
        self.log_text.pack(side="left", fill="both", expand=True)
        scrollbar = ttk.Scrollbar(log_frame, orient="vertical", command=self.log_text.yview)
        scrollbar.pack(side="right", fill="y")
        self.log_text.configure(yscrollcommand=scrollbar.set)

    def _labeled_entry(self, parent: ttk.Widget, label: str, variable: tk.StringVar, row: int) -> None:
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=(0, 10), pady=6)
        ttk.Entry(parent, textvariable=variable).grid(row=row, column=1, sticky="ew", pady=6)
        parent.columnconfigure(1, weight=1)

    def _labeled_file_row(self, parent: ttk.Widget, label: str, variable: tk.StringVar, row: int, filetypes: list[tuple[str, str]]) -> None:
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=(0, 10), pady=6)
        row_frame = ttk.Frame(parent)
        row_frame.grid(row=row, column=1, sticky="ew", pady=6)
        ttk.Entry(row_frame, textvariable=variable).pack(side="left", fill="x", expand=True)
        ttk.Button(row_frame, text="Browse", command=lambda: self._pick_file(variable, filetypes)).pack(side="left", padx=(8, 0))
        parent.columnconfigure(1, weight=1)

    def _labeled_dir_row(self, parent: ttk.Widget, label: str, variable: tk.StringVar, row: int) -> None:
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=(0, 10), pady=6)
        row_frame = ttk.Frame(parent)
        row_frame.grid(row=row, column=1, sticky="ew", pady=6)
        ttk.Entry(row_frame, textvariable=variable).pack(side="left", fill="x", expand=True)
        ttk.Button(row_frame, text="Browse", command=lambda: self._pick_directory(variable)).pack(side="left", padx=(8, 0))
        parent.columnconfigure(1, weight=1)

    def _pick_file(self, variable: tk.StringVar, filetypes: list[tuple[str, str]]) -> None:
        selected = filedialog.askopenfilename(initialdir=str(ROOT), filetypes=filetypes)
        if selected:
            variable.set(selected)

    def _pick_directory(self, variable: tk.StringVar) -> None:
        selected = filedialog.askdirectory(initialdir=str(ROOT))
        if selected:
            variable.set(selected)

    def _toggle_query_source(self) -> None:
        use_sequence = self.query_source.get() == "sequence"
        self.query_sequence_text.configure(state="normal" if use_sequence else "disabled")
        self.query_fasta_entry.configure(state="disabled" if use_sequence else "normal")

    def _append_log(self, text: str) -> None:
        self.log_text.insert("end", text)
        self.log_text.see("end")

    def _set_text_widget(self, widget: tk.Text, content: str) -> None:
        widget.configure(state="normal")
        widget.delete("1.0", "end")
        widget.insert("1.0", content)
        widget.configure(state="disabled")

    def _build_common_args(self) -> list[str]:
        args = []
        if self.threads_var.get().strip():
            args.extend(["--threads", self.threads_var.get().strip()])
        if self.minigraph_var.get().strip():
            args.extend(["--minigraph-bin", self.minigraph_var.get().strip()])
        if self.panaligner_var.get().strip():
            args.extend(["--panaligner-bin", self.panaligner_var.get().strip()])
        return args

    def run_pipeline_mode(self) -> None:
        command = [sys.executable, str(MAIN_SCRIPT), self.pipeline_mode.get()]
        command.extend(self._build_common_args())

        if self.split_seed_var.get().strip():
            command.extend(["--split-seed", self.split_seed_var.get().strip()])
        if self.test_fraction_var.get().strip():
            command.extend(["--test-fraction", self.test_fraction_var.get().strip()])
        if self.theory_graph_var.get().strip():
            command.extend(["--theory-graph", self.theory_graph_var.get().strip()])

        fasta_paths = [part.strip() for part in self.input_fastas_var.get().split(";") if part.strip()]
        if fasta_paths:
            command.append("--input-fastas")
            command.extend(fasta_paths)

        self._start_process(command, "Running project mode")

    def run_custom_query(self) -> None:
        command = [sys.executable, str(MAIN_SCRIPT), "--custom-query-analysis"]
        command.extend(self._build_common_args())

        graph_manifest = self.graph_manifest_var.get().strip()
        if graph_manifest:
            command.extend(["--graph-manifest", graph_manifest])

        query_argument = self.query_argument_var.get().strip()
        if query_argument:
            command.extend(["--query-argument", query_argument])

        gene_value = self.query_gene_var.get().strip()
        if gene_value:
            command.extend(["--gene", gene_value])

        output_dir = self.query_output_var.get().strip()
        if output_dir:
            command.extend(["--output-dir", output_dir])

        if self.query_source.get() == "sequence":
            sequence = self.query_sequence_text.get("1.0", "end").strip()
            if not sequence:
                messagebox.showerror("Missing sequence", "Enter a DNA query sequence before running the analysis.")
                return
            command.extend(["--query-sequence", sequence])
        else:
            query_fasta = self.query_fasta_var.get().strip()
            if not query_fasta:
                messagebox.showerror("Missing FASTA", "Choose a query FASTA file before running the analysis.")
                return
            command.extend(["--query-fasta", query_fasta])

        self._start_process(command, "Running custom query analysis")

    def _start_process(self, command: list[str], label: str) -> None:
        if self.process is not None:
            messagebox.showinfo("Run in progress", "Wait for the current run to finish or stop it before starting another one.")
            return

        self.log_text.delete("1.0", "end")
        self.output_buffer.clear()
        self.last_summary = None
        self.preview_image = None
        self.image_label.configure(image="", text="Run started. Waiting for generated visuals...")
        self._set_text_widget(self.summary_text, "Working...\n")
        self.status_var.set(label)

        self._append_log("Command:\n")
        self._append_log(" ".join(f'"{part}"' if " " in part else part for part in command) + "\n\n")

        def worker() -> None:
            try:
                self.process = subprocess.Popen(
                    command,
                    cwd=str(ROOT),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                )
                assert self.process.stdout is not None
                for line in self.process.stdout:
                    self.output_queue.put(("line", line))
                return_code = self.process.wait()
                self.output_queue.put(("done", str(return_code)))
            except Exception as exc:  # pragma: no cover
                self.output_queue.put(("error", str(exc)))

        threading.Thread(target=worker, daemon=True).start()

    def stop_run(self) -> None:
        if self.process is None:
            return
        self.process.terminate()
        self.status_var.set("Stopping current run...")

    def _poll_output(self) -> None:
        try:
            while True:
                kind, payload = self.output_queue.get_nowait()
                if kind == "line":
                    self.output_buffer.append(payload)
                    self._append_log(payload)
                elif kind == "done":
                    code = int(payload)
                    self.process = None
                    self._handle_completion(code)
                elif kind == "error":
                    self.process = None
                    self.status_var.set("Run failed before launch.")
                    messagebox.showerror("Launch error", payload)
        except queue.Empty:
            pass
        self.root.after(150, self._poll_output)

    def _handle_completion(self, return_code: int) -> None:
        output_text = "".join(self.output_buffer).strip()
        self.status_var.set("Run completed successfully." if return_code == 0 else f"Run failed with exit code {return_code}.")
        self.last_summary = self._extract_json_summary(output_text)
        self._refresh_outputs()
        self._render_summary(self.last_summary, return_code, output_text)
        self._render_preview(self.last_summary)

    def _extract_json_summary(self, output_text: str) -> dict | None:
        stripped = output_text.strip()
        if not stripped:
            return None
        start = stripped.find("{")
        end = stripped.rfind("}")
        if start == -1 or end == -1 or end <= start:
            return None
        try:
            return json.loads(stripped[start : end + 1])
        except json.JSONDecodeError:
            return None

    def _render_summary(self, summary: dict | None, return_code: int, output_text: str) -> None:
        lines: list[str] = []
        if summary:
            lines.append("Latest Run Summary")
            lines.append("===================")
            for key, value in summary.items():
                lines.append(f"{key.replace('_', ' ').title()}: {value}")
            lines.append("")

        metrics = self._load_json(OUTPUTS_DIR / "evaluation" / "alignment_metrics.json")
        if isinstance(metrics, dict) and "overall" in metrics:
            overall = metrics["overall"]
            lines.append("Current Evaluation Snapshot")
            lines.append("===========================")
            lines.append(f"Queries: {overall.get('query_count', 0)}")
            lines.append(f"Aligned: {overall.get('aligned_query_count', 0)}")
            lines.append(f"Alignment rate: {overall.get('alignment_rate', 0.0):.4f}")
            lines.append(f"Mean identity: {overall.get('mean_identity', 0.0):.4f}")
            lines.append(f"Mean coverage: {overall.get('mean_coverage', 0.0):.4f}")
            lines.append(f"Mean MAPQ: {overall.get('mean_mapq', 0.0):.4f}")
            lines.append(f"Mean score: {overall.get('mean_alignment_score', 0.0):.4f}")
            lines.append("")

        prediction_path = Path(summary["prediction_json"]) if summary and "prediction_json" in summary else DEFAULT_CUSTOM_OUTPUT / "prediction.json"
        prediction = self._load_json(prediction_path)
        if isinstance(prediction, dict):
            lines.append("Custom Query Snapshot")
            lines.append("=====================")
            lines.append(f"Selected gene: {prediction.get('selected_gene', '')}")
            lines.append(f"Prediction: {prediction.get('prediction', '')}")
            lines.append(f"Confidence: {prediction.get('confidence', 0.0):.4f}")
            lines.append(f"Healthy score: {prediction.get('healthy_score', 0.0):.4f}")
            lines.append(f"Unhealthy score: {prediction.get('unhealthy_score', 0.0):.4f}")
            lines.append(f"Alignment detected: {prediction.get('alignment_detected', False)}")
            lines.append(f"Argument: {prediction.get('query_argument', '')}")
            explanation = prediction.get("explanation")
            if explanation:
                lines.append("")
                lines.append(explanation)

        if not lines:
            if return_code == 0:
                lines = ["Run completed. Refresh outputs or open the generated files to inspect details."]
            else:
                trimmed = output_text[-800:] if output_text else "No console output captured."
                lines = ["Run failed.", "", trimmed]

        self._set_text_widget(self.summary_text, "\n".join(lines) + "\n")

    def _render_preview(self, summary: dict | None) -> None:
        self.preview_image = None
        preferred_paths: list[Path] = []
        if summary and isinstance(summary.get("score_plot"), str):
            preferred_paths.append(Path(summary["score_plot"]))
        preferred_paths.extend(self.image_paths[:1])
        selected = next((path for path in preferred_paths if path.exists()), None)
        if not selected:
            self.image_label.configure(image="", text="No graph preview available yet.")
            return
        self._show_image(selected)

    def _show_image(self, image_path: Path) -> None:
        try:
            image = tk.PhotoImage(file=str(image_path))
        except tk.TclError:
            self.image_label.configure(image="", text=f"Unable to preview {image_path.name} in Tk.")
            self.preview_image = None
            return
        scale = max(1, image.width() // 760, image.height() // 540)
        self.preview_image = image.subsample(scale, scale) if scale > 1 else image
        self.image_label.configure(image=self.preview_image, text="")

    def _refresh_outputs(self) -> None:
        self._refresh_metrics()
        self._refresh_images()
        self._refresh_files()
        self._refresh_artifacts()

    def _refresh_metrics(self) -> None:
        metrics = self._load_json(OUTPUTS_DIR / "evaluation" / "alignment_metrics.json")
        if not isinstance(metrics, dict):
            self._set_text_widget(self.metrics_overall_text, "No evaluation metrics found yet.\nRun evaluation or full pipeline to populate this panel.\n")
            self._clear_tree(self.metrics_tree)
            return

        overall = metrics.get("overall", {})
        overall_lines = [
            "Overall Evaluation Metrics",
            "==========================",
            f"Query count: {overall.get('query_count', 0)}",
            f"Aligned query count: {overall.get('aligned_query_count', 0)}",
            f"Alignment rate: {overall.get('alignment_rate', 0.0):.4f}",
            f"Mean identity: {overall.get('mean_identity', 0.0):.4f}",
            f"Mean coverage: {overall.get('mean_coverage', 0.0):.4f}",
            f"Mean MAPQ: {overall.get('mean_mapq', 0.0):.4f}",
            f"Mean alignment score: {overall.get('mean_alignment_score', 0.0):.4f}",
        ]
        self._set_text_widget(self.metrics_overall_text, "\n".join(overall_lines) + "\n")

        self._clear_tree(self.metrics_tree)
        for gene, summary in sorted(metrics.get("per_gene", {}).items()):
            self.metrics_tree.insert(
                "",
                "end",
                values=(
                    gene,
                    summary.get("sequence_count", 0),
                    summary.get("aligned_sequences", 0),
                    f"{summary.get('mean_identity', 0.0):.4f}",
                    f"{summary.get('mean_coverage', 0.0):.4f}",
                    f"{summary.get('mean_mapq', 0.0):.4f}",
                    f"{summary.get('mean_alignment_score', 0.0):.4f}",
                ),
            )

    def _refresh_images(self) -> None:
        self.image_paths = sorted(OUTPUTS_DIR.rglob("*.png"))
        self.images_listbox.delete(0, "end")
        for path in self.image_paths:
            self.images_listbox.insert("end", str(path.relative_to(ROOT)))
        if self.image_paths:
            self.images_listbox.selection_set(0)
            self.images_listbox.event_generate("<<ListboxSelect>>")
        else:
            self.image_label.configure(image="", text="No generated graphs found yet.")

    def _refresh_files(self) -> None:
        wanted_suffixes = {".txt", ".json", ".csv", ".gaf"}
        self.file_paths = sorted(path for path in OUTPUTS_DIR.rglob("*") if path.is_file() and path.suffix.lower() in wanted_suffixes)
        self.files_listbox.delete(0, "end")
        for path in self.file_paths:
            self.files_listbox.insert("end", str(path.relative_to(ROOT)))
        if self.file_paths:
            self.files_listbox.selection_set(0)
            self.files_listbox.event_generate("<<ListboxSelect>>")
        else:
            self._set_text_widget(self.file_viewer, "No generated text or JSON outputs found yet.\n")

    def _refresh_artifacts(self) -> None:
        self._clear_tree(self.artifacts_tree)
        artifact_paths: list[Path] = []
        if self.last_summary:
            for value in self.last_summary.values():
                if isinstance(value, str):
                    candidate = Path(value)
                    if candidate.exists():
                        artifact_paths.append(candidate)
        for path in artifact_paths:
            kind = path.suffix.lower().lstrip(".") or "file"
            self.artifacts_tree.insert("", "end", values=(str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path), kind))

    def _clear_tree(self, tree: ttk.Treeview) -> None:
        for item in tree.get_children():
            tree.delete(item)

    def _load_json(self, path: Path) -> dict | list | None:
        if not path.exists():
            return None
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return None

    def _on_image_select(self, _event: object) -> None:
        selection = self.images_listbox.curselection()
        if not selection:
            return
        image_path = self.image_paths[selection[0]]
        self._show_image(image_path)

    def _on_file_select(self, _event: object) -> None:
        selection = self.files_listbox.curselection()
        if not selection:
            return
        file_path = self.file_paths[selection[0]]
        try:
            content = file_path.read_text(encoding="utf-8")
        except OSError as exc:
            content = f"Unable to read file:\n{exc}\n"
        self._set_text_widget(self.file_viewer, content)

    def _open_primary_output(self) -> None:
        preferred: list[Path] = []
        if self.last_summary:
            for key in ("score_plot", "prediction_json", "final_report", "theory_report"):
                value = self.last_summary.get(key)
                if isinstance(value, str):
                    preferred.append(Path(value))
        preferred.extend([OUTPUTS_DIR / "reports" / "final_project_summary.txt", OUTPUTS_DIR / "evaluation" / "alignment_metrics.json"])
        target = next((path for path in preferred if path.exists()), None)
        if target is None:
            messagebox.showinfo("No output yet", "Run something first so there is an output artifact to open.")
            return
        open_path(target)

    def _open_output_dir(self) -> None:
        output_dir = Path(self.query_output_var.get().strip() or DEFAULT_CUSTOM_OUTPUT)
        output_dir.mkdir(parents=True, exist_ok=True)
        open_path(output_dir)

    def _open_selected_image(self) -> None:
        selection = self.images_listbox.curselection()
        if not selection:
            return
        open_path(self.image_paths[selection[0]])

    def _open_selected_file(self) -> None:
        selection = self.files_listbox.curselection()
        if not selection:
            return
        open_path(self.file_paths[selection[0]])

    def _open_selected_artifact(self) -> None:
        selection = self.artifacts_tree.selection()
        if not selection:
            return
        artifact_path = Path(self.artifacts_tree.item(selection[0], "values")[0])
        if not artifact_path.is_absolute():
            artifact_path = ROOT / artifact_path
        if artifact_path.exists():
            open_path(artifact_path)


def main() -> None:
    root = tk.Tk()
    PanAlignerGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
