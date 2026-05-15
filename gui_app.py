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
DEFAULT_GRAPH_MANIFEST = ROOT / "data" / "metadata" / "graph_manifest.json"
DEFAULT_CUSTOM_OUTPUT = ROOT / "outputs" / "alignments" / "custom_query"


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
        self.root.geometry("1280x860")
        self.root.minsize(1100, 760)

        self.process: subprocess.Popen[str] | None = None
        self.output_queue: queue.Queue[tuple[str, str]] = queue.Queue()
        self.output_buffer: list[str] = []
        self.preview_image: tk.PhotoImage | None = None
        self.last_summary: dict | None = None

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
        self.summary_var = tk.StringVar(value="Run a pipeline mode or custom query analysis to view a summary here.")

        self._configure_style()
        self._build_layout()
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
            text="Run the core project pipeline, execute custom query alignment checks, and review outputs from one place.",
            background="#F5EFE6",
            foreground="#5A4A42",
        ).pack(anchor="w", pady=(4, 0))

        body = ttk.Panedwindow(self.root, orient="horizontal")
        body.pack(fill="both", expand=True, padx=16, pady=(0, 16))

        controls = ttk.Frame(body, padding=8)
        results = ttk.Frame(body, padding=8)
        body.add(controls, weight=5)
        body.add(results, weight=6)

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

        summary_frame = ttk.LabelFrame(parent, text="Run Summary", style="Section.TLabelframe", padding=12)
        summary_frame.pack(fill="x", pady=(0, 10))
        ttk.Label(summary_frame, textvariable=self.summary_var, wraplength=520, justify="left").pack(anchor="w")

        preview_actions = ttk.Frame(summary_frame)
        preview_actions.pack(fill="x", pady=(8, 0))
        ttk.Button(preview_actions, text="Open Referenced Output", command=self._open_primary_output).pack(side="left")

        preview_frame = ttk.LabelFrame(parent, text="Plot Preview", style="Section.TLabelframe", padding=12)
        preview_frame.pack(fill="both", expand=False, pady=(0, 10))
        self.preview_label = ttk.Label(preview_frame, text="Custom query score plots will appear here after a successful run.", anchor="center", justify="center")
        self.preview_label.pack(fill="both", expand=True)

        log_frame = ttk.LabelFrame(parent, text="Live Console Output", style="Section.TLabelframe", padding=10)
        log_frame.pack(fill="both", expand=True)
        self.log_text = tk.Text(log_frame, height=18, wrap="word", font=("Consolas", 10), bg="#FAFAFA")
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
        text_state = "normal" if use_sequence else "disabled"
        entry_state = "disabled" if use_sequence else "normal"
        self.query_sequence_text.configure(state=text_state)
        self.query_fasta_entry.configure(state=entry_state)

    def _append_log(self, text: str) -> None:
        self.log_text.insert("end", text)
        self.log_text.see("end")

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
        self.preview_label.configure(image="", text="Run started. Waiting for generated outputs...")
        self.summary_var.set("Working...")
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
            except Exception as exc:  # pragma: no cover - GUI fallback path
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
        if return_code == 0:
            self.status_var.set("Run completed successfully.")
        else:
            self.status_var.set(f"Run failed with exit code {return_code}.")

        summary = self._extract_json_summary(output_text)
        self.last_summary = summary
        self._render_summary(summary, return_code, output_text)
        self._render_preview(summary)

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
        if summary:
            lines = [f"{key.replace('_', ' ').title()}: {value}" for key, value in summary.items()]
            self.summary_var.set("\n".join(lines))
            return
        if return_code == 0:
            self.summary_var.set("Run completed. No JSON summary was returned, so review the console output and generated files.")
            return
        trimmed = output_text[-600:] if output_text else "No console output captured."
        self.summary_var.set(f"Run failed.\n\n{trimmed}")

    def _render_preview(self, summary: dict | None) -> None:
        self.preview_image = None
        plot_path = None
        if summary and isinstance(summary.get("score_plot"), str):
            plot_path = Path(summary["score_plot"])
        if not plot_path or not plot_path.exists():
            self.preview_label.configure(image="", text="No plot preview available for this run.")
            return

        image = tk.PhotoImage(file=str(plot_path))
        scale = max(1, image.width() // 520, image.height() // 320)
        self.preview_image = image.subsample(scale, scale) if scale > 1 else image
        self.preview_label.configure(image=self.preview_image, text="")

    def _open_primary_output(self) -> None:
        if not self.last_summary:
            messagebox.showinfo("No output yet", "Run something first so the GUI knows which result to open.")
            return
        for key in ("score_plot", "prediction_json", "final_report", "theory_report"):
            value = self.last_summary.get(key)
            if isinstance(value, str) and value:
                target = Path(value)
                if target.exists():
                    open_path(target)
                    return
        messagebox.showinfo("No file output", "The latest run did not expose a directly openable result path.")

    def _open_output_dir(self) -> None:
        output_dir = Path(self.query_output_var.get().strip() or DEFAULT_CUSTOM_OUTPUT)
        output_dir.mkdir(parents=True, exist_ok=True)
        open_path(output_dir)


def main() -> None:
    root = tk.Tk()
    app = PanAlignerGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
