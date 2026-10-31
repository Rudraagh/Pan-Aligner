"""Static checks on the pilot's PBS jobs. Their real behaviour is only known after running on the cluster."""

from __future__ import annotations

import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
JOBS = sorted((ROOT / "scripts" / "hpc").glob("pbs_pilot_*.sh"))


class PilotJobs(unittest.TestCase):
    def test_all_stages_exist(self):
        self.assertEqual([re.search(r"pilot_(\d)", p.name).group(1) for p in JOBS], list("012345"))

    def test_jobs_declare_resources_and_no_gpu(self):
        for path in JOBS:
            text = path.read_text(encoding="utf-8")
            self.assertTrue(text.startswith("#!/bin/bash"), path.name)
            self.assertIn("#PBS -q cpuq", text, path.name)
            self.assertRegex(text, r"#PBS -l select=1:ncpus=\d+:mem=\d+gb", path.name)
            self.assertRegex(text, r"#PBS -l walltime=", path.name)
            self.assertIn("UNVERIFIED", text, path.name)
            self.assertIn("pilot_finish", text, path.name)
            self.assertNotIn("gpuq", text, path.name)

    def test_scripts_use_unix_line_endings(self):
        paths = [*JOBS, ROOT / "scripts" / "hpc" / "pilot_common.sh", ROOT / "scripts" / "ml" / "wg_truth" / "synthetic_smoke.sh"]
        for path in paths:
            self.assertNotIn(b"\r", path.read_bytes(), path.name)

    def test_downloads_only_with_explicit_flag(self):
        for path in JOBS:
            for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
                if "curl " in line and not line.lstrip().startswith("#"):
                    self.assertIn("1_inputs", path.name, f"{path.name}:{number}")
        self.assertIn('if [[ "$FETCH" == "1" ]]', (ROOT / "scripts" / "hpc" / "pbs_pilot_1_inputs.sh").read_text(encoding="utf-8"))

    def test_bundle_includes_third_party_sources(self):
        text = (ROOT / "scripts" / "hpc" / "make_staging_bundle.sh").read_text(encoding="utf-8")
        self.assertIn("tools/third_party/pbsim2", text)
        self.assertIn("tools/third_party/minimap2", text)

    def test_bundle_export_cannot_convert_line_endings(self):
        text = (ROOT / "scripts" / "hpc" / "make_staging_bundle.sh").read_text(encoding="utf-8")
        archives = [line for line in text.splitlines() if " archive --format" in line]
        self.assertEqual(len(archives), 2)
        for line in archives:
            self.assertIn("core.autocrlf=false", line)

    @unittest.skipUnless((ROOT / ".git").exists(), "needs the Git index")
    def test_k8_is_executable_in_git(self):
        relative = "tools/third_party/k8-bin/k8-1.2/k8-x86_64-Linux"
        entry = subprocess.run(["git", "ls-files", "-s", relative], cwd=ROOT, capture_output=True, text=True, check=True).stdout
        self.assertTrue(entry.startswith("100755 "), entry)

    def test_build_makes_k8_executable_and_runs_it(self):
        # Bundles built on Windows lose the mode when re-archived, so the node has to set it.
        text = (ROOT / "scripts" / "hpc" / "build_tools_linux.sh").read_text(encoding="utf-8")
        self.assertIn('chmod u+x "$K8"', text)
        self.assertIn('"$K8" "$PAFTOOLS" version', text)


@unittest.skipUnless(os.name == "posix", "runs the job under bash")
class AlignJobReads(unittest.TestCase):
    def run_job(self, pilot: Path, bin_dir: Path) -> subprocess.CompletedProcess:
        env = {**os.environ, "PROJECT_DIR": str(ROOT), "PILOT_DIR": str(pilot), "BIN_DIR": str(bin_dir), "PYTHON_BIN": sys.executable, "ALIGNER": "minigraph", "NCPUS": "1"}
        return subprocess.run(["bash", str(ROOT / "scripts" / "hpc" / "pbs_pilot_4_align.sh")], env=env, capture_output=True, text=True)

    def setup_tools(self, root: Path) -> Path:
        bin_dir = root / "bin"
        bin_dir.mkdir()
        aligner = bin_dir / "minigraph"
        aligner.write_text('#!/bin/sh\nprintf "read\\t1\\n"\n', encoding="ascii", newline="\n")
        aligner.chmod(0o755)
        return bin_dir

    def test_missing_reads_fail_and_leave_no_gaf(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            setting = root / "pilot" / "sim" / "acc95_len10000"
            setting.mkdir(parents=True)
            (setting / "mat.reads.fa").write_text(">a\nACGT\n", encoding="ascii")
            result = self.run_job(root / "pilot", self.setup_tools(root))
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertFalse((setting / "minigraph.gaf").exists())
            self.assertFalse((setting / "all.reads.fa").exists())

    def test_reads_from_both_haplotypes_are_aligned(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            setting = root / "pilot" / "sim" / "acc95_len10000"
            setting.mkdir(parents=True)
            (setting / "mat.reads.fa").write_text(">a\nACGT\n", encoding="ascii")
            (setting / "pat.reads.fa").write_text(">b\nACGT\n", encoding="ascii")
            result = self.run_job(root / "pilot", self.setup_tools(root))
            self.assertEqual(result.returncode, 0, result.stdout)
            self.assertTrue((setting / "minigraph.gaf").exists())
            self.assertEqual((setting / "all.reads.fa").read_text(encoding="ascii"), ">a\nACGT\n>b\nACGT\n")


if __name__ == "__main__":
    unittest.main()
