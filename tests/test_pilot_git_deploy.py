"""The pilot must work from a plain Git checkout. Static checks; real behaviour is only known on the cluster."""

from __future__ import annotations

import configparser
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HPC = ROOT / "scripts" / "hpc"
RESOURCES = {
    "pbs_pilot_0_build_smoke.sh": ("4", "8gb", "01:00:00"),
    "pbs_pilot_1_inputs.sh": ("16", "32gb", "04:00:00"),
    "pbs_pilot_2_loadtest.sh": ("32", "64gb", "01:00:00"),
    "pbs_pilot_3_simulate.sh": ("4", "8gb", "02:00:00"),
    "pbs_pilot_4_align.sh": ("32", "64gb", "08:00:00"),
    "pbs_pilot_5_label.sh": ("4", "32gb", "04:00:00"),
}
# The commits inside the canonical bundle panaligner_c018c46b339f.tar.gz.
TOOL_COMMITS = {
    "PanAligner": "bdc00cfaf13578597e76583a47a58eeefde509e6",
    "minigraph": "2f569ebe3071fcb242f1bf0eabbb917941b47239",
    "zlib": "f9dd6009be3ed32415edf1e89d1bc38380ecb95d",
    "tools/third_party/pbsim2": "eeb5a19420534a0f672c81db2670117e62a9ee38",
    "tools/third_party/minimap2": "e5066c79765f78fad9a288331f860d78670cc45f",
}


def git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=True).stdout


class Resources(unittest.TestCase):
    def test_every_job_keeps_its_resources(self):
        for name, (cpus, memory, walltime) in RESOURCES.items():
            text = (HPC / name).read_text(encoding="utf-8")
            self.assertIn(f"#PBS -l select=1:ncpus={cpus}:mem={memory}\n", text, name)
            self.assertIn(f"#PBS -l walltime={walltime}\n", text, name)

    def test_no_job_requests_a_gpu(self):
        for name in RESOURCES:
            text = (HPC / name).read_text(encoding="utf-8")
            self.assertNotRegex(text, r"ngpus|gpuq|cupy", name)


@unittest.skipUnless((ROOT / ".git").exists(), "needs Git metadata")
class ToolPins(unittest.TestCase):
    def test_gitmodules_maps_each_tool(self):
        config = configparser.ConfigParser()
        config.read(ROOT / ".gitmodules", encoding="utf-8")
        paths = {config[section]["path"]: config[section]["url"] for section in config.sections()}
        self.assertEqual(set(paths), set(TOOL_COMMITS))
        for url in paths.values():
            self.assertTrue(url.startswith("https://"), url)

    def test_tree_pins_the_bundle_commits(self):
        for path, commit in TOOL_COMMITS.items():
            entry = git("ls-files", "-s", "--", path).split()
            self.assertEqual(entry[:2], ["160000", commit], path)

    def test_tools_are_not_ignored_and_run_data_is(self):
        text = (ROOT / ".gitignore").read_text(encoding="utf-8")
        for name in ("PanAligner/", "minigraph/", "zlib/"):
            self.assertNotIn(f"\n{name}\n", text)
        self.assertIn("/pilot_data/", text)

    def test_manifest_matches_the_index(self):
        subprocess.run([sys.executable, "-B", str(HPC / "make_pilot_manifest.py"), "--check"], cwd=ROOT, check=True, capture_output=True)


class CheckoutIdentity(unittest.TestCase):
    def test_runs_record_the_git_commit(self):
        text = (HPC / "pilot_common.sh").read_text(encoding="utf-8")
        self.assertIn('git -C "$PROJECT_DIR" rev-parse HEAD', text)

    def test_wheels_for_python_311_are_listed_and_present(self):
        folder = ROOT / "wheelhouse" / "py311"
        wheels = (folder / "WHEELS.txt").read_text(encoding="utf-8").split()
        self.assertGreater(len(wheels), 0)
        for name in wheels:
            self.assertTrue((folder / name).is_file(), name)


class Interpreter(unittest.TestCase):
    def test_default_is_the_named_python_only(self):
        self.assertIn("/apps/anaconda3/envs/gcc", (HPC / "hpc_env.sh").read_text(encoding="utf-8"))
        setup = (HPC / "pbs_setup_venv.sh").read_text(encoding="utf-8")
        self.assertNotIn("for candidate in python3", setup)
        self.assertIn("HPC_PYTHON_VERSION", setup)
        self.assertNotIn('PYTHON_BIN="${PYTHON_BIN:-python3}"', (HPC / "pilot_common.sh").read_text(encoding="utf-8"))

    def test_nothing_installs_outside_the_project(self):
        for path in HPC.glob("*.sh"):
            code = "\n".join(line for line in path.read_text(encoding="utf-8").splitlines() if not line.lstrip().startswith("#"))
            self.assertNotRegex(code, r"conda (install|create)|pip install --user|sudo ", path.name)


class Preflight(unittest.TestCase):
    PATH = HPC / "pilot_preflight.sh"

    def test_it_never_submits_or_downloads(self):
        lines = [line for line in self.PATH.read_text(encoding="utf-8").splitlines() if not line.lstrip().startswith("#")]
        code = "\n".join(lines)
        self.assertNotRegex(code, r"\bcurl\b.*https?:|\bwget\b|\bgit (fetch|pull|clone)\b|pip install")
        self.assertNotRegex(code, r"(^|[;&|(]\s*)qsub\s")

    def test_guide_limits_authorisation_to_job_0(self):
        text = (ROOT / "docs" / "18_hpc_git_deployment.md").read_text(encoding="utf-8")
        self.assertIn("Jobs 1–5 are not authorised by completing Job 0", text)
        self.assertEqual(re.findall(r"^qsub .*", text, flags=re.M), [line for line in re.findall(r"^qsub .*", text, flags=re.M) if "pbs_pilot_0_build_smoke.sh" in line])

    @unittest.skipUnless(os.name == "posix" and (ROOT / ".git").exists(), "runs under bash in a checkout")
    def test_wrong_commit_fails_clearly(self):
        result = subprocess.run(["bash", str(self.PATH), "--expected-commit", "0" * 40, "--before-venv", "--no-pbs"], cwd=ROOT, capture_output=True, text=True, env={**os.environ, "PYTHON_BIN": sys.executable})
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("FAIL", result.stdout)
        self.assertIn("expected commit", result.stdout)

    @unittest.skipUnless(os.name == "posix", "runs under bash")
    def test_commit_argument_is_required(self):
        result = subprocess.run(["bash", str(self.PATH)], cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)


if __name__ == "__main__":
    unittest.main()
