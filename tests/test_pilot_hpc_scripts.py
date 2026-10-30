"""Static checks on the pilot's PBS jobs. Their real behaviour is only known after running on the cluster."""

from __future__ import annotations

import re
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


if __name__ == "__main__":
    unittest.main()
