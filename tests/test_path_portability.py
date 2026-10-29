import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "tools"))

from common import project_relative, read_json, resolve_project_path
from predictor import resolve_graph_path
from wsl_exec import WORKSPACE_ROOT, to_wsl_arg


class PathPortabilityTests(unittest.TestCase):
    def test_manifest_graph_path_resolves_from_project_root(self):
        self.assertEqual(resolve_graph_path("graphs/psen2.healthy.gfa"), ROOT / "graphs" / "psen2.healthy.gfa")

    def test_legacy_absolute_paths_from_other_machines_remap_into_project(self):
        expected = (ROOT / "data" / "raw" / "app_reference.fa").resolve()
        windows_path = r"C:\Users\someone\Elsewhere\Project Copy\data\raw\app_reference.fa"
        posix_path = "/mnt/c/Users/someone/Project Copy/data/raw/app_reference.fa"
        self.assertEqual(resolve_project_path(windows_path), expected)
        self.assertEqual(resolve_project_path(posix_path), expected)

    def test_project_relative_is_posix_and_round_trips(self):
        path = ROOT / "data" / "raw" / "app_reference.fa"
        self.assertEqual(project_relative(path), "data/raw/app_reference.fa")
        self.assertEqual(resolve_project_path(project_relative(path)), path.resolve())

    def test_committed_manifests_contain_no_machine_specific_paths(self):
        for manifest in sorted((ROOT / "data" / "metadata").glob("*.json")):
            with self.subTest(manifest=manifest.name):
                text = manifest.read_text(encoding="utf-8")
                self.assertNotIn(":\\\\", text)
                self.assertNotIn("/mnt/", text)
                read_json(manifest)

    def test_wsl_argument_is_project_relative_without_a_username(self):
        graph_path = WORKSPACE_ROOT / "graphs" / "psen2.healthy.gfa"
        self.assertEqual(to_wsl_arg(str(graph_path), "/unused"), "graphs/psen2.healthy.gfa")

