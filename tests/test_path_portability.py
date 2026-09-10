import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "tools"))

from predictor import resolve_graph_path
from wsl_exec import WORKSPACE_ROOT, to_wsl_arg


class PathPortabilityTests(unittest.TestCase):
    def test_manifest_graph_path_resolves_from_project_root(self):
        self.assertEqual(resolve_graph_path("graphs/psen2.healthy.gfa"), ROOT / "graphs" / "psen2.healthy.gfa")

    def test_wsl_argument_is_project_relative_without_a_username(self):
        graph_path = WORKSPACE_ROOT / "graphs" / "psen2.healthy.gfa"
        self.assertEqual(to_wsl_arg(str(graph_path), "/unused"), "graphs/psen2.healthy.gfa")

