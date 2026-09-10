from __future__ import annotations

import os
import shlex
import subprocess
import sys
from pathlib import Path


WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
WSL_TOOL_SPECS = {
    "panaligner": ("PANALIGNER_WSL_BIN", "PanAligner"),
    "minigraph": ("MINIGRAPH_WSL_BIN", "minigraph"),
}


def windows_path_to_wsl(path: Path) -> str:
    completed = subprocess.run(
        ["wsl.exe", "wslpath", "-a", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"Unable to convert Windows path to WSL path: {path}\n{completed.stderr}")
    return completed.stdout.strip()


def resolve_wsl_root() -> str:
    return windows_path_to_wsl(WORKSPACE_ROOT)


def discover_wsl_binary(tool: str) -> str:
    """Find a WSL tool without assuming a particular WSL username or home path."""
    environment_variable, executable_name = WSL_TOOL_SPECS[tool]
    shell_script = (
        f'if [ -n "${{{environment_variable}:-}}" ]; then printf "%s\\n" "${{{environment_variable}}}"; '
        f'elif command -v {shlex.quote(executable_name)} >/dev/null 2>&1; then command -v {shlex.quote(executable_name)}; '
        f'else find "$HOME" -type f -name {shlex.quote(executable_name)} -perm -111 -print -quit 2>/dev/null; fi'
    )
    completed = subprocess.run(
        ["wsl.exe", "bash", "-lc", shell_script],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    discovered = completed.stdout.strip()
    if completed.returncode != 0 or not discovered:
        raise RuntimeError(
            f"Unable to locate {executable_name} in WSL. Set {environment_variable} to its executable path.\n"
            f"{completed.stderr}"
        )
    return discovered


def resolve_wsl_binary_path(wsl_binary: str, workspace_wsl_root: str) -> str:
    if wsl_binary.startswith("@"):
        tool = wsl_binary[1:]
        if tool not in WSL_TOOL_SPECS:
            raise ValueError(f"Unknown WSL tool selector: {wsl_binary}")
        return discover_wsl_binary(tool)
    return wsl_binary


def to_wsl_arg(argument: str, workspace_wsl_root: str) -> str:
    if not argument:
        return argument
    if argument.startswith("-"):
        return argument

    path = Path(argument)
    if not path.is_absolute():
        candidate = (Path.cwd() / path).resolve()
    else:
        candidate = path.resolve()

    try:
        candidate.relative_to(WORKSPACE_ROOT)
    except ValueError:
        return argument

    relative = candidate.relative_to(WORKSPACE_ROOT).as_posix()
    return relative if relative else "."


def main() -> int:
    if len(sys.argv) < 2:
        print("Usage: python tools/wsl_exec.py <wsl-binary> [args...]", file=sys.stderr)
        return 1

    workspace_wsl_root = resolve_wsl_root()
    wsl_binary = resolve_wsl_binary_path(sys.argv[1], workspace_wsl_root)
    converted_args = [to_wsl_arg(arg, workspace_wsl_root) for arg in sys.argv[2:]]
    command_parts = [wsl_binary, *converted_args]
    bash_command = "cd " + shlex.quote(workspace_wsl_root) + " && exec " + " ".join(shlex.quote(part) for part in command_parts)
    completed = subprocess.run(["wsl.exe", "bash", "-lc", bash_command], check=False)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
