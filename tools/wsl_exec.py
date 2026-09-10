from __future__ import annotations

import os
import shlex
import subprocess
import sys
from pathlib import Path


WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
WSL_ROOT = "/tmp/fyp"


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


def resolve_wsl_binary_path(wsl_binary: str, workspace_wsl_root: str) -> str:
    if wsl_binary == WSL_ROOT:
        return "."
    if wsl_binary.startswith(f"{WSL_ROOT}/"):
        suffix = wsl_binary[len(WSL_ROOT) :]
        return f".{suffix}"
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
