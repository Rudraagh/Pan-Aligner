from __future__ import annotations

import os
import shlex
import subprocess
import sys
from pathlib import Path


WORKSPACE_ROOT = Path(__file__).resolve().parents[1]
WSL_ROOT = "/tmp/fyp"


def to_wsl_arg(argument: str) -> str:
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
    return f"{WSL_ROOT}/{relative}" if relative else WSL_ROOT


def main() -> int:
    if len(sys.argv) < 2:
        print("Usage: python tools/wsl_exec.py <wsl-binary> [args...]", file=sys.stderr)
        return 1

    wsl_binary = sys.argv[1]
    converted_args = [to_wsl_arg(arg) for arg in sys.argv[2:]]
    command_parts = [wsl_binary, *converted_args]
    bash_command = "cd /tmp/fyp && exec " + " ".join(shlex.quote(part) for part in command_parts)
    completed = subprocess.run(["wsl.exe", "bash", "-lc", bash_command], check=False)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())

