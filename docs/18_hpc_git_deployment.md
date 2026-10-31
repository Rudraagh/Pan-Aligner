# 18. Git-based HPC deployment of the Model A pilot

This replaces uploading a bundle for the pilot. `docs/16_hpc_altair_workflow.md` still describes the jobs, and `outputs/model_a_pilot_plan/20261030_plan_01/PILOT_PLAN.md` holds the gates.

**Jobs 1–5 are not authorised by completing Job 0.** Each needs its own explicit approval, after reviewing the previous job's output. Job 1 downloads about 7 GB.

**Status:** tested only on synthetic data in WSL. None of this has run on the cluster.

## What comes from where

| Need | Source |
|---|---|
| Scripts, jobs, tests, this guide, `k8` binary | Tracked in the branch |
| PanAligner, minigraph, zlib, PBSIM2, minimap2 sources | Git submodules pinned to the commits in the canonical bundle. Needs access to `github.com` from where you clone |
| Python packages | `wheelhouse/py311/`, tracked. Installs with no network |
| Python 3.11.0, GCC 11.3.0 | `/apps/anaconda3/envs/gcc`. Used read-only, never activated or changed |
| GIAB data | **Not in Git.** Job 1 only, after approval |

`hpc_env.sh` puts that environment's `bin` first on `PATH`, so the default batch `python3` is never used. Binaries record its `lib` folder as a run-time path, so later jobs need no `LD_LIBRARY_PATH`.

To use the GCC 12.3.0 module instead, set `MODULES=<its name>` on every step, starting from venv setup. The module name is not recorded here; find it with `module avail gcc`. Scripts load `MODULES` after setting `PATH`, so the module takes priority for jobs. For the steps below that run directly on the login node, run `module load <its name>` first.

## Before you start

- Confirm you may use the login node for `git clone`, a small venv install and short checks.
- Confirm the queue, working storage and the `Qlist` question with the administrator. See `outputs/model_a_hpc_readiness/20261031_readiness_01/ADMIN_QUESTIONS.md`, which is not in this branch.
- Never put a password or token in a command, a remote URL or a file.

## 1. Clone and check out the exact branch

Choose the folder the administrator approves, then:

```bash
git clone --branch hpc/model-a-giab-pilot --single-branch https://github.com/Rudraagh/Pan-Aligner.git panaligner_pilot
cd panaligner_pilot
```

If it was already cloned:

```bash
git fetch origin hpc/model-a-giab-pilot
git checkout hpc/model-a-giab-pilot
git merge --ff-only origin/hpc/model-a-giab-pilot
```

## 2. Fetch the pinned tool sources

Name the five paths. The repository holds other tool links that the pilot doesn't use and that have no URL, so don't use `--recursive`.

```bash
git submodule update --init -- PanAligner minigraph zlib tools/third_party/pbsim2 tools/third_party/minimap2
```

## 3. Verify the commit and files

Set the commit you were given when this branch was approved. It can't be written in this file, because a commit can't contain its own hash.

```bash
export EXPECTED_COMMIT=<full 40-character commit hash>
test "$(git rev-parse HEAD)" = "$EXPECTED_COMMIT" && echo COMMIT_OK
sha256sum --quiet -c scripts/hpc/pilot_manifest.sha256 && echo FILES_OK
```

## 4. First preflight

```bash
bash scripts/hpc/pilot_preflight.sh --expected-commit "$EXPECTED_COMMIT" --before-venv
```

Fix every `FAIL` before continuing. Read each `WARN`. `.hpc-venv` isn't needed yet.

## 5. Create the isolated Python environment

Installs only into `.hpc-venv/` in the checkout, from `wheelhouse/py311/`, with no network. It's short, and it takes the interpreter from `hpc_env.sh`.

```bash
bash scripts/hpc/pbs_setup_venv.sh
```

It writes `hpc_runs/<date>_setup_venv_local/`. Exit codes: 3 wrong or missing Python, 4 no wheels, 5 no venv, 6 install failed, 7 imports failed.

If login-node use isn't allowed, submit it instead, with `qsub scripts/hpc/pbs_setup_venv.sh`, as permitted by the administrator. That is environment setup, not a pilot job.

## 6. Second preflight

```bash
bash scripts/hpc/pilot_preflight.sh --expected-commit "$EXPECTED_COMMIT"
```

All checks must pass. This doesn't test compute nodes.

## 7. Building tools and the synthetic run

**Job 0 does both, on a compute node.** PanAligner compiles with `-march=native`, so a binary built on the login node may crash on a compute node of another CPU type. No separate build is needed.

Only if the administrator allows compiling on the login node, you can try this first. It uses its own folders, which Job 0 doesn't read:

```bash
JOBS=2 bash scripts/hpc/build_tools_linux.sh "$PWD/.hpc-build-login"
bash scripts/ml/wg_truth/synthetic_smoke.sh "$PWD/hpc_runs/login_synthetic" .hpc-venv/bin/python .hpc-build-login/bin/PanAligner .hpc-build-login/bin/minigraph .hpc-build-login/bin/minimap2 .hpc-build-login/bin/pbsim .hpc-build-login/pbsim2_data/P6C4.model
```

Its verdict is `NO_GO`, which is expected on synthetic data.

## 8. Submit Job 0 only

Job 0 requests 4 CPUs, 8 GB and 1 hour on `cpuq`, and uses no real data and no GPU. Submit from the top of the checkout, so that `PBS_O_WORKDIR` is right:

```bash
qsub scripts/hpc/pbs_pilot_0_build_smoke.sh
```

Check on it with `qstat -u "$USER"`. If the `Qlist` error appears, stop and ask the administrator; don't change resources by guesswork.

## 9. Inspect Job 0

Results are in `pilot_data/runs/*_build_smoke_*/`, plus `.hpc-build/build_info.txt`. The build log is only in the PBS output file, `pa_pilot_smoke.o<jobid>`.

- `run_manifest.json`: `overall_exit` is 0 and `commit` is `$EXPECTED_COMMIT`.
- `job.log`: 55 unit tests pass.
- `synthetic/pilot_report.json`: verdict `NO_GO`, expected on synthetic data. A `GO` here would be suspicious.
- `build_info.txt`: all tools built, including a `paftools version under k8` line, and it names the compute node's CPU.
- `timings.tsv`: `peak_kb` holds numbers, not `unknown`.

Bring these files back for review. **Stop here.** Passing shows the tools build and the plumbing runs on the cluster. It says nothing about real data, memory on the whole graph, or GPUs.

## Limits

- Cloning needs `github.com` from the login node, including the upstream tool repositories. If it's blocked, this route doesn't work, and `docs/16` remains the fallback.
- Compute nodes need the checkout's filesystem and `/apps`. They don't need internet for Job 0.
- Git must exist on compute nodes for `commit` to be recorded in manifests. Otherwise it shows `unknown`.
- Only Python 3.11 wheels are tracked.
