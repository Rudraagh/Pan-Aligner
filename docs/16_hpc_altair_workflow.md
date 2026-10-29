# 16. Laptop → Altair Access → PBS Workflow

## Purpose

The laptop holds the authoritative Git repository. The university cluster (`hpccluster`, PBS) is a remote experiment runner reached only through **Altair Access** on a lab workstation. Claude Code cannot reach the cluster; files move by upload/download through Altair (or another institution-approved mechanism), never through a public remote by default.

Nothing in this document has been verified on the cluster yet. Queue access, memory units, module names and the Python environment are **unverified** until the steps below produce real output.

## 0. Before the first job — confirm with the HPC administrator

- [ ] The account you log in with is your own, authorized account.
- [ ] Which queue students should use (`cpuq` ran a job before; `workq` returned "Queue is not enabled").
- [ ] What `Insufficient amount of resource: Qlist` means and what the job must request.
- [ ] Units of the Altair "Amount of Memory (MB)" field (once shown as `mem=1024gb`, once as 10,240 MB for input 10).
- [ ] Where job files should live (home vs scratch), quotas, and the file-transfer size limit.
- [ ] Whether a Python ≥ 3.10 environment exists (module/CondaEnv name). Default `/bin/python3` is 3.6.8 and cannot run this code.
- [ ] Whether `gpuq` is open to you, and its GPU model, driver and CUDA version (only needed later).

## 1. On the laptop — build the bundle

Commit the code you want to run first; the bundle is built from `HEAD`, not the working tree.

```bash
bash scripts/hpc/make_staging_bundle.sh
```

Output goes to `hpc_staging/` (gitignored):

| File | Purpose |
|---|---|
| `panaligner_<commit>.tar.gz` | committed repo (without `outputs/`) + pristine PanAligner, minigraph and zlib sources |
| `panaligner_<commit>.tar.gz.sha256` | checksum checked after upload |
| `bundle_manifest.json` | commit, tool commits, size, file count (keep on the laptop) |
| `pbs_env_probe.sh` | read-only environment probe job |
| `pbs_unpack_bundle.sh` | verify-and-extract job |

Use `--with-outputs` only if a later experiment needs the committed outputs (adds ~200 MB before compression).

## 2. In Altair Access — upload

Upload to one working folder on the cluster (e.g. a `panaligner_work` folder in your home or scratch):

- `panaligner_<commit>.tar.gz`
- `panaligner_<commit>.tar.gz.sha256`
- `pbs_unpack_bundle.sh`
- `pbs_env_probe.sh`

## 3. Job 1 — environment probe (no bundle needed)

Submit `pbs_env_probe.sh` with the **Shell Script** launcher from the working folder. Use the smallest resources the form allows (1 CPU, ~1 GB, 10 min) and the queue the admin confirmed.

It writes `env_probe_<date>_<jobid>.txt`: OS/glibc, CPU, storage, compilers, zlib headers, every Python on `PATH`, `module avail`, Conda envs and whether `nvidia-smi` exists. Download it and share it with Claude; it decides `MODULES` and `PYTHON_BIN` for the next jobs.

## 4. Job 2 — unpack

Submit `pbs_unpack_bundle.sh` from the same folder. It:

1. checks the upload against the `.sha256` file;
2. refuses to overwrite an existing `panaligner_<commit>/` folder;
3. extracts the bundle and checks every file against `STAGED_FILES.sha256`.

Its log is `unpack_<date>_<jobid>.log`. Stop if it reports any checksum failure.

## 5. Job 3 — CPU smoke test

Submit `panaligner_<commit>/scripts/hpc/pbs_cpu_smoke.sh` with the working directory set to `panaligner_<commit>/`. If the probe showed Python ≥ 3.10 only after loading a module, set `MODULES` (and if needed `PYTHON_BIN`) through the launcher's environment/arguments field, or edit the defaults at the top of the uploaded copy.

Set `RUN_BUILD=1` on a later run to rebuild PanAligner/minigraph for the node (`scripts/hpc/build_tools_linux.sh`, needs a C++20 compiler such as a `gcc11` module).

Results go to `panaligner_<commit>/hpc_runs/<date>_cpu_smoke_<jobid>/`:

| File | Contents |
|---|---|
| `run_manifest.json` | commit, job id, queue, per-step exit codes, `gpu_validated: false` |
| `environment.txt`, `qstat.txt`, `modules.txt`, `pip_freeze.txt` | environment record |
| `checksums.sha256` | input data/manifests as seen on the cluster |
| `imports.log`, `unit_tests.log`, `job.log` | step output |

A missing Python ≥ 3.10 writes `FAILED` and exits with code 3; that is an environment result, not a code failure.

## 6. Retrieve and record

Download the whole `hpc_runs/<run>/` folder (and the probe/unpack logs) into the laptop's `hpc_runs/` (gitignored). Check:

- `run_manifest.json` `commit` matches the bundle's `bundle_manifest.json`;
- `checksums.sha256` matches the laptop (`sha256sum -c` from the repository root).

Report failed and skipped steps as they are. A CPU smoke pass says nothing about GPU execution.

## 7. Later — GPU validation (only after the admin confirms `gpuq`)

Probe on `gpuq` first (allocated GPU, `nvidia-smi`, CUDA version), then install the matching CuPy wheel into a per-user environment, run `tests/test_proposed_minimizer_engine.py` and confirm the CUDA parity test **runs** (not skipped) and passes. Benchmarks come after parity, with warm-up, repetitions, and separate kernel/transfer/end-to-end timing.
