from __future__ import annotations

import hashlib
from dataclasses import asdict

from .models import MinimizerRecord


def _stable_hash(kmer: str) -> int:
    """The reference key: the first 64 bits of the SHA-1 hex digest."""
    return int(hashlib.sha1(kmer.encode("ascii")).hexdigest()[:16], 16)


def _compute_cpu_minimizers(sequence: str, k: int, window_size: int, assembly_id: str) -> list[MinimizerRecord]:
    if k <= 0:
        raise ValueError("k must be positive.")
    if window_size <= 0:
        raise ValueError("window_size must be positive.")
    if len(sequence) < k:
        return []

    kmers = [(sequence[index : index + k], _stable_hash(sequence[index : index + k]), index) for index in range(0, len(sequence) - k + 1)]
    minimizers: list[MinimizerRecord] = []
    last_signature: tuple[int, int] | None = None
    for window_start in range(0, len(kmers) - window_size + 1):
        kmer, hash_value, position = min(kmers[window_start : window_start + window_size], key=lambda item: (item[1], item[2], item[0]))
        signature = (hash_value, position)
        if signature != last_signature:
            minimizers.append(MinimizerRecord(kmer, hash_value, position, window_start, assembly_id))
            last_signature = signature
    # Original boundary behaviour: if there are fewer k-mers than a window,
    # return one global minimum rather than an empty sketch.
    if not minimizers:
        kmer, hash_value, position = min(kmers, key=lambda item: (item[1], item[2], item[0]))
        minimizers.append(MinimizerRecord(kmer, hash_value, position, 0, assembly_id))
    return minimizers


def _load_usable_cupy():
    """Return CuPy only after proving a CUDA device can allocate memory."""
    try:
        import cupy as cp  # type: ignore

        if cp.cuda.runtime.getDeviceCount() < 1:
            return None, "CuPy is installed but no CUDA device is visible."
        cp.zeros(1, dtype=cp.uint8)
        return cp, "CUDA device is available."
    except Exception as exc:
        return None, f"No usable CuPy/CUDA backend: {exc}"


_CUDA_KERNELS = r'''
extern "C" __global__
void sha1_prefix64(const unsigned char* sequence, const int sequence_length,
                   const int k, unsigned long long* hashes) {
    int index = blockDim.x * blockIdx.x + threadIdx.x;
    int count = sequence_length - k + 1;
    if (index >= count) return;
    unsigned char block[64];
    #pragma unroll
    for (int i = 0; i < 64; ++i) block[i] = 0;
    for (int i = 0; i < k; ++i) block[i] = sequence[index + i];
    block[k] = 0x80;
    unsigned long long bit_length = ((unsigned long long)k) * 8ULL;
    for (int i = 0; i < 8; ++i) block[63 - i] = (unsigned char)(bit_length >> (8 * i));
    unsigned int w[80];
    for (int i = 0; i < 16; ++i) w[i] = ((unsigned int)block[4*i] << 24) | ((unsigned int)block[4*i+1] << 16) | ((unsigned int)block[4*i+2] << 8) | (unsigned int)block[4*i+3];
    for (int i = 16; i < 80; ++i) { unsigned int x = w[i-3]^w[i-8]^w[i-14]^w[i-16]; w[i] = (x << 1) | (x >> 31); }
    unsigned int a=0x67452301U, b=0xEFCDAB89U, c=0x98BADCFEU, d=0x10325476U, e=0xC3D2E1F0U;
    for (int i = 0; i < 80; ++i) {
        unsigned int f, constant;
        if (i < 20) { f = (b & c) | ((~b) & d); constant = 0x5A827999U; }
        else if (i < 40) { f = b ^ c ^ d; constant = 0x6ED9EBA1U; }
        else if (i < 60) { f = (b & c) | (b & d) | (c & d); constant = 0x8F1BBCDCU; }
        else { f = b ^ c ^ d; constant = 0xCA62C1D6U; }
        unsigned int temp = ((a << 5) | (a >> 27)) + f + e + constant + w[i];
        e=d; d=c; c=(b << 30) | (b >> 2); b=a; a=temp;
    }
    unsigned int h0 = 0x67452301U + a, h1 = 0xEFCDAB89U + b;
    hashes[index] = ((unsigned long long)h0 << 32) | (unsigned long long)h1;
}

extern "C" __global__
void window_minima(const unsigned long long* hashes, const int kmer_count,
                   const int window_size, const int window_count, int* positions) {
    int start = blockDim.x * blockIdx.x + threadIdx.x;
    if (start >= window_count) return;
    int width = window_count == 1 && kmer_count < window_size ? kmer_count : window_size;
    unsigned long long best = 0xffffffffffffffffULL;
    int best_position = 0;
    for (int offset = 0; offset < width; ++offset) {
        unsigned long long candidate = hashes[start + offset];
        // Strict comparison intentionally preserves the earliest position on ties.
        if (candidate < best) { best = candidate; best_position = start + offset; }
    }
    positions[start] = best_position;
}
'''


def _compute_gpu_minimizers(sequence: str, k: int, window_size: int, assembly_id: str, cp, timings: dict | None = None) -> list[MinimizerRecord]:
    """Run SHA-1-prefix hashing, window minima, and duplicate compaction on CUDA.

    If `timings` is a dict, CUDA-event times (ms) for each device phase are added
    to it; Python record construction afterwards is not included.
    """
    import numpy as np

    if k <= 0:
        raise ValueError("k must be positive.")
    if window_size <= 0:
        raise ValueError("window_size must be positive.")
    if len(sequence) < k:
        return []
    try:
        encoded = sequence.encode("ascii")
    except UnicodeEncodeError as exc:
        raise ValueError("GPU minimizer mode supports ASCII sequence characters only.") from exc
    if k > 55:
        raise ValueError("GPU minimizer mode currently supports k <= 55 because its SHA-1 kernel uses one message block.")
    if len(sequence) >= 2**31:
        raise ValueError("GPU minimizer mode supports sequences shorter than 2^31 bases (32-bit kernel indices).")

    marks: list[tuple[str, object]] = []

    def mark(phase: str) -> None:
        if timings is not None:
            event = cp.cuda.Event()
            event.record()
            marks.append((phase, event))

    mark("start")
    sequence_device = cp.asarray(bytearray(encoded), dtype=cp.uint8)
    mark("host_to_device")
    kmer_count = len(sequence) - k + 1
    hashes = cp.empty(kmer_count, dtype=cp.uint64)
    threads = 128
    # Scalars are passed as int32 to match the kernels' `int` parameters; CuPy
    # would otherwise pass Python ints as long long without checking.
    cp.RawKernel(_CUDA_KERNELS, "sha1_prefix64")(((kmer_count + threads - 1) // threads,), (threads,), (sequence_device, np.int32(len(sequence)), np.int32(k), hashes))
    mark("hash_kernel")
    window_count = kmer_count - window_size + 1 if kmer_count >= window_size else 1
    positions = cp.empty(window_count, dtype=cp.int32)
    cp.RawKernel(_CUDA_KERNELS, "window_minima")(((window_count + threads - 1) // threads,), (threads,), (hashes, np.int32(kmer_count), np.int32(window_size), np.int32(window_count), positions))
    mark("window_kernel")

    # Keep compaction on the device; only compacted minimizers cross to Python.
    keep = cp.empty(window_count, dtype=cp.bool_)
    keep[0] = True
    if window_count > 1:
        keep[1:] = positions[1:] != positions[:-1]
    selected_positions = positions[keep]
    selected_hashes = hashes[selected_positions]
    selected_window_starts = cp.arange(window_count, dtype=cp.int32)[keep]
    mark("compaction")
    host_positions = cp.asnumpy(selected_positions).tolist()
    host_hashes = cp.asnumpy(selected_hashes).tolist()
    host_window_starts = cp.asnumpy(selected_window_starts).tolist()
    mark("device_to_host")
    if timings is not None:
        marks[-1][1].synchronize()
        for (_, start), (phase, end) in zip(marks, marks[1:]):
            timings[f"{phase}_ms"] = cp.cuda.get_elapsed_time(start, end)
    return [
        MinimizerRecord(sequence[position : position + k], int(hash_value), int(position), int(window_start), assembly_id)
        for position, hash_value, window_start in zip(host_positions, host_hashes, host_window_starts)
    ]


def compute_minimizers(sequence: str, k: int, window_size: int, assembly_id: str, backend: str = "auto") -> dict:
    """Return reference-compatible minimizers using CPU or genuine CUDA execution."""
    normalized_backend = backend.lower()
    if normalized_backend not in {"auto", "cpu", "gpu"}:
        raise ValueError("backend must be one of: auto, cpu, gpu.")
    backend_used, backend_note = "cpu", "CPU reference implementation"
    if normalized_backend == "cpu":
        records = _compute_cpu_minimizers(sequence, k, window_size, assembly_id)
    else:
        cp, availability_note = _load_usable_cupy()
        if cp is None:
            if normalized_backend == "gpu":
                raise RuntimeError(f"GPU backend was requested but is unavailable. {availability_note} Install a CUDA-matched CuPy package and verify the NVIDIA driver.")
            records = _compute_cpu_minimizers(sequence, k, window_size, assembly_id)
            backend_note = f"Auto mode selected CPU. {availability_note}"
        else:
            try:
                records = _compute_gpu_minimizers(sequence, k, window_size, assembly_id, cp)
                backend_used, backend_note = "gpu", "CUDA SHA-1 k-mer hashing, sliding-window minima, and duplicate compaction executed on the GPU."
            except Exception as exc:
                if normalized_backend == "gpu":
                    raise RuntimeError(f"GPU minimizer execution failed: {exc}") from exc
                records = _compute_cpu_minimizers(sequence, k, window_size, assembly_id)
                backend_note = f"Auto mode selected CPU because GPU minimizer execution was unusable: {exc}"
    return {"assembly_id": assembly_id, "backend_requested": normalized_backend, "backend_used": backend_used, "backend_note": backend_note, "k": k, "window_size": window_size, "minimizer_count": len(records), "records": [asdict(item) for item in records]}
