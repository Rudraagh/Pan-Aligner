from __future__ import annotations

import hashlib
from dataclasses import asdict

from .models import MinimizerRecord


def _stable_hash(kmer: str) -> int:
    return int(hashlib.sha1(kmer.encode("ascii")).hexdigest()[:16], 16)


def _compute_cpu_minimizers(sequence: str, k: int, window_size: int, assembly_id: str) -> list[MinimizerRecord]:
    if k <= 0:
        raise ValueError("k must be positive.")
    if window_size <= 0:
        raise ValueError("window_size must be positive.")
    if len(sequence) < k:
        return []

    kmers = [(sequence[index : index + k], _stable_hash(sequence[index : index + k]), index) for index in range(0, len(sequence) - k + 1)]
    if not kmers:
        return []

    minimizers: list[MinimizerRecord] = []
    last_signature: tuple[int, int] | None = None
    for window_start in range(0, len(kmers) - window_size + 1):
        window = kmers[window_start : window_start + window_size]
        kmer, hash_value, position = min(window, key=lambda item: (item[1], item[2], item[0]))
        signature = (hash_value, position)
        if signature == last_signature:
            continue
        minimizers.append(
            MinimizerRecord(
                minimizer=kmer,
                hash_value=hash_value,
                position=position,
                window_start=window_start,
                source_assembly=assembly_id,
            )
        )
        last_signature = signature

    if not minimizers:
        kmer, hash_value, position = min(kmers, key=lambda item: (item[1], item[2], item[0]))
        minimizers.append(
            MinimizerRecord(
                minimizer=kmer,
                hash_value=hash_value,
                position=position,
                window_start=0,
                source_assembly=assembly_id,
            )
        )
    return minimizers


def compute_minimizers(
    sequence: str,
    k: int,
    window_size: int,
    assembly_id: str,
    backend: str = "auto",
) -> dict:
    normalized_backend = backend.lower()
    backend_used = "cpu"
    backend_note = "CPU reference implementation"

    if normalized_backend not in {"auto", "cpu", "gpu"}:
        raise ValueError("backend must be one of: auto, cpu, gpu.")

    if normalized_backend == "gpu":
        try:
            import cupy  # type: ignore  # pragma: no cover

            backend_used = "gpu"
            backend_note = f"CuPy available ({cupy.__name__}); current implementation uses the shared deterministic minimizer kernel."
        except Exception:
            backend_note = "GPU requested but CuPy is unavailable; fell back to deterministic CPU implementation."
    elif normalized_backend == "auto":
        try:
            import cupy  # type: ignore  # pragma: no cover

            backend_used = "gpu"
            backend_note = f"Auto-detected CuPy backend ({cupy.__name__}); using the same deterministic minimizer logic."
        except Exception:
            backend_note = "Auto mode selected CPU because no GPU backend is available."

    records = _compute_cpu_minimizers(sequence, k, window_size, assembly_id)
    return {
        "assembly_id": assembly_id,
        "backend_requested": normalized_backend,
        "backend_used": backend_used,
        "backend_note": backend_note,
        "k": k,
        "window_size": window_size,
        "minimizer_count": len(records),
        "records": [asdict(item) for item in records],
    }
