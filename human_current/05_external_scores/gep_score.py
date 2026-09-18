"""Independent fixed-H GEP scorer with a sparse-safe GPU path."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
import scipy.sparse as sp
import torch

from contracts import validate_overlap


@dataclass
class GepScoreResult:
    raw_usage: np.ndarray
    normalized_usage: np.ndarray
    reconstruction_error: np.ndarray
    explained_fraction: np.ndarray
    converged: bool
    iterations: int
    kkt_residual: np.ndarray
    gene_sd: np.ndarray | None = None
    overlap_genes: int | None = None
    identity_spectrum_mass: np.ndarray | None = None
    identity_mass: np.ndarray | None = None
    functional_mass: np.ndarray | None = None


def _torch_dtype(dtype: str | torch.dtype) -> torch.dtype:
    if isinstance(dtype, torch.dtype):
        return dtype
    choices = {
        "float32": torch.float32,
        "float64": torch.float64,
        "fp32": torch.float32,
        "fp64": torch.float64,
    }
    try:
        return choices[str(dtype).lower()]
    except KeyError as exc:
        raise ValueError(f"unsupported dtype: {dtype}") from exc


def _resolve_device(device: str) -> torch.device:
    if device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA is unavailable")
    return torch.device(device)


def positive_seeded_initialization(
    rows: int,
    programs: int,
    seed: int,
    device: str | torch.device,
    dtype: torch.dtype,
) -> torch.Tensor:
    target = torch.device(device)
    generator = torch.Generator(device=target)
    generator.manual_seed(int(seed))
    return torch.rand(
        (rows, programs), device=target, dtype=dtype, generator=generator
    ).add_(0.1)


def _solve_fixed_h(
    x: torch.Tensor,
    h: torch.Tensor,
    *,
    seed: int,
    max_iter: int,
    tol: float,
    initial: torch.Tensor | None = None,
) -> tuple[torch.Tensor, bool, int, torch.Tensor, torch.Tensor]:
    if x.ndim != 2 or h.ndim != 2 or x.shape[1] != h.shape[1]:
        raise ValueError("fixed-H matrix shapes are incompatible")
    if not torch.isfinite(x).all() or not torch.isfinite(h).all():
        raise ValueError("fixed-H input contains non-finite values")
    if torch.min(x).item() < 0 or torch.min(h).item() < 0:
        raise ValueError("fixed-H input must be nonnegative")
    numerator = x @ h.transpose(0, 1)
    gram = h @ h.transpose(0, 1)
    usage = initial
    if usage is None:
        avg = torch.sqrt(torch.mean(x) / h.shape[0])
        usage = torch.full(
            (x.shape[0], h.shape[0]),
            avg.item(),
            device=x.device,
            dtype=x.dtype,
        )
    epsilon = torch.finfo(x.dtype).eps
    initial_error = torch.linalg.vector_norm(usage @ h - x)
    previous_error = initial_error
    converged = False
    iterations = 0
    for iteration in range(int(max_iter)):
        denominator = usage @ gram + epsilon
        updated = usage * numerator / denominator
        usage = updated
        iterations = iteration + 1
        if iterations % 10 == 0:
            current_error = torch.linalg.vector_norm(usage @ h - x)
            if initial_error.item() == 0.0 or (
                (previous_error - current_error) / initial_error
            ).item() < tol:
                converged = True
                break
            previous_error = current_error
    kkt = torch.abs(usage @ gram - numerator).amax(dim=1)
    return usage, converged, iterations, kkt, gram


def _result_from_tensors(
    x: torch.Tensor,
    h: torch.Tensor,
    usage: torch.Tensor,
    converged: bool,
    iterations: int,
    kkt: torch.Tensor,
    gram: torch.Tensor,
) -> GepScoreResult:
    epsilon = torch.finfo(x.dtype).eps
    reconstruction = usage @ h
    residual = torch.linalg.vector_norm(reconstruction - x, dim=1)
    x_norm = torch.linalg.vector_norm(x, dim=1)
    explained = 1.0 - (residual / (x_norm + epsilon)).square()
    explained = torch.clamp(explained, min=-1.0, max=1.0)
    row_sum = usage.sum(dim=1, keepdim=True)
    normalized = torch.where(row_sum > epsilon, usage / row_sum, torch.zeros_like(usage))
    return GepScoreResult(
        raw_usage=usage.detach().cpu().numpy(),
        normalized_usage=normalized.detach().cpu().numpy(),
        reconstruction_error=residual.detach().cpu().numpy(),
        explained_fraction=explained.detach().cpu().numpy(),
        converged=converged,
        iterations=iterations,
        kkt_residual=kkt.detach().cpu().numpy(),
    )


def score_dense_fixed_h(
    x,
    h,
    device: str = "cuda",
    dtype: str | torch.dtype = "float32",
    seed: int = 20260822,
    max_iter: int = 1000,
    tol: float = 1e-4,
    initial_mean: float | None = None,
) -> GepScoreResult:
    """Solve ``min_W ||X-WH||²`` for an already scaled dense X matrix."""
    target = _resolve_device(device)
    torch_dtype = _torch_dtype(dtype)
    x_tensor = torch.as_tensor(np.asarray(x), dtype=torch_dtype, device=target)
    h_tensor = torch.as_tensor(np.asarray(h), dtype=torch_dtype, device=target)
    initial = None
    if initial_mean is not None:
        avg = float(np.sqrt(max(float(initial_mean), 0.0) / h_tensor.shape[0]))
        initial = torch.full(
            (x_tensor.shape[0], h_tensor.shape[0]),
            avg,
            device=target,
            dtype=torch_dtype,
        )
    usage, converged, iterations, kkt, gram = _solve_fixed_h(
        x_tensor,
        h_tensor,
        seed=seed,
        max_iter=max_iter,
        tol=tol,
        initial=initial,
    )
    return _result_from_tensors(
        x_tensor, h_tensor, usage, converged, iterations, kkt, gram
    )


def _gene_sd(counts) -> np.ndarray:
    rows = counts.shape[0]
    if rows <= 0:
        raise ValueError("query must contain at least one observation")
    if sp.issparse(counts):
        numeric = counts.astype(np.float64, copy=False)
        total = np.asarray(numeric.sum(axis=0)).ravel()
        square_total = np.asarray(numeric.multiply(numeric).sum(axis=0)).ravel()
    else:
        numeric = np.asarray(counts, dtype=np.float64)
        total = numeric.sum(axis=0)
        square_total = np.square(numeric).sum(axis=0)
    variance = square_total / rows - np.square(total / rows)
    variance = np.maximum(variance, 0.0)
    return np.sqrt(variance)


def _index_overlap(query_gene_ids: Iterable[str], reference_gene_ids: Iterable[str]):
    query = [str(x) for x in query_gene_ids]
    reference = [str(x) for x in reference_gene_ids]
    positions = {gene: i for i, gene in enumerate(query)}
    query_idx = []
    reference_idx = []
    for ref_i, gene in enumerate(reference):
        query_i = positions.get(gene)
        if query_i is not None:
            query_idx.append(query_i)
            reference_idx.append(ref_i)
    return np.asarray(query_idx, dtype=np.int64), np.asarray(reference_idx, dtype=np.int64)


def _sparse_scaled_batch(counts, query_idx: np.ndarray, inverse_sd: np.ndarray):
    batch = counts[:, query_idx].astype(np.float32, copy=False)
    return batch.multiply(inverse_sd[query_idx])


def score_counts_fixed_h(
    counts,
    query_gene_ids: Iterable[str],
    h,
    reference_gene_ids: Iterable[str],
    *,
    identity_indices: Iterable[int],
    functional_indices: Iterable[int],
    device: str = "cuda",
    dtype: str | torch.dtype = "float32",
    seed: int = 20260822,
    max_iter: int = 1000,
    tol: float = 1e-4,
    batch_size: int = 4096,
) -> GepScoreResult:
    """Score raw counts with two-pass ddof=0 scaling and batched fixed-H MU."""
    query_gene_ids = [str(x) for x in query_gene_ids]
    reference_gene_ids = [str(x) for x in reference_gene_ids]
    if counts.ndim != 2:
        raise ValueError("query counts must be a two-dimensional matrix")
    if sp.issparse(counts):
        counts = counts.tocsr()
    else:
        counts = np.asarray(counts)
    if counts.shape[1] != len(query_gene_ids):
        raise ValueError("query gene axis does not match count matrix")
    query_idx, reference_idx = _index_overlap(query_gene_ids, reference_gene_ids)
    h_array = np.asarray(h, dtype=np.float64)
    if h_array.ndim != 2 or h_array.shape[1] != len(reference_gene_ids):
        raise ValueError("invalid fixed-H shape")
    identity_indices = np.asarray(list(identity_indices), dtype=np.int64)
    functional_indices = np.asarray(list(functional_indices), dtype=np.int64)
    retained_mass = h_array[:, reference_idx].sum(axis=1) / np.maximum(
        h_array.sum(axis=1), np.finfo(np.float64).eps
    )
    identity_retained = retained_mass[identity_indices]
    validate_overlap(query_idx.size, identity_retained)
    from contracts import assert_raw_counts

    assert_raw_counts(counts)
    sd = _gene_sd(counts)
    inverse_sd = np.divide(
        1.0,
        sd,
        out=np.zeros_like(sd),
        where=sd > 0,
    )
    target = _resolve_device(device)
    torch_dtype = _torch_dtype(dtype)
    h_sub = torch.as_tensor(
        h_array[:, reference_idx], dtype=torch_dtype, device=target
    )
    gram = h_sub @ h_sub.transpose(0, 1)
    n_rows = counts.shape[0]
    n_programs = h_array.shape[0]
    raw = np.zeros((n_rows, n_programs), dtype=np.float64 if torch_dtype == torch.float64 else np.float32)
    normalized = np.zeros_like(raw)
    reconstruction_error = np.zeros(n_rows, dtype=raw.dtype)
    explained = np.zeros(n_rows, dtype=raw.dtype)
    kkt = np.zeros(n_rows, dtype=raw.dtype)
    any_converged = True
    total_iterations = 0
    inverse_full = np.zeros(counts.shape[1], dtype=np.float64)
    inverse_full[query_idx] = inverse_sd[query_idx]
    if sp.issparse(counts):
        global_sum = float(np.asarray(counts.dot(inverse_full)).sum())
    else:
        global_sum = float(np.asarray(counts, dtype=np.float64).dot(inverse_full).sum())
    global_mean = max(global_sum / (n_rows * query_idx.size), 0.0)
    avg = float(np.sqrt(global_mean / n_programs))
    for start in range(0, n_rows, int(batch_size)):
        stop = min(n_rows, start + int(batch_size))
        if sp.issparse(counts):
            scaled = _sparse_scaled_batch(counts[start:stop], query_idx, inverse_sd)
            numerator_np = scaled @ h_array[:, reference_idx].T
            x_norm_sq = np.asarray(scaled.multiply(scaled).sum(axis=1)).ravel()
        else:
            scaled = np.asarray(counts[start:stop][:, query_idx], dtype=np.float64)
            scaled = scaled * inverse_sd[query_idx]
            numerator_np = scaled @ h_array[:, reference_idx].T
            x_norm_sq = np.square(scaled).sum(axis=1)
        x_tensor = torch.as_tensor(numerator_np, dtype=torch_dtype, device=target)
        usage = torch.full(
            (stop - start, n_programs),
            avg,
            device=target,
            dtype=torch_dtype,
        )
        epsilon = torch.finfo(torch_dtype).eps
        x_norm_sq_tensor = torch.as_tensor(
            x_norm_sq, dtype=torch_dtype, device=target
        )

        def frobenius_error(current_usage):
            fit_norm_sq = (current_usage @ gram * current_usage).sum(dim=1)
            cross = (current_usage * x_tensor).sum(dim=1)
            residual_sq = torch.clamp(
                x_norm_sq_tensor - 2.0 * cross + fit_norm_sq,
                min=0.0,
            )
            return torch.sqrt(residual_sq.sum())

        initial_error = frobenius_error(usage)
        previous_error = initial_error
        converged = False
        iterations = 0
        for iteration in range(int(max_iter)):
            updated = usage * x_tensor / (usage @ gram + epsilon)
            usage = updated
            iterations = iteration + 1
            if iterations % 10 == 0:
                current_error = frobenius_error(usage)
                if initial_error.item() == 0.0 or (
                    (previous_error - current_error) / initial_error
                ).item() < tol:
                    converged = True
                    break
                previous_error = current_error
        usage_cpu = usage.detach().cpu().numpy()
        row_sum = usage_cpu.sum(axis=1, keepdims=True)
        normalized_cpu = np.divide(
            usage_cpu,
            row_sum,
            out=np.zeros_like(usage_cpu),
            where=row_sum > np.finfo(usage_cpu.dtype).eps,
        )
        fitted = usage @ h_sub
        fit_norm_sq = torch.square(fitted).sum(dim=1).detach().cpu().numpy()
        cross = (usage * x_tensor).sum(dim=1).detach().cpu().numpy()
        residual_sq = np.maximum(x_norm_sq - 2.0 * cross + fit_norm_sq, 0.0)
        x_norm = np.sqrt(np.maximum(x_norm_sq, 0.0))
        error = np.sqrt(residual_sq)
        explained_batch = 1.0 - np.divide(
            residual_sq,
            np.maximum(x_norm_sq, np.finfo(np.float64).eps),
        )
        batch_kkt = torch.abs(usage @ gram - x_tensor).amax(dim=1).detach().cpu().numpy()
        raw[start:stop] = usage_cpu
        normalized[start:stop] = normalized_cpu
        reconstruction_error[start:stop] = error
        explained[start:stop] = np.clip(explained_batch, -1.0, 1.0)
        kkt[start:stop] = batch_kkt
        any_converged = any_converged and converged
        total_iterations = max(total_iterations, iterations)
    identity_mass = normalized[:, identity_indices].sum(axis=1)
    functional_mass = normalized[:, functional_indices].sum(axis=1)
    return GepScoreResult(
        raw_usage=raw,
        normalized_usage=normalized,
        reconstruction_error=reconstruction_error,
        explained_fraction=explained,
        converged=any_converged,
        iterations=total_iterations,
        kkt_residual=kkt,
        gene_sd=sd,
        overlap_genes=int(query_idx.size),
        identity_spectrum_mass=identity_retained,
        identity_mass=identity_mass,
        functional_mass=functional_mass,
    )
