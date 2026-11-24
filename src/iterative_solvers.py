"""Jacobi, Gauss-Seidel, SOR and SSOR solvers for CSR matrices."""
from __future__ import annotations

from dataclasses import dataclass
import time
import numpy as np
from scipy import sparse


@dataclass
class SolverResult:
    method: str
    x: np.ndarray
    converged: bool
    iterations: int
    residual_relative: float
    elapsed_seconds: float
    history: list[tuple[int, float]]
    omega: float | None = None
    message: str = ""


def _validate(A, b, x0, tol, max_iter):
    A = sparse.csr_matrix(A, dtype=float)
    b = np.asarray(b, dtype=float).ravel()
    if A.shape[0] != A.shape[1] or A.shape[0] != b.size:
        raise ValueError("Kích thước A và b không tương thích.")
    if not np.all(np.isfinite(A.data)) or not np.all(np.isfinite(b)):
        raise ValueError("A và b phải chứa các giá trị hữu hạn.")
    diag = A.diagonal()
    if np.any(np.abs(diag) < 1e-14):
        raise ValueError("Ma trận có phần tử đường chéo bằng 0.")
    if tol <= 0:
        raise ValueError("tol phải dương.")
    if max_iter < 1:
        raise ValueError("max_iter phải ít nhất bằng 1.")
    x = np.zeros_like(b) if x0 is None else np.asarray(x0, dtype=float).copy()
    if x.shape != b.shape:
        raise ValueError("x0 không tương thích với b.")
    bnorm = np.linalg.norm(b)
    scale = bnorm if bnorm > 0 else 1.0
    return A, b, x, diag, scale


def _relative_residual(A, b, x, scale):
    return float(np.linalg.norm(b - A @ x) / scale)


def jacobi(A, b, x0=None, tol=1e-8, max_iter=200000, check_every=1) -> SolverResult:
    A, b, x, diag, scale = _validate(A, b, x0, tol, max_iter)
    offdiag_x = None
    history = []
    start = time.perf_counter()
    residual = _relative_residual(A, b, x, scale)
    history.append((0, residual))
    if residual <= tol:
        return SolverResult("Jacobi", x, True, 0, residual, time.perf_counter()-start, history)

    for k in range(1, max_iter + 1):
        # A@x contains diagonal and off-diagonal contributions.
        x_new = (b - (A @ x - diag * x)) / diag
        x = x_new
        if k % check_every == 0 or k == max_iter:
            residual = _relative_residual(A, b, x, scale)
            history.append((k, residual))
            if not np.isfinite(residual):
                break
            if residual <= tol:
                return SolverResult("Jacobi", x, True, k, residual, time.perf_counter()-start, history)
    return SolverResult("Jacobi", x, False, k, residual, time.perf_counter()-start, history,
                        message="Đạt max_iter hoặc phần dư không hữu hạn.")


def _sor_sweep_inplace(data, indices, indptr, diag, b, x, omega, reverse=False):
    n = len(b)
    order = range(n - 1, -1, -1) if reverse else range(n)
    for i in order:
        s = b[i]
        for p in range(indptr[i], indptr[i + 1]):
            j = indices[p]
            if j != i:
                s -= data[p] * x[j]
        gs_value = s / diag[i]
        x[i] = (1.0 - omega) * x[i] + omega * gs_value


def gauss_seidel(A, b, x0=None, tol=1e-8, max_iter=200000, check_every=1) -> SolverResult:
    return sor(A, b, omega=1.0, x0=x0, tol=tol, max_iter=max_iter,
               check_every=check_every, method_name="Gauss-Seidel")


def sor(A, b, omega=1.2, x0=None, tol=1e-8, max_iter=200000,
        check_every=1, method_name="SOR") -> SolverResult:
    A, b, x, diag, scale = _validate(A, b, x0, tol, max_iter)
    if not (0.0 < omega < 2.0):
        raise ValueError("Với ma trận SPD, SOR/SSOR yêu cầu 0 < omega < 2.")
    data, indices, indptr = A.data, A.indices, A.indptr
    history = []
    start = time.perf_counter()
    residual = _relative_residual(A, b, x, scale)
    history.append((0, residual))
    if residual <= tol:
        return SolverResult(method_name, x, True, 0, residual, time.perf_counter()-start,
                            history, omega=omega)

    for k in range(1, max_iter + 1):
        _sor_sweep_inplace(data, indices, indptr, diag, b, x, omega, reverse=False)
        if k % check_every == 0 or k == max_iter:
            residual = _relative_residual(A, b, x, scale)
            history.append((k, residual))
            if not np.isfinite(residual):
                break
            if residual <= tol:
                return SolverResult(method_name, x, True, k, residual,
                                    time.perf_counter()-start, history, omega=omega)
    return SolverResult(method_name, x, False, k, residual, time.perf_counter()-start,
                        history, omega=omega, message="Đạt max_iter hoặc phần dư không hữu hạn.")


def ssor(A, b, omega=1.2, x0=None, tol=1e-8, max_iter=200000,
         check_every=1) -> SolverResult:
    A, b, x, diag, scale = _validate(A, b, x0, tol, max_iter)
    if not (0.0 < omega < 2.0):
        raise ValueError("Với ma trận SPD, SSOR yêu cầu 0 < omega < 2.")
    data, indices, indptr = A.data, A.indices, A.indptr
    history = []
    start = time.perf_counter()
    residual = _relative_residual(A, b, x, scale)
    history.append((0, residual))
    if residual <= tol:
        return SolverResult("SSOR", x, True, 0, residual, time.perf_counter()-start,
                            history, omega=omega)

    for k in range(1, max_iter + 1):
        # SSOR iteration = one forward SOR sweep followed by one backward SOR sweep.
        _sor_sweep_inplace(data, indices, indptr, diag, b, x, omega, reverse=False)
        _sor_sweep_inplace(data, indices, indptr, diag, b, x, omega, reverse=True)
        if k % check_every == 0 or k == max_iter:
            residual = _relative_residual(A, b, x, scale)
            history.append((k, residual))
            if not np.isfinite(residual):
                break
            if residual <= tol:
                return SolverResult("SSOR", x, True, k, residual,
                                    time.perf_counter()-start, history, omega=omega)
    return SolverResult("SSOR", x, False, k, residual, time.perf_counter()-start,
                        history, omega=omega, message="Đạt max_iter hoặc phần dư không hữu hạn.")
