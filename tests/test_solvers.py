import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from iterative_solvers import jacobi, gauss_seidel, sor, ssor


def test_all_iterative_solvers_on_spd_matrix():
    A = sparse.csr_matrix([[4.0, -1.0, 0.0], [-1.0, 4.0, -1.0], [0.0, -1.0, 3.0]])
    b = np.array([15.0, 10.0, 10.0])
    exact = spsolve(A, b)
    methods = [
        jacobi(A, b, tol=1e-11),
        gauss_seidel(A, b, tol=1e-11),
        sor(A, b, omega=1.15, tol=1e-11),
        ssor(A, b, omega=1.15, tol=1e-11),
    ]
    for result in methods:
        assert result.converged, (result.method, result.residual_relative)
        assert np.linalg.norm(result.x - exact, ord=np.inf) < 1e-9


def test_matpower_ieee14_dc_model():
    root = Path(__file__).resolve().parents[1]
    from matpower_parser import load_matpower_case
    from dc_pf import build_dc_model, solve_direct_reference
    case = load_matpower_case(root / "data" / "case14.m")
    model = build_dc_model(case)
    x = solve_direct_reference(model)
    rel = np.linalg.norm(model.rhs - model.Bred @ x) / np.linalg.norm(model.rhs)
    assert model.Bred.shape == (13, 13)
    assert rel < 1e-10
