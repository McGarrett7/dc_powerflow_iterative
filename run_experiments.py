"""Run iterative DC power-flow experiments and export CSV/plots."""
from __future__ import annotations

import argparse
import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import sparse

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
from matpower_parser import load_matpower_case
from dc_pf import build_dc_model, solve_direct_reference, calculate_outputs
from iterative_solvers import jacobi, gauss_seidel, sor, ssor


def run_repeated(fn, repeats: int):
    results = [fn() for _ in range(repeats)]
    # Iteration counts and residuals are deterministic; report median runtime.
    chosen = results[-1]
    times = [r.elapsed_seconds for r in results]
    chosen.elapsed_seconds = float(np.median(times))
    return chosen, times


def omega_grid(start: float, stop: float, step: float):
    vals = np.arange(start, stop + step * 0.25, step)
    vals = [float(v) for v in vals if 0.0 < v < 2.0]
    if not any(abs(v - 1.0) < 1e-12 for v in vals):
        vals.append(1.0)
    return sorted(set(round(v, 8) for v in vals))


def sweep_omega(A, b, method, values, tol, max_iter, check_every):
    records = []
    solver = sor if method == "SOR" else ssor
    for w in values:
        result = solver(A, b, omega=w, tol=tol, max_iter=max_iter, check_every=check_every)
        records.append({
            "method": method,
            "omega": w,
            "converged": result.converged,
            "iterations": result.iterations,
            "residual_relative": result.residual_relative,
            "time_seconds": result.elapsed_seconds,
        })
    candidates = [r for r in records if r["converged"]]
    if not candidates:
        return 1.0, records
    # Choose the lowest iteration count; break ties by measured elapsed time.
    best = min(candidates, key=lambda r: (r["iterations"], r["time_seconds"]))
    return float(best["omega"]), records


def matrix_storage_bytes_csr(A):
    return int(A.data.nbytes + A.indices.nbytes + A.indptr.nbytes)


def make_case_plots(case_name, sweep_df, history_map, outdir):
    if not sweep_df.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        for method in ("SOR", "SSOR"):
            sub = sweep_df[sweep_df["method"] == method]
            if len(sub):
                ax.plot(sub["omega"], sub["iterations"], marker="o", label=method)
        ax.set_xlabel("Tham số nới lỏng ω")
        ax.set_ylabel("Số vòng lặp")
        ax.set_title(f"Khảo sát tham số nới lỏng – {case_name}")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / f"omega_sweep_{case_name}.png", dpi=180)
        plt.close(fig)

    if history_map:
        fig, ax = plt.subplots(figsize=(8, 5))
        for method, history in history_map.items():
            hist = np.asarray(history, dtype=float)
            if hist.ndim == 2 and len(hist):
                ax.semilogy(hist[:, 0], np.maximum(hist[:, 1], 1e-18), label=method)
        ax.set_xlabel("Số vòng lặp")
        ax.set_ylabel("Phần dư tương đối")
        ax.set_title(f"Quá trình hội tụ – {case_name}")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / f"convergence_{case_name}.png", dpi=180)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Thực nghiệm Jacobi/GS/SOR/SSOR cho DC Power Flow.")
    parser.add_argument("--cases", nargs="+", default=["case14", "case30", "case57", "case118", "case300"])
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--out-dir", type=Path, default=ROOT / "results")
    parser.add_argument("--tol", type=float, default=1e-8)
    parser.add_argument("--max-iter", type=int, default=200000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--check-every", type=int, default=1)
    parser.add_argument("--omega-min", type=float, default=0.80)
    parser.add_argument("--omega-max", type=float, default=1.95)
    parser.add_argument("--omega-step", type=float, default=0.05)
    parser.add_argument("--no-omega-sweep", action="store_true")
    args = parser.parse_args()

    if args.repeats < 1 or args.check_every < 1:
        parser.error("--repeats và --check-every phải >= 1")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    omega_values = omega_grid(args.omega_min, args.omega_max, args.omega_step)

    summary_rows, matrix_rows, all_sweep_rows = [], [], []
    for case_name in args.cases:
        case_file = args.data_dir / f"{case_name}.m"
        if not case_file.exists():
            print(f"[SKIP] Thiếu {case_file}. Chạy python download_cases.py để tải bộ dữ liệu MATPOWER.")
            continue
        print(f"\n=== {case_name} ===")
        case = load_matpower_case(case_file)
        model = build_dc_model(case)
        A, b = model.Bred, model.rhs
        exact = solve_direct_reference(model)
        n = A.shape[0]
        nnz = int(A.nnz)
        csr_bytes = matrix_storage_bytes_csr(A)
        dense_bytes = int(n * n * 8)
        matrix_rows.append({
            "case": case_name, "n_bus_in_service": len(model.bus_ids),
            "n_unknown_angles": n, "n_branches_in_service": len(model.branch_from),
            "nnz_B_reduced": nnz, "density_percent": 100.0 * nnz / (n*n),
            "csr_storage_bytes": csr_bytes, "dense_storage_bytes": dense_bytes,
            "csr_storage_MB": csr_bytes / 1e6, "dense_storage_MB": dense_bytes / 1e6,
            "reference_bus_ids": ";".join(str(model.bus_ids[i]) for i in model.ref_local),
        })

        best_omega = {"SOR": 1.0, "SSOR": 1.0}
        case_sweep_rows = []
        if not args.no_omega_sweep:
            for method in ("SOR", "SSOR"):
                w, records = sweep_omega(A, b, method, omega_values, args.tol,
                                         args.max_iter, args.check_every)
                best_omega[method] = w
                for rec in records:
                    rec["case"] = case_name
                    case_sweep_rows.append(rec)
                    all_sweep_rows.append(rec)
                print(f"[SWEEP] {method}: omega được chọn theo số vòng lặp = {w:.3f}")

        method_specs = [
            ("Jacobi", lambda: jacobi(A, b, tol=args.tol, max_iter=args.max_iter,
                                      check_every=args.check_every)),
            ("Gauss-Seidel", lambda: gauss_seidel(A, b, tol=args.tol,
                         max_iter=args.max_iter, check_every=args.check_every)),
            ("SOR", lambda: sor(A, b, omega=best_omega["SOR"], tol=args.tol,
                         max_iter=args.max_iter, check_every=args.check_every)),
            ("SSOR", lambda: ssor(A, b, omega=best_omega["SSOR"], tol=args.tol,
                          max_iter=args.max_iter, check_every=args.check_every)),
        ]
        history_map = {}
        for method, fn in method_specs:
            result, times = run_repeated(fn, args.repeats)
            angle_error_rad = result.x - exact
            angle_error_deg = np.rad2deg(angle_error_rad)
            row = {
                "case": case_name, "method": method,
                "omega": result.omega if result.omega is not None else np.nan,
                "converged": result.converged,
                "iterations": result.iterations,
                "residual_relative": result.residual_relative,
                "time_seconds_median": result.elapsed_seconds,
                "time_seconds_runs": ";".join(f"{v:.8g}" for v in times),
                "angle_error_l2_deg": float(np.linalg.norm(angle_error_deg)),
                "angle_error_inf_deg": float(np.max(np.abs(angle_error_deg))),
                "matrix_n": n, "matrix_nnz": nnz,
            }
            summary_rows.append(row)
            history_map[method] = result.history
            print(f"[{method:13s}] converged={result.converged} iter={result.iterations} "
                  f"res={result.residual_relative:.3e} time={result.elapsed_seconds:.4f}s "
                  f"angle_err_inf={row['angle_error_inf_deg']:.3e} deg")

        # Use best tuned SOR/SSOR result in exported physical outputs.
        best_solver = sor(A, b, omega=best_omega["SOR"], tol=args.tol,
                          max_iter=args.max_iter, check_every=args.check_every)
        if not best_solver.converged:
            best_solver = solve_result = gauss_seidel(A, b, tol=args.tol,
                            max_iter=args.max_iter, check_every=args.check_every)
        bus_rows, branch_rows = calculate_outputs(model, best_solver.x)
        pd.DataFrame(bus_rows).to_csv(args.out_dir / f"bus_angles_{case_name}.csv", index=False)
        pd.DataFrame(branch_rows).to_csv(args.out_dir / f"branch_flows_{case_name}.csv", index=False)
        pd.DataFrame(history_map[next(iter(history_map))] if False else
                     [{"method": name, "iteration": k, "relative_residual": rr}
                      for name, hist in history_map.items() for k, rr in hist]
                     ).to_csv(args.out_dir / f"convergence_history_{case_name}.csv", index=False)
        make_case_plots(case_name, pd.DataFrame(case_sweep_rows), history_map, args.out_dir)

    if not summary_rows:
        raise SystemExit("Không có case nào được chạy. Hãy kiểm tra thư mục data/.")
    summary_df = pd.DataFrame(summary_rows)
    matrix_df = pd.DataFrame(matrix_rows)
    summary_df.to_csv(args.out_dir / "solver_summary.csv", index=False)
    matrix_df.to_csv(args.out_dir / "matrix_statistics.csv", index=False)
    if all_sweep_rows:
        pd.DataFrame(all_sweep_rows).to_csv(args.out_dir / "omega_sweep.csv", index=False)

    # Overview charts
    for metric, ylabel, filename, title in [
        ("iterations", "Số vòng lặp", "iterations_by_case.png", "So sánh số vòng lặp"),
        ("time_seconds_median", "Thời gian (s)", "runtime_by_case.png", "So sánh thời gian thực thi"),
    ]:
        pivot = summary_df.pivot(index="case", columns="method", values=metric)
        ax = pivot.plot(kind="bar", figsize=(9, 5))
        ax.set_xlabel("Hệ thống IEEE")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.3)
        plt.tight_layout()
        plt.savefig(args.out_dir / filename, dpi=180)
        plt.close()

    metadata = {
        "python": sys.version,
        "platform": platform.platform(),
        "tolerance_relative_residual": args.tol,
        "max_iterations": args.max_iter,
        "repeats": args.repeats,
        "check_every": args.check_every,
        "omega_values": omega_values if not args.no_omega_sweep else [],
        "case_files": [str((args.data_dir / f"{c}.m").resolve()) for c in args.cases
                       if (args.data_dir / f"{c}.m").exists()],
    }
    (args.out_dir / "experiment_config.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\nĐã xuất kết quả tại: {args.out_dir.resolve()}")
    print("Các file chính: solver_summary.csv, matrix_statistics.csv, omega_sweep.csv, "
          "bus_angles_*.csv, branch_flows_*.csv và các biểu đồ PNG.")


if __name__ == "__main__":
    main()
