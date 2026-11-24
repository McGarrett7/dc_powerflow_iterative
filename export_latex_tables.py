"""Create LaTeX table fragments from experiment CSV outputs."""
from pathlib import Path
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=Path, default=Path(__file__).resolve().parent / "results")
    args = parser.parse_args()
    stats_path = args.results_dir / "matrix_statistics.csv"
    summary_path = args.results_dir / "solver_summary.csv"
    if not stats_path.exists() or not summary_path.exists():
        raise SystemExit("Thiếu matrix_statistics.csv hoặc solver_summary.csv. Hãy chạy run_experiments.py trước.")
    stats = pd.read_csv(stats_path)
    summary = pd.read_csv(summary_path)

    out = []
    out.append(r"% Tự động sinh bởi export_latex_tables.py; kiểm tra số liệu và caption trước khi đưa vào báo cáo.")
    out += [r"\begin{table}[htbp]", r"\centering",
            r"\caption{Đặc điểm cấu trúc ma trận DC Power Flow trên các hệ thống IEEE}",
            r"\label{tab:dc-matrix-statistics}", r"\begin{tabular}{lrrrrr}", r"\hline",
            r"Hệ thống & Số ẩn & Số nhánh & $nnz(B')$ & Mật độ (phần trăm) & CSR (KB) \\", r"\hline"]
    for _, row in stats.iterrows():
        out.append(
            f"{row['case']} & {int(row['n_unknown_angles'])} & {int(row['n_branches_in_service'])} & "
            f"{int(row['nnz_B_reduced'])} & {row['density_percent']:.2f} & "
            f"{row['csr_storage_bytes']/1024:.2f} " + r"\\"
        )
    out += [r"\hline", r"\end{tabular}", r"\end{table}", ""]

    out += [r"\begin{table}[htbp]", r"\centering",
            r"\caption{So sánh kết quả các phương pháp lặp trên bài toán DC Power Flow}",
            r"\label{tab:dc-solver-comparison}", r"\small", r"\begin{tabular}{llrrrr}", r"\hline",
            r"Hệ thống & Phương pháp & $\omega$ & Số vòng & Phần dư tương đối & Thời gian (s) \\", r"\hline"]
    for _, row in summary.iterrows():
        omega = "--" if pd.isna(row["omega"]) else f"{row['omega']:.3f}"
        out.append(
            f"{row['case']} & {row['method']} & {omega} & {int(row['iterations'])} & "
            f"{row['residual_relative']:.3e} & {row['time_seconds_median']:.6f} " + r"\\"
        )
    out += [r"\hline", r"\end{tabular}", r"\end{table}", ""]

    target = args.results_dir / "tables_chapter5.tex"
    target.write_text("\n".join(out), encoding="utf-8")
    print(f"Đã tạo {target}")


if __name__ == "__main__":
    main()
