"""Minimal reader for standard MATPOWER case files (numeric bus/gen/branch matrices)."""
from __future__ import annotations

import re
from pathlib import Path
import numpy as np


def _read_matrix(text: str, field: str) -> np.ndarray:
    pattern = rf"mpc\.{re.escape(field)}\s*=\s*\[(.*?)\]\s*;"
    match = re.search(pattern, text, flags=re.IGNORECASE | re.DOTALL)
    if not match:
        raise ValueError(f"Không tìm thấy mpc.{field} = [...] trong file MATPOWER.")

    body = match.group(1)
    rows: list[list[float]] = []
    for raw_line in body.splitlines():
        line = raw_line.split("%", 1)[0].strip()
        if not line:
            continue
        # MATPOWER case data uses numeric literals separated by whitespace or commas.
        line = line.replace("...", " ")
        line = line.rstrip(";").strip()
        if not line:
            continue
        tokens = [tok for tok in re.split(r"[\s,]+", line) if tok]
        try:
            rows.append([float(tok.replace("D", "E").replace("d", "e")) for tok in tokens])
        except ValueError as exc:
            raise ValueError(
                f"Không đọc được dữ liệu số trong mpc.{field}: {raw_line.strip()}\n"
                "Bộ đọc hiện hỗ trợ ma trận dữ liệu số chuẩn của MATPOWER, không hỗ trợ biểu thức MATLAB."
            ) from exc

    if not rows:
        raise ValueError(f"Ma trận mpc.{field} rỗng.")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError(f"Các hàng trong mpc.{field} không có cùng số cột.")
    return np.asarray(rows, dtype=float)


def load_matpower_case(path: str | Path) -> dict:
    """Load baseMVA, bus, gen and branch arrays from a MATPOWER .m file."""
    path = Path(path)
    text = path.read_text(encoding="utf-8", errors="replace")
    base_match = re.search(r"mpc\.baseMVA\s*=\s*([-+0-9.eEdD]+)\s*;", text)
    if not base_match:
        raise ValueError(f"Không tìm thấy mpc.baseMVA trong {path}.")
    base_mva = float(base_match.group(1).replace("D", "E").replace("d", "e"))
    bus = _read_matrix(text, "bus")
    gen = _read_matrix(text, "gen")
    branch = _read_matrix(text, "branch")
    if bus.shape[1] < 13 or gen.shape[1] < 10 or branch.shape[1] < 11:
        raise ValueError("File case không có đủ cột theo định dạng MATPOWER phiên bản 2.")
    return {"name": path.stem, "path": str(path), "baseMVA": base_mva,
            "bus": bus, "gen": gen, "branch": branch}
