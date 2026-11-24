"""Build and solve the MATPOWER-style DC power-flow linear system."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


@dataclass
class DCModel:
    case_name: str
    base_mva: float
    bus_ids: np.ndarray
    bus_types: np.ndarray
    ref_local: np.ndarray
    nonref_local: np.ndarray
    theta_ref_rad: np.ndarray
    Bbus: sparse.csr_matrix
    Bred: sparse.csr_matrix
    P_spec_pu: np.ndarray
    Pbusinj_pu: np.ndarray
    rhs: np.ndarray
    Bf: sparse.csr_matrix
    Pfinj_pu: np.ndarray
    branch_from: np.ndarray
    branch_to: np.ndarray
    branch_x: np.ndarray
    branch_tap: np.ndarray
    branch_shift_deg: np.ndarray


def build_dc_model(case: dict) -> DCModel:
    """Build Bbus/Bf and reduced system, including tap ratios and phase shifts.

    MATPOWER DC convention:
        Pbus_pu = Bbus @ Va_rad + Pbusinj_pu
        Pf_pu   = Bf @ Va_rad + Pfinj_pu
    """
    base = float(case["baseMVA"])
    bus = np.asarray(case["bus"], dtype=float)
    gen = np.asarray(case["gen"], dtype=float)
    branch = np.asarray(case["branch"], dtype=float)

    if base <= 0:
        raise ValueError("baseMVA phải dương.")

    # MATPOWER bus columns: 0=bus_i, 1=type, 2=Pd, 8=Va(deg).
    active_bus_mask = bus[:, 1] != 4  # exclude isolated buses
    active_bus = bus[active_bus_mask]
    if active_bus.size == 0:
        raise ValueError("Case không có nút đang vận hành.")
    bus_ids = active_bus[:, 0].astype(int)
    id_to_local = {int(bus_id): i for i, bus_id in enumerate(bus_ids)}
    nb = len(bus_ids)

    # Keep online branches whose terminals are both active.
    f_ids_all = branch[:, 0].astype(int)
    t_ids_all = branch[:, 1].astype(int)
    status_all = branch[:, 10] > 0
    active_branch_mask = np.array([
        bool(status_all[k]) and int(f_ids_all[k]) in id_to_local and int(t_ids_all[k]) in id_to_local
        for k in range(len(branch))
    ])
    br = branch[active_branch_mask]
    if len(br) == 0:
        raise ValueError("Không có nhánh đang vận hành kết nối các nút hợp lệ.")
    f = np.array([id_to_local[int(v)] for v in br[:, 0]], dtype=int)
    t = np.array([id_to_local[int(v)] for v in br[:, 1]], dtype=int)
    x = br[:, 3].astype(float)
    if np.any(np.abs(x) < 1e-14):
        bad = np.where(np.abs(x) < 1e-14)[0].tolist()
        raise ValueError(f"Có nhánh đang vận hành có điện kháng x=0 (vị trí trong danh sách nhánh đang chạy: {bad}).")

    # MATPOWER ratio=0 means a transmission line, equivalent to tap=1.
    tap = br[:, 8].astype(float).copy()
    tap[np.isclose(tap, 0.0)] = 1.0
    shift_deg = br[:, 9].astype(float)
    shift_rad = np.deg2rad(shift_deg)
    b = 1.0 / (x * tap)  # branch status already filtered to in-service

    rows = np.repeat(np.arange(len(br)), 2)
    cols = np.column_stack([f, t]).ravel()
    vals = np.column_stack([np.ones(len(br)), -np.ones(len(br))]).ravel()
    Cft = sparse.coo_matrix((vals, (rows, cols)), shape=(len(br), nb)).tocsr()
    Bf = (sparse.diags(b) @ Cft).tocsr()
    Bbus = (Cft.T @ Bf).tocsr()
    Pfinj = -b * shift_rad
    Pbusinj = np.asarray(Cft.T @ Pfinj).ravel()

    # Specified net injections: sum online generator Pg - bus Pd, in p.u.
    p_mw = -active_bus[:, 2].astype(float).copy()
    for g in gen:
        if g[7] > 0 and int(g[0]) in id_to_local:
            p_mw[id_to_local[int(g[0])]] += float(g[1])
    P_spec = p_mw / base

    # One reference per connected component is required for a unique DC angle solution.
    bus_type = active_bus[:, 1].astype(int)
    ref = np.where(bus_type == 3)[0]
    if len(ref) == 0:
        raise ValueError("Không tìm thấy nút cân bằng (bus type = 3) trong case.")
    ncomp, labels = sparse.csgraph.connected_components(
        sparse.csr_matrix((np.ones(2 * len(br)), (np.r_[f, t], np.r_[t, f])), shape=(nb, nb)),
        directed=False,
    )
    for component in range(ncomp):
        refs_here = np.sum(labels[ref] == component)
        if refs_here != 1:
            raise ValueError(
                f"Thành phần liên thông {component} có {refs_here} nút cân bằng; "
                "cần đúng một nút cân bằng cho mỗi thành phần."
            )

    nonref = np.array([i for i in range(nb) if i not in set(ref.tolist())], dtype=int)
    theta_ref = np.deg2rad(active_bus[ref, 8].astype(float))
    # Reduced equation: B_nn * theta_n = P_n - Pbusinj_n - B_nr * theta_ref.
    Bred = Bbus[nonref, :][:, nonref].tocsr()
    Bnr = Bbus[nonref, :][:, ref]
    rhs = P_spec[nonref] - Pbusinj[nonref] - np.asarray(Bnr @ theta_ref).ravel()

    if Bred.shape[0] == 0:
        raise ValueError("Không có nút không cân bằng để giải.")
    if np.any(np.abs(Bred.diagonal()) < 1e-14):
        raise ValueError("Ma trận B' có phần tử đường chéo bằng 0; cần kiểm tra dữ liệu lưới.")

    return DCModel(
        case_name=case.get("name", "case"), base_mva=base,
        bus_ids=bus_ids, bus_types=bus_type, ref_local=ref,
        nonref_local=nonref, theta_ref_rad=theta_ref,
        Bbus=Bbus, Bred=Bred, P_spec_pu=P_spec,
        Pbusinj_pu=Pbusinj, rhs=rhs, Bf=Bf, Pfinj_pu=Pfinj,
        branch_from=br[:, 0].astype(int), branch_to=br[:, 1].astype(int),
        branch_x=x, branch_tap=tap, branch_shift_deg=shift_deg,
    )


def solve_direct_reference(model: DCModel) -> np.ndarray:
    """Sparse direct reference solution for validation only."""
    x = np.asarray(spsolve(model.Bred.tocsc(), model.rhs), dtype=float)
    if not np.all(np.isfinite(x)):
        raise RuntimeError("Nghiệm tham chiếu từ sparse direct solver không hữu hạn.")
    return x


def expand_angles(model: DCModel, theta_nonref: np.ndarray) -> np.ndarray:
    theta = np.zeros(len(model.bus_ids), dtype=float)
    theta[model.ref_local] = model.theta_ref_rad
    theta[model.nonref_local] = theta_nonref
    return theta


def calculate_outputs(model: DCModel, theta_nonref: np.ndarray):
    theta = expand_angles(model, theta_nonref)
    bus_inj_pu = np.asarray(model.Bbus @ theta).ravel() + model.Pbusinj_pu
    flow_from_mw = (np.asarray(model.Bf @ theta).ravel() + model.Pfinj_pu) * model.base_mva
    bus_rows = []
    for i, bus_id in enumerate(model.bus_ids):
        bus_rows.append({
            "bus_id": int(bus_id),
            "bus_type": int(model.bus_types[i]),
            "angle_rad": float(theta[i]),
            "angle_deg": float(np.rad2deg(theta[i])),
            "net_injection_dc_MW": float(bus_inj_pu[i] * model.base_mva),
            "specified_injection_MW_nonref": (float(model.P_spec_pu[i] * model.base_mva)
                                            if i not in set(model.ref_local.tolist()) else np.nan),
            "is_reference": bool(i in set(model.ref_local.tolist())),
        })
    branch_rows = []
    for k in range(len(flow_from_mw)):
        branch_rows.append({
            "from_bus": int(model.branch_from[k]),
            "to_bus": int(model.branch_to[k]),
            "reactance_pu": float(model.branch_x[k]),
            "tap_ratio": float(model.branch_tap[k]),
            "phase_shift_deg": float(model.branch_shift_deg[k]),
            "P_from_to_MW": float(flow_from_mw[k]),
            "P_to_from_MW": float(-flow_from_mw[k]),
        })
    return bus_rows, branch_rows
