#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import os
import sys
import json
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Legacy/atomic-unit constants used in velocity_correlation.py
T_AU_FS = 2.4188843265857e-2          # fs in 1 a.u. of time
BOHR_ANG = 0.529177210903             # Angstrom in 1 bohr
BOHR_M = 5.29177210903e-11            # m
T_AU_S = 2.4188843265857e-17          # s
V_AFS_TO_AU = T_AU_FS / BOHR_ANG      # (Angstrom/fs) -> a.u. velocity
AU_DIFF_TO_MILLI_CM2_S = (BOHR_M**2 / T_AU_S) * 1e7  # a.u. diffusivity -> (10^-3) cm^2/s


class ProgressTracker:
    def __init__(self, label: str, total: Optional[int] = None, update_every: int = 10):
        self.label = label
        self.total = total if total is not None and total > 0 else None
        self.update_every = max(1, int(update_every))
        self.last_len = 0
        self.start = time.perf_counter()

    def update(self, current: int) -> None:
        if current <= 0:
            return
        if (current % self.update_every) != 0 and (self.total is None or current != self.total):
            return
        elapsed = time.perf_counter() - self.start
        if self.total is not None:
            frac = min(1.0, current / self.total)
            filled = int(frac * 24)
            bar = "#" * filled + "-" * (24 - filled)
            msg = f"\r[{self.label}] [{bar}] {current}/{self.total} ({100.0 * frac:5.1f}%) {elapsed:6.1f}s"
        else:
            msg = f"\r[{self.label}] frames read: {current} {elapsed:6.1f}s"
        pad = " " * max(0, self.last_len - len(msg))
        sys.stderr.write(msg + pad)
        sys.stderr.flush()
        self.last_len = len(msg)

    def close(self, current: int) -> None:
        self.update(current)
        if self.last_len:
            sys.stderr.write("\n")
            sys.stderr.flush()


def _next_pow_two(n: int) -> int:
    n2 = 1
    while n2 < n:
        n2 <<= 1
    return n2


def autocorr_fft(a: np.ndarray) -> np.ndarray:
    n = len(a)
    nfft = _next_pow_two(2 * n)
    A = np.fft.rfft(a, n=nfft)
    c = np.fft.irfft(A * np.conjugate(A), n=nfft)[:n]
    denom = np.arange(n, 0, -1, dtype=float)
    return c / denom


def vacf_from_velocities(vel: np.ndarray) -> np.ndarray:
    n_frames = vel.shape[0]
    # Remove only the instantaneous ensemble drift at each frame.
    # Subtracting each particle's time-average velocity can bias the long-time
    # VACF tail downward and underestimate diffusion on finite trajectories.
    v = vel - vel.mean(axis=1, keepdims=True)
    # Vectorized FFT over all particle components (much faster than Python loops).
    vm = v.reshape(n_frames, -1)  # (T, N*3)
    nfft = _next_pow_two(2 * n_frames)
    A = np.fft.rfft(vm, n=nfft, axis=0)
    c = np.fft.irfft(A * np.conjugate(A), n=nfft, axis=0)[:n_frames, :]
    denom = np.arange(n_frames, 0, -1, dtype=float)
    vacf = (c.sum(axis=1) / denom) / float(vm.shape[1])
    return vacf.astype(float, copy=False)


def configure_numpy_threads(n_threads: Optional[int]) -> None:
    if n_threads is None:
        return
    n = str(int(n_threads))
    for k in [
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ]:
        os.environ[k] = n
    print(f"[INFO] Requested numerical backend threads: {n}")


def integrate_trapz(y: np.ndarray, dt: float) -> float:
    integrate = getattr(np, "trapezoid", None)
    if integrate is None:
        integrate = np.trapz
    return float(integrate(y, dx=dt))


def cumulative_integrate_trapz(y: np.ndarray, dt: float) -> np.ndarray:
    out = np.zeros_like(y, dtype=float)
    if len(y) > 1:
        out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * dt, dtype=float)
    return out


def _save_vacf_csv(path: Path, times_fs: np.ndarray, vacf: np.ndarray, vacf_norm: np.ndarray) -> None:
    pd.DataFrame(
        {
            "time_fs": np.asarray(times_fs, dtype=float),
            "vacf": np.asarray(vacf, dtype=float),
            "vacf_normalized": np.asarray(vacf_norm, dtype=float),
        }
    ).to_csv(path, index=False)


def _with_backup_suffix(path: Path) -> Path:
    return path


def _format_time_tag(value_fs: Optional[float]) -> str:
    if value_fs is None:
        return ""
    value = float(value_fs)
    if value.is_integer():
        return f"_{int(value)}"
    text = f"{value:.6f}".rstrip("0").rstrip(".")
    return f"_{text.replace('.', 'p')}"


def _vacf_output_csv(out_dir: Path, species_label: str, max_time_fs: Optional[float]) -> Path:
    return _with_backup_suffix(out_dir / f"vacf_{species_label}{_format_time_tag(max_time_fs)}.csv")


def _msd_output_csv(out_dir: Path, species_label: str, max_time_fs: Optional[float]) -> Path:
    return _with_backup_suffix(out_dir / f"msd_{species_label}{_format_time_tag(max_time_fs)}.csv")


def load_velocities_from_long_csv(
    path: Path,
    species_col="species",
    id_col="id",
    t_col="iter",
    vx_col="vx",
    vy_col="vy",
    vz_col="vz",
    species_filter: Optional[str] = None,
) -> Tuple[np.ndarray, float]:
    header = pd.read_csv(path, nrows=0)
    columns = list(header.columns)
    need = {species_col, id_col, t_col, vx_col, vy_col, vz_col}
    if not need.issubset(set(columns)):
        raise ValueError(f"CSV is missing required columns. Found {sorted(columns)}, need at least {sorted(list(need))}")
    usecols = [c for c in [species_col, id_col, t_col, vx_col, vy_col, vz_col, "time_fs"] if c in columns]
    df = pd.read_csv(path, usecols=usecols)
    if species_filter is not None:
        df = df[df[species_col] == species_filter].copy()
        if df.empty:
            raise ValueError(f"No rows for species='{species_filter}' in {path}")
    steps = np.sort(df[t_col].unique())
    if len(steps) < 2:
        raise ValueError("Not enough time steps to compute dt")
    if "time_fs" in df.columns:
        times = np.sort(df["time_fs"].unique())
        dt_fs = float(np.median(np.diff(times)))
    else:
        dt_fs = 1.0
    ids = np.sort(df[id_col].unique())
    n_particles = len(ids)
    n_frames = len(steps)

    t_values = df[t_col].to_numpy(copy=False)
    id_values = df[id_col].to_numpy(copy=False)
    t_codes = np.searchsorted(steps, t_values)
    id_codes = np.searchsorted(ids, id_values)
    if (
        np.any(t_codes >= n_frames)
        or np.any(id_codes >= n_particles)
        or np.any(steps[t_codes] != t_values)
        or np.any(ids[id_codes] != id_values)
    ):
        raise ValueError("Found rows with unknown timestep or id while building velocity array")

    vel = np.zeros((n_frames, n_particles, 3), dtype=float)
    vel_values = df[[vx_col, vy_col, vz_col]].to_numpy(dtype=float, copy=False)
    vel[t_codes, id_codes, :] = vel_values
    return vel, dt_fs


def parse_lammps_type_map(spec: Optional[str]) -> Dict[int, str]:
    if spec is None or not spec.strip():
        return {}
    out: Dict[int, str] = {}
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if ":" not in chunk:
            raise ValueError(f"Invalid --lammps-type-map entry '{chunk}'. Expected format like '1:O,2:H'")
        k, v = chunk.split(":", 1)
        out[int(k.strip())] = v.strip()
    return out


def _build_water_triplets(
    records: List[Tuple[int, Optional[int], str]],
    has_mol: bool,
    water_o_species: str,
    water_h_species: str,
    context: str,
) -> List[Tuple[int, int, int]]:
    if has_mol:
        mol_map: Dict[int, List[Tuple[int, str]]] = {}
        for atom_id, mol, sp in records:
            if mol is None:
                raise ValueError("Internal parsing error on mol field")
            mol_map.setdefault(mol, []).append((atom_id, sp))
        triplets: List[Tuple[int, int, int]] = []
        for _, vals in sorted(mol_map.items()):
            o_ids = [a for a, s in vals if s == water_o_species]
            h_ids = [a for a, s in vals if s == water_h_species]
            if len(o_ids) != 1 or len(h_ids) != 2:
                raise ValueError(
                    f"Cannot build H2O molecules{context}: expected 1 O + 2 H per molecule."
                )
            h_ids.sort()
            triplets.append((o_ids[0], h_ids[0], h_ids[1]))
        return triplets

    seq = [(atom_id, sp) for atom_id, _, sp in records]
    if len(seq) % 3 != 0:
        raise ValueError(f"Cannot infer H2O molecules{context} without 'mol' column")
    triplets: List[Tuple[int, int, int]] = []
    for i in range(0, len(seq), 3):
        (id0, sp0), (id1, sp1), (id2, sp2) = seq[i:i + 3]
        if not (sp0 == water_o_species and sp1 == water_h_species and sp2 == water_h_species):
            raise ValueError(f"Expected O-H-H atom-id pattern to infer molecules{context}")
        triplets.append((id0, id1, id2))
    return triplets


def load_velocities_from_lammpstrj(
    path: Path,
    species: List[str],
    user_type_map: Optional[Dict[int, str]] = None,
    compute_water_self: bool = False,
    water_o_species: str = "O",
    water_h_species: str = "H",
    water_self_target: str = "com",
    mass_o: float = 15.999,
    mass_h: float = 1.008,
    max_frames: Optional[int] = None,
) -> Tuple[Dict[str, np.ndarray], float, Optional[np.ndarray], Optional[float]]:
    if not path.exists():
        raise FileNotFoundError(path)
    type_to_species = dict(user_type_map or {})
    if not type_to_species:
        type_to_species = {i + 1: sp for i, sp in enumerate(species)}

    by_species: Dict[str, List[np.ndarray]] = {sp: [] for sp in species}
    species_set = set(species)
    species_counts: Dict[str, int] = {}
    species_id_to_index: Dict[str, Dict[int, int]] = {}
    water_triplets: Optional[List[Tuple[int, int, int]]] = None
    water_id_to_index: Optional[Dict[int, int]] = None
    water_triplet_indices: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]] = None
    water_self_frames: List[np.ndarray] = []
    timesteps: List[int] = []
    box_lengths_ang: List[float] = []
    m_tot = mass_o + 2.0 * mass_h
    progress = ProgressTracker("LAMMPS VACF read", total=max_frames, update_every=10)

    with path.open("r", encoding="utf-8") as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue

            ts_line = f.readline()
            if not ts_line:
                break
            timestep = int(ts_line.strip())
            timesteps.append(timestep)

            n_hdr = f.readline().strip()
            if n_hdr != "ITEM: NUMBER OF ATOMS":
                raise ValueError(f"Unexpected LAMMPS dump structure in {path}: expected 'ITEM: NUMBER OF ATOMS', got '{n_hdr}'")
            n_atoms = int(f.readline().strip())

            box_hdr = f.readline().strip()
            if not box_hdr.startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"Unexpected LAMMPS dump structure in {path}: expected 'ITEM: BOX BOUNDS', got '{box_hdr}'")
            xb = f.readline().split()
            yb = f.readline().split()
            zb = f.readline().split()
            if len(xb) >= 2 and len(yb) >= 2 and len(zb) >= 2:
                lx = float(xb[1]) - float(xb[0])
                ly = float(yb[1]) - float(yb[0])
                lz = float(zb[1]) - float(zb[0])
                if lx > 0 and ly > 0 and lz > 0:
                    box_lengths_ang.append((lx * ly * lz) ** (1.0 / 3.0))

            atoms_hdr = f.readline().strip()
            if not atoms_hdr.startswith("ITEM: ATOMS"):
                raise ValueError(f"Unexpected LAMMPS dump structure in {path}: expected 'ITEM: ATOMS ...', got '{atoms_hdr}'")
            cols = atoms_hdr.split()[2:]
            idx = {c: i for i, c in enumerate(cols)}
            required = {"id", "type", "vx", "vy", "vz"}
            if not required.issubset(idx):
                raise ValueError(
                    f"LAMMPS ATOMS line in {path} must include {sorted(required)}; found columns {cols}"
                )

            has_mol = "mol" in idx
            id_idx = idx["id"]
            type_idx = idx["type"]
            vx_idx = idx["vx"]
            vy_idx = idx["vy"]
            vz_idx = idx["vz"]
            mol_idx = idx["mol"] if has_mol else None

            if not species_id_to_index:
                frame_rows: Dict[str, List[Tuple[int, float, float, float]]] = {sp: [] for sp in species}
                water_rows: List[Tuple[int, Optional[int], str, float, float, float]] = []
                for _ in range(n_atoms):
                    atom_line = f.readline()
                    if not atom_line:
                        raise ValueError(f"Unexpected EOF while reading atoms in {path}")
                    toks = atom_line.split()
                    atom_type = int(float(toks[type_idx]))
                    sp = type_to_species.get(atom_type)
                    if sp not in species_set:
                        continue
                    atom_id = int(float(toks[id_idx]))
                    vx = float(toks[vx_idx])
                    vy = float(toks[vy_idx])
                    vz = float(toks[vz_idx])
                    frame_rows[sp].append((atom_id, vx, vy, vz))
                    if compute_water_self and sp in {water_o_species, water_h_species}:
                        mol = int(float(toks[mol_idx])) if mol_idx is not None else None
                        water_rows.append((atom_id, mol, sp, vx, vy, vz))

                for sp in species:
                    rows = frame_rows[sp]
                    if not rows:
                        raise ValueError(
                            f"No atoms found for species '{sp}' in one frame of {path}. "
                            f"Check --species and --lammps-type-map."
                        )
                    rows.sort(key=lambda x: x[0])
                    species_counts[sp] = len(rows)
                    species_id_to_index[sp] = {row[0]: i for i, row in enumerate(rows)}
                    by_species[sp].append(np.array([[r[1], r[2], r[3]] for r in rows], dtype=float))

                if compute_water_self:
                    water_rows.sort(key=lambda x: x[0])
                    water_triplets = _build_water_triplets(
                        [(atom_id, mol, sp) for atom_id, mol, sp, _, _, _ in water_rows],
                        has_mol,
                        water_o_species,
                        water_h_species,
                        "",
                    )
                    water_ids = [atom_id for atom_id, _, _, _, _, _ in water_rows]
                    water_id_to_index = {atom_id: i for i, atom_id in enumerate(water_ids)}
                    o_idx = np.fromiter((water_id_to_index[oid] for oid, _, _ in water_triplets), dtype=int)
                    h1_idx = np.fromiter((water_id_to_index[h1id] for _, h1id, _ in water_triplets), dtype=int)
                    h2_idx = np.fromiter((water_id_to_index[h2id] for _, _, h2id in water_triplets), dtype=int)
                    water_triplet_indices = (o_idx, h1_idx, h2_idx)
                    water_vel = np.array([[vx, vy, vz] for _, _, _, vx, vy, vz in water_rows], dtype=float)
                    if water_self_target == "com":
                        water_target = (mass_o * water_vel[o_idx] + mass_h * water_vel[h1_idx] + mass_h * water_vel[h2_idx]) / m_tot
                    elif water_self_target == "o":
                        water_target = water_vel[o_idx]
                    else:
                        raise ValueError(f"Unsupported water_self_target='{water_self_target}'")
                    water_self_frames.append(water_target)
            else:
                frame_arrays = {
                    sp: np.empty((species_counts[sp], 3), dtype=float) for sp in species
                }
                filled_counts = {sp: 0 for sp in species}
                water_vel = (
                    np.empty((len(water_id_to_index), 3), dtype=float)
                    if compute_water_self and water_id_to_index is not None
                    else None
                )
                water_count = 0
                for _ in range(n_atoms):
                    atom_line = f.readline()
                    if not atom_line:
                        raise ValueError(f"Unexpected EOF while reading atoms in {path}")
                    toks = atom_line.split()
                    atom_type = int(float(toks[type_idx]))
                    sp = type_to_species.get(atom_type)
                    if sp not in species_set:
                        continue
                    atom_id = int(float(toks[id_idx]))
                    sp_index = species_id_to_index[sp].get(atom_id)
                    if sp_index is None:
                        raise ValueError(
                            f"Species '{sp}' has inconsistent atom IDs across frames in {path}"
                        )
                    vx = float(toks[vx_idx])
                    vy = float(toks[vy_idx])
                    vz = float(toks[vz_idx])
                    frame_arrays[sp][sp_index] = (vx, vy, vz)
                    filled_counts[sp] += 1
                    if water_vel is not None:
                        water_index = water_id_to_index.get(atom_id)
                        if water_index is not None:
                            water_vel[water_index] = (vx, vy, vz)
                            water_count += 1

                for sp in species:
                    if filled_counts[sp] != species_counts[sp]:
                        raise ValueError(
                            f"Species '{sp}' has inconsistent particle count across frames in {path} "
                            f"({species_counts[sp]} vs {filled_counts[sp]})."
                        )
                    by_species[sp].append(frame_arrays[sp])

                if water_vel is not None:
                    if water_count != water_vel.shape[0]:
                        raise ValueError("Water atom IDs are inconsistent across frames in the trajectory")
                    o_idx, h1_idx, h2_idx = water_triplet_indices
                    if water_self_target == "com":
                        water_target = (mass_o * water_vel[o_idx] + mass_h * water_vel[h1_idx] + mass_h * water_vel[h2_idx]) / m_tot
                    elif water_self_target == "o":
                        water_target = water_vel[o_idx]
                    else:
                        raise ValueError(f"Unsupported water_self_target='{water_self_target}'")
                    water_self_frames.append(water_target)

            if max_frames is not None and len(timesteps) >= int(max_frames):
                break
            progress.update(len(timesteps))
    progress.close(len(timesteps))

    if len(timesteps) < 2:
        raise ValueError(f"Need at least 2 frames in {path} to compute VACF")
    dt_fs = float(np.median(np.diff(np.asarray(timesteps, dtype=float))))

    vel_by_species: Dict[str, np.ndarray] = {}
    for sp, frames in by_species.items():
        if not frames:
            raise ValueError(f"No frames parsed for species '{sp}' in {path}")
        n_part = frames[0].shape[0]
        for arr in frames:
            if arr.shape[0] != n_part:
                raise ValueError(
                    f"Species '{sp}' has inconsistent particle count across frames in {path} "
                    f"({n_part} vs {arr.shape[0]})."
                )
        vel_by_species[sp] = np.stack(frames, axis=0)
    water_self_vel = np.stack(water_self_frames, axis=0) if compute_water_self else None
    box_len_ang = float(np.median(np.asarray(box_lengths_ang))) if box_lengths_ang else None
    return vel_by_species, dt_fs, water_self_vel, box_len_ang


def read_avvct_csv(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    header = pd.read_csv(path, nrows=0)
    columns = list(header.columns)
    vacf_col = None
    for cand in ["vacf", "AvVCT", "AvVCT_slice", "value", "acf"]:
        if cand in columns:
            vacf_col = cand
            break
    if vacf_col is None:
        vacf_col = columns[-1]
    t_col = "iter" if "iter" in columns else columns[0]
    df = pd.read_csv(path, usecols=[t_col, vacf_col])
    return df[t_col].to_numpy(copy=False), df[vacf_col].to_numpy(copy=False)


def read_saved_vacf_csv(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    header = pd.read_csv(path, nrows=0)
    columns = list(header.columns)
    if "time_fs" not in columns or "vacf" not in columns:
        raise ValueError(f"Saved VACF CSV {path} must contain at least 'time_fs' and 'vacf' columns")
    usecols = ["time_fs", "vacf"] + (["vacf_normalized"] if "vacf_normalized" in columns else [])
    df = pd.read_csv(path, usecols=usecols)
    times_fs = df["time_fs"].to_numpy(dtype=float, copy=False)
    vacf = df["vacf"].to_numpy(dtype=float, copy=False)
    if "vacf_normalized" in df.columns:
        vacf_norm = df["vacf_normalized"].to_numpy(dtype=float, copy=False)
    else:
        vacf_norm = vacf / (vacf[0] if vacf[0] != 0 else 1.0)
    return times_fs, vacf, vacf_norm


def ensure_dir(p: Path):
    p.parent.mkdir(parents=True, exist_ok=True)


def plot_vacf_normalized2_with_inset(
    output_png: Path,
    times_fs: np.ndarray,
    vacfs_norm: Dict[str, np.ndarray],
    vacfs_raw: Dict[str, np.ndarray],
    inset_xlim: Tuple[float, float] = (0, 250.0),
    inset_ylim: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    water_species_label: Optional[str] = None,
):
    is_water_self = (
        water_species_label is not None
        and water_species_label in vacfs_norm
        and water_species_label in vacfs_raw
    )
    if is_water_self:
        vacfs_norm = {water_species_label: vacfs_norm[water_species_label]}
        vacfs_raw = {water_species_label: vacfs_raw[water_species_label]}

    fig, ax_main = plt.subplots(figsize=(8.2, 5.2))
    times_ps = np.asarray(times_fs, dtype=float) / 1000.0
    dt_fs = float(np.median(np.diff(np.asarray(times_fs, dtype=float)))) if len(times_fs) >= 2 else None
    label = "$H_2O$"
    for sp, y in vacfs_norm.items():
        ax_main.plot(times_ps, y, label=label, linewidth=2)
    ax_main.axhline(0.0, linestyle="--", linewidth=1)
    ax_main.set_xlabel("Time [ps]", fontsize=16)
    ax_main.set_ylabel("VACF(t) / VACF(0)", fontsize=16)
    if title:
        ax_main.set_title(title)
    ax_main.grid(True, alpha=0.25)

    ax_in = inset_axes(
        ax_main,
        width="55%",
        height="46%",
        loc="upper right",
        borderpad=1.2,
    )

    if dt_fs is not None:
        dt_au = dt_fs / T_AU_FS
        for sp, y in vacfs_raw.items():
            d_t = cumulative_integrate_trapz(y, dt_au)
            d_t = d_t * AU_DIFF_TO_MILLI_CM2_S * 1e2
            ax_in.plot(times_ps, d_t, linewidth=2, label=label)

    ax_in.set_xlim(inset_xlim[0] / 1000.0, inset_xlim[1] / 1000.0)
    if inset_ylim is not None:
        ax_in.set_ylim(*inset_ylim)
    if is_water_self:
        ax_in.axhline(2.32, color="tab:red", linestyle=":", linewidth=2, label="D = $2.32 \\cdot 10^{-5}$ cm$^2$/s")
    ax_in.grid(True, alpha=0.25)
    ax_in.set_xlabel("Time [ps]", fontsize=16)
    ax_in.set_ylabel(r"D [$10^{-5}$ cm$^2$/s]", fontsize=16)
    ax_in.tick_params(axis="both", labelsize=8)
    ax_in.legend(fontsize=13, frameon=False, loc="upper right")
    ensure_dir(output_png)
    fig.tight_layout()
    fig.savefig(output_png, dpi=300)
    plt.close(fig)


def compute_diffusion_from_vacf(vacf: np.ndarray, dt_fs: float, dim: int = 3) -> Dict[str, float]:
    dt_au = dt_fs / T_AU_FS
    integral_au = integrate_trapz(vacf, dt_au)
    D_milli_cm2_per_s = integral_au * AU_DIFF_TO_MILLI_CM2_S
    D_cm2_per_s = D_milli_cm2_per_s * 1e-3
    D_1e5_cm2_per_s = D_cm2_per_s * 1e5
    return {
        "D_1e-3_cm^2/s": D_milli_cm2_per_s,
        "D_1e-5_cm^2/s": D_1e5_cm2_per_s,
        "cm^2_over_s": D_cm2_per_s,
    }


def _diffusion_stats_for_block_size(
    vel: np.ndarray,
    dt_fs: float,
    blk: int,
    corr_frames: int,
    dim: int = 3,
) -> Optional[Dict[str, float]]:
    n_frames = vel.shape[0]
    if blk < max(8, corr_frames):
        return None
    n_blocks_avail = n_frames // blk
    if n_blocks_avail < 2:
        return None

    dvals = []
    for b in range(n_blocks_avail):
        i0 = b * blk
        i1 = i0 + blk
        if (i1 - i0) < max(8, corr_frames):
            continue
        vacf_b = vacf_from_velocities(vel[i0:i1])
        vacf_b = vacf_b[:corr_frames]
        d_b = compute_diffusion_from_vacf(vacf_b, dt_fs, dim=dim)["D_1e-5_cm^2/s"]
        dvals.append(d_b)

    if len(dvals) < 2:
        return None

    arr = np.asarray(dvals, dtype=float)
    std = float(np.std(arr, ddof=1))
    sem = float(std / np.sqrt(arr.size))
    return {
        "blk_frames": int(blk),
        "blk_time_fs": float(blk * dt_fs),
        "n_blocks": int(arr.size),
        "std": std,
        "sem": sem,
    }


def _candidate_block_lengths(n_frames: int, min_blk: int) -> List[int]:
    max_blk = n_frames // 2
    if max_blk < min_blk:
        return []

    values = {int(min_blk), int(max_blk)}
    if max_blk > min_blk:
        geom = np.geomspace(min_blk, max_blk, num=min(18, max_blk - min_blk + 1))
        values.update(int(round(x)) for x in geom)
    candidates = sorted(v for v in values if min_blk <= v <= max_blk)
    return candidates


def _select_block_plateau(scan: List[Dict[str, float]]) -> Dict[str, float]:
    if not scan:
        return {}
    if len(scan) == 1:
        chosen = dict(scan[0])
        chosen["plateau_method"] = "single-point"
        return chosen

    sems = np.asarray([row["sem"] for row in scan], dtype=float)
    n_blocks = np.asarray([row["n_blocks"] for row in scan], dtype=int)

    for i in range(len(scan) - 2):
        tail = sems[i:i + 3]
        denom = max(float(np.mean(np.abs(tail))), 1e-16)
        rel_span = float((np.max(tail) - np.min(tail)) / denom)
        if rel_span <= 0.15 and np.all(n_blocks[i:i + 3] >= 4):
            chosen = dict(scan[i + 2])
            chosen["plateau_method"] = "stable-three-point"
            return chosen

    for row in reversed(scan):
        if int(row["n_blocks"]) >= 4:
            chosen = dict(row)
            chosen["plateau_method"] = "largest-with-four-blocks"
            return chosen

    chosen = dict(scan[-1])
    chosen["plateau_method"] = "largest-available"
    return chosen


def estimate_diffusion_error_blocks(
    vel: np.ndarray,
    dt_fs: float,
    n_blocks: int,
    max_time_fs: Optional[float] = None,
    dim: int = 3,
) -> Dict[str, float]:
    if n_blocks < 2:
        return {}
    if max_time_fs is None:
        return {}

    n_frames = vel.shape[0]
    corr_frames = int(max_time_fs / dt_fs)
    if corr_frames < 8:
        return {}
    min_blk = corr_frames
    scan = []
    for blk in _candidate_block_lengths(n_frames, min_blk):
        stats = _diffusion_stats_for_block_size(vel, dt_fs, blk, corr_frames, dim=dim)
        if stats is not None:
            scan.append(stats)

    if len(scan) < 1:
        return {}
    chosen = _select_block_plateau(scan)
    if not chosen:
        return {}

    return {
        "D_1e-5_cm^2/s_std_blocks": float(chosen["std"]),
        "D_1e-5_cm^2/s_sem_blocks": float(chosen["sem"]),
        "D_n_blocks_used": int(chosen["n_blocks"]),
        "D_block_length_frames_used": int(chosen["blk_frames"]),
        "D_block_length_fs_used": float(chosen["blk_time_fs"]),
        "D_correlation_frames_used": int(corr_frames),
        "D_correlation_time_fs_used": float(corr_frames * dt_fs),
        "D_block_plateau_method": str(chosen["plateau_method"]),
        "D_block_scan_n": int(len(scan)),
    }


KB_SI = 1.380649e-23
XI_YEHHUMMER = 2.837297


def yeh_hummer_correction_cm2_s(T_K: float, eta_mpas: float, L_ang: float) -> float:
    eta_pa_s = eta_mpas * 1e-3
    L_m = L_ang * 1e-10
    delta_m2_s = KB_SI * T_K * XI_YEHHUMMER / (6.0 * np.pi * eta_pa_s * L_m)
    return float(delta_m2_s * 1e4)


def load_water_com_positions_from_lammpstrj(
    path: Path,
    species: List[str],
    user_type_map: Optional[Dict[int, str]] = None,
    water_o_species: str = "O",
    water_h_species: str = "H",
    water_self_target: str = "com",
    mass_o: float = 15.999,
    mass_h: float = 1.008,
    max_frames: Optional[int] = None,
) -> Tuple[np.ndarray, float, float]:
    if not path.exists():
        raise FileNotFoundError(path)
    type_to_species = dict(user_type_map or {})
    if not type_to_species:
        type_to_species = {i + 1: sp for i, sp in enumerate(species)}

    timesteps: List[int] = []
    box_lengths: List[np.ndarray] = []
    wrapped_frames: List[np.ndarray] = []
    water_triplets: Optional[List[Tuple[int, int, int]]] = None
    water_id_to_index: Optional[Dict[int, int]] = None
    water_triplet_indices: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]] = None
    water_species_set = {water_o_species, water_h_species}
    m_tot = mass_o + 2.0 * mass_h
    progress = ProgressTracker("LAMMPS MSD read", total=max_frames, update_every=10)

    with path.open("r", encoding="utf-8") as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            ts_line = f.readline()
            if not ts_line:
                break
            timesteps.append(int(ts_line.strip()))

            if f.readline().strip() != "ITEM: NUMBER OF ATOMS":
                raise ValueError("Unexpected dump format: missing 'ITEM: NUMBER OF ATOMS'")
            n_atoms = int(f.readline().strip())

            box_hdr = f.readline().strip()
            if not box_hdr.startswith("ITEM: BOX BOUNDS"):
                raise ValueError("Unexpected dump format: missing 'ITEM: BOX BOUNDS'")
            xb = f.readline().split()
            yb = f.readline().split()
            zb = f.readline().split()
            lx = float(xb[1]) - float(xb[0])
            ly = float(yb[1]) - float(yb[0])
            lz = float(zb[1]) - float(zb[0])
            box_lengths.append(np.array([lx, ly, lz], dtype=float))

            atoms_hdr = f.readline().strip()
            if not atoms_hdr.startswith("ITEM: ATOMS"):
                raise ValueError("Unexpected dump format: missing 'ITEM: ATOMS ...'")
            cols = atoms_hdr.split()[2:]
            idx = {c: i for i, c in enumerate(cols)}
            required = {"id", "type", "x", "y", "z"}
            if not required.issubset(idx):
                raise ValueError(f"LAMMPS ATOMS must include {sorted(required)} for MSD")
            has_mol = "mol" in idx
            id_idx = idx["id"]
            type_idx = idx["type"]
            x_idx = idx["x"]
            y_idx = idx["y"]
            z_idx = idx["z"]
            mol_idx = idx["mol"] if has_mol else None

            if water_id_to_index is None:
                all_rows: List[Tuple[int, Optional[int], str, np.ndarray]] = []
                for _ in range(n_atoms):
                    toks = f.readline().split()
                    atom_type = int(float(toks[type_idx]))
                    sp = type_to_species.get(atom_type)
                    if sp not in water_species_set:
                        continue
                    atom_id = int(float(toks[id_idx]))
                    mol = int(float(toks[mol_idx])) if mol_idx is not None else None
                    r = np.array([float(toks[x_idx]), float(toks[y_idx]), float(toks[z_idx])], dtype=float)
                    all_rows.append((atom_id, mol, sp, r))
                all_rows.sort(key=lambda x: x[0])
                water_triplets = _build_water_triplets(
                    [(atom_id, mol, sp) for atom_id, mol, sp, _ in all_rows],
                    has_mol,
                    water_o_species,
                    water_h_species,
                    " for MSD",
                )
                water_ids = [atom_id for atom_id, _, _, _ in all_rows]
                water_id_to_index = {atom_id: i for i, atom_id in enumerate(water_ids)}
                o_idx = np.fromiter((water_id_to_index[oid] for oid, _, _ in water_triplets), dtype=int)
                h1_idx = np.fromiter((water_id_to_index[h1id] for _, h1id, _ in water_triplets), dtype=int)
                h2_idx = np.fromiter((water_id_to_index[h2id] for _, _, h2id in water_triplets), dtype=int)
                water_triplet_indices = (o_idx, h1_idx, h2_idx)
                coords = np.array([r for _, _, _, r in all_rows], dtype=float)
            else:
                coords = np.empty((len(water_id_to_index), 3), dtype=float)
                water_count = 0
                for _ in range(n_atoms):
                    toks = f.readline().split()
                    atom_type = int(float(toks[type_idx]))
                    sp = type_to_species.get(atom_type)
                    if sp not in water_species_set:
                        continue
                    atom_id = int(float(toks[id_idx]))
                    water_index = water_id_to_index.get(atom_id)
                    if water_index is None:
                        raise ValueError("Water atom IDs are inconsistent across frames in the trajectory")
                    coords[water_index] = (float(toks[x_idx]), float(toks[y_idx]), float(toks[z_idx]))
                    water_count += 1
                if water_count != coords.shape[0]:
                    raise ValueError("Water atom count is inconsistent across frames in the trajectory")

            o_idx, h1_idx, h2_idx = water_triplet_indices
            box_vec = box_lengths[-1]
            ro = coords[o_idx]
            rh1 = coords[h1_idx]
            rh2 = coords[h2_idx]
            dr1 = rh1 - ro
            dr2 = rh2 - ro
            dr1 -= box_vec * np.round(dr1 / box_vec)
            dr2 -= box_vec * np.round(dr2 / box_vec)
            rh1_img = ro + dr1
            rh2_img = ro + dr2
            if water_self_target == "com":
                target_pos = (mass_o * ro + mass_h * rh1_img + mass_h * rh2_img) / m_tot
            elif water_self_target == "o":
                target_pos = ro
            else:
                raise ValueError(f"Unsupported water_self_target='{water_self_target}'")
            wrapped_frames.append(target_pos)

            if max_frames is not None and len(timesteps) >= int(max_frames):
                break
            progress.update(len(timesteps))
    progress.close(len(timesteps))

    if len(timesteps) < 2:
        raise ValueError("Need at least 2 frames to compute MSD")
    dt_fs = float(np.median(np.diff(np.asarray(timesteps, dtype=float))))

    wrapped = np.stack(wrapped_frames, axis=0)
    box_arr = np.stack(box_lengths, axis=0)
    unwrapped = np.empty_like(wrapped)
    unwrapped[0] = wrapped[0]
    for t in range(1, wrapped.shape[0]):
        d = wrapped[t] - wrapped[t - 1]
        L = box_arr[t][None, :]
        d -= L * np.round(d / L)
        unwrapped[t] = unwrapped[t - 1] + d
    L_ang = float(np.median(np.cbrt(np.prod(box_arr, axis=1))))
    return unwrapped, dt_fs, L_ang


def msd_multi_origin(
    pos: np.ndarray,
    max_lag: Optional[int] = None,
    n_origins: int = 50,
    origin_mode: str = "uniform",
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    n_frames = pos.shape[0]
    if n_frames < 2:
        return np.zeros(1, dtype=float)
    if max_lag is None or max_lag > (n_frames - 1):
        max_lag = n_frames - 1
    if max_lag < 1:
        return np.zeros(1, dtype=float)

    max_origin = n_frames - max_lag - 1
    if max_origin < 0:
        max_origin = 0
    n_available = max_origin + 1
    if origin_mode == "all":
        origins = np.arange(n_available, dtype=int)
    elif origin_mode == "random":
        n_use = min(max(1, n_origins), n_available)
        if rng is None:
            rng = np.random.default_rng()
        origins = np.sort(rng.choice(n_available, size=n_use, replace=False).astype(int))
    elif origin_mode == "uniform":
        n_use = min(max(1, n_origins), n_available)
        origins = np.linspace(0, max_origin, n_use, dtype=int)
    else:
        raise ValueError(f"Unsupported origin_mode='{origin_mode}'")
    max_valid_stop = n_frames - origins

    msd = np.zeros(max_lag + 1, dtype=float)
    for lag in range(1, max_lag + 1):
        valid_count = np.searchsorted(max_valid_stop, lag, side="left")
        if valid_count >= origins.size:
            break
        valid = origins[:origins.size - valid_count]
        dr = pos[valid + lag] - pos[valid]
        msd[lag] = float(np.mean(np.sum(dr * dr, axis=2)))
    return msd


def diffusion_from_msd(
    msd: np.ndarray,
    dt_fs: float,
    dim: int = 3,
    fit_start_fs: Optional[float] = None,
    fit_end_fs: Optional[float] = None,
) -> Dict[str, float]:
    t = np.arange(len(msd), dtype=float) * dt_fs
    if len(t) < 3:
        return {}
    tmax = float(t[-1])
    t0 = 0.2 * tmax if fit_start_fs is None else float(fit_start_fs)
    t1 = 0.8 * tmax if fit_end_fs is None else float(fit_end_fs)
    mask = (t >= t0) & (t <= t1)
    if int(mask.sum()) < 2:
        return {}
    t_fit = t[mask]
    msd_fit = msd[mask]
    dt_centered = t_fit - t_fit.mean()
    slope = float(np.dot(dt_centered, msd_fit - msd_fit.mean()) / np.dot(dt_centered, dt_centered))
    d_cm2_s = slope / (2.0 * float(dim)) * 0.1
    return {
        "D_msd_1e-5_cm^2/s": d_cm2_s * 1e5,
        "D_msd_fit_start_fs": t0,
        "D_msd_fit_end_fs": t1,
    }


def estimate_msd_error_blocks(
    pos: np.ndarray,
    dt_fs: float,
    n_blocks: int,
    max_time_fs: Optional[float] = None,
    fit_start_fs: Optional[float] = None,
    fit_end_fs: Optional[float] = None,
    n_origins: int = 50,
    origin_mode: str = "uniform",
    rng_seed: Optional[int] = None,
) -> Dict[str, float]:
    if n_blocks < 2:
        return {}
    if max_time_fs is None:
        return {}
    n_frames = pos.shape[0]
    blk = int(max_time_fs / dt_fs)
    if blk < 16:
        return {}
    n_blocks_avail = n_frames // blk
    if n_blocks_avail < 2:
        return {}
    vals = []
    rng = np.random.default_rng(rng_seed) if origin_mode == "random" else None
    for b in range(n_blocks_avail):
        i0 = b * blk
        i1 = i0 + blk
        if (i1 - i0) < 16:
            continue
        msd_b = msd_multi_origin(pos[i0:i1], n_origins=n_origins, origin_mode=origin_mode, rng=rng)
        d_b = diffusion_from_msd(msd_b, dt_fs, fit_start_fs=fit_start_fs, fit_end_fs=fit_end_fs).get("D_msd_1e-5_cm^2/s", None)
        if d_b is not None:
            vals.append(float(d_b))
    if len(vals) < 2:
        return {}
    arr = np.asarray(vals, dtype=float)
    std = float(np.std(arr, ddof=1))
    sem = float(std / np.sqrt(arr.size))
    return {
        "D_msd_1e-5_cm^2/s_sem_blocks": sem,
        "D_msd_n_blocks_used": int(arr.size),
    }


def build_compact_diffusion_table(df: pd.DataFrame) -> pd.DataFrame:
    if "species" in df.columns:
        water_mask = df["species"].astype(str).str.startswith("H2O_")
        if water_mask.any():
            df = df[water_mask].copy()
    cols = ["species", "D_1e-5_cm^2/s"]
    if "D_1e-5_cm^2/s_sem_blocks" in df.columns:
        cols.append("D_1e-5_cm^2/s_sem_blocks")
    if "D_n_blocks_used" in df.columns:
        cols.append("D_n_blocks_used")
    if "D_msd_1e-5_cm^2/s" in df.columns:
        cols.append("D_msd_1e-5_cm^2/s")
    if "D_msd_1e-5_cm^2/s_sem_blocks" in df.columns:
        cols.append("D_msd_1e-5_cm^2/s_sem_blocks")
    if "D_msd_n_blocks_used" in df.columns:
        cols.append("D_msd_n_blocks_used")
    out = df.loc[:, [c for c in cols if c in df.columns]].copy()
    rename_map = {
        "D_1e-5_cm^2/s_sem_blocks": "err_D",
        "D_n_blocks_used": "n_blocks",
        "D_msd_1e-5_cm^2/s": "D_msd",
        "D_msd_1e-5_cm^2/s_sem_blocks": "err_D_msd",
        "D_msd_n_blocks_used": "n_blocks_msd",
    }
    out = out.rename(columns=rename_map)
    for c in out.columns:
        if c == "species":
            continue
        if c in {"n_blocks", "n_blocks_msd"}:
            out[c] = out[c].astype("Int64")
        else:
            out[c] = pd.to_numeric(out[c], errors="coerce").round(5)
    return out


def run_unified(
    input_path: str,
    species: List[str],
    mode: str = "from-vel",
    input_format: str = "auto",
    dt_fs: Optional[float] = None,
    output_dir: Optional[str] = None,
    max_time_fs: Optional[float] = None,
    msd_max_time_fs: Optional[float] = None,
    inset_xlim: Tuple[float, float] = (0.0, 250.0),
    title: Optional[str] = None,
    species_col: str = "species",
    id_col: str = "id",
    t_col: str = "iter",
    vx_col: str = "vx",
    vy_col: str = "vy",
    vz_col: str = "vz",
    lammps_type_map: Optional[Dict[int, str]] = None,
    compute_water_self: bool = False,
    water_o_species: str = "O",
    water_h_species: str = "H",
    water_self_target: str = "com",
    mass_o: float = 15.999,
    mass_h: float = 1.008,
    n_blocks_err: int = 0,
    compute_water_self_msd: bool = False,
    msd_n_origins: int = 50,
    msd_origin_mode: str = "uniform",
    msd_random_seed: Optional[int] = None,
    msd_fit_start_fs: Optional[float] = None,
    msd_fit_end_fs: Optional[float] = None,
    max_frames: Optional[int] = None,
    reuse_vacf_csv: bool = False,
    compute_vacf: bool = True,
    compute_msd: bool = False,
    diffusion_output_name: Optional[str] = None,
    diffusion_full_output_name: Optional[str] = None,
):
    if not compute_vacf and not compute_msd:
        raise ValueError("Nothing to do: enable at least one of VACF or MSD.")
    vacf_max_time_fs = max_time_fs
    water_species_label = "H2O_COM" if water_self_target == "com" else "H2O_O"
    input_path = Path(input_path)
    out_dir = Path(output_dir) if output_dir else (input_path if input_path.is_dir() else input_path.parent)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_times = None
    vacf_norm_by_species = {}
    vacf_raw_by_species = {}
    diff_rows = []

    if mode == "from-vacf":
        if not compute_vacf:
            raise ValueError("mode='from-vacf' requires VACF enabled. Remove --skip-vacf.")
        for sp in species:
            f = input_path / sp / "AvVCT_slice.csv"
            if not f.exists():
                raise FileNotFoundError(f"Expected {f} for species '{sp}'")
            iters, vacf_raw = read_avvct_csv(f)
            if dt_fs is None:
                dt_meta = input_path / "dt.json"
                if dt_meta.exists():
                    dt_fs_local = float(json.loads(dt_meta.read_text()).get("dt_fs", None))
                else:
                    raise ValueError("dt_fs not provided. Pass --dt-fs when using mode=from-vacf unless dt.json is present")
            else:
                dt_fs_local = dt_fs
            times_fs = iters * dt_fs_local if np.issubdtype(np.asarray(iters).dtype, np.integer) else iters
            vacf_norm = vacf_raw / (vacf_raw[0] if vacf_raw[0] != 0 else 1.0)
            if vacf_max_time_fs is not None:
                n_keep = int(vacf_max_time_fs / dt_fs_local)
                vacf_norm = vacf_norm[:n_keep]
                times_fs = times_fs[:n_keep]
                vacf_raw = vacf_raw[:n_keep]
            if all_times is None:
                all_times = times_fs
            vacf_norm_by_species[sp] = vacf_norm
            vacf_raw_by_species[sp] = vacf_raw
            D = compute_diffusion_from_vacf(vacf_raw, dt_fs_local, dim=3)
            diff_rows.append({"species": sp, **D})
            vacf_csv = _vacf_output_csv(out_dir, sp, vacf_max_time_fs)
            _save_vacf_csv(vacf_csv, times_fs, vacf_raw, vacf_norm)

    elif mode == "from-vel":
        is_lammps = (
            input_format == "lammps"
            or (input_format == "auto" and input_path.is_file() and input_path.suffix.lower() in {".lammpstrj", ".dump"})
        )
        compute_only_water_com = is_lammps and compute_water_self
        analysis_species = [] if (compute_only_water_com or not compute_vacf) else species

        cached_vacf: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        can_reuse_cached_species = compute_vacf and reuse_vacf_csv and n_blocks_err < 2
        if compute_vacf and reuse_vacf_csv and n_blocks_err >= 2:
            print(f"[WARN] Recomputing VACF: cached reuse disabled because --n-blocks-err={n_blocks_err} requires velocities.")
        for sp in analysis_species:
            vacf_csv = _vacf_output_csv(out_dir, sp, vacf_max_time_fs)
            if can_reuse_cached_species and vacf_csv.exists():
                times_fs, vacf, vacf_norm = read_saved_vacf_csv(vacf_csv)
                if vacf_max_time_fs is not None:
                    keep = times_fs <= vacf_max_time_fs
                    times_fs = times_fs[keep]
                    vacf = vacf[keep]
                    vacf_norm = vacf_norm[keep]
                if len(vacf) == 0:
                    raise ValueError(f"Cached VACF CSV {vacf_csv} is empty after applying max_time_fs")
                cached_vacf[sp] = (times_fs, vacf, vacf_norm)
            elif reuse_vacf_csv:
                print(f"[WARN] Recomputing VACF for species '{sp}': cached file not found at {vacf_csv}")

        water_cache = None
        can_reuse_water = compute_vacf and reuse_vacf_csv and compute_water_self and not compute_msd and n_blocks_err < 2
        water_vacf_csv = _vacf_output_csv(out_dir, water_species_label, vacf_max_time_fs)
        if compute_vacf and reuse_vacf_csv and compute_water_self and compute_msd:
            print(f"[WARN] Recomputing VACF for {water_species_label}: cached reuse disabled because --compute-water-self-msd requires trajectory data.")
        if can_reuse_water and water_vacf_csv.exists():
            times_fs, vacf, vacf_norm = read_saved_vacf_csv(water_vacf_csv)
            if vacf_max_time_fs is not None:
                keep = times_fs <= vacf_max_time_fs
                times_fs = times_fs[keep]
                vacf = vacf[keep]
                vacf_norm = vacf_norm[keep]
            if len(vacf) == 0:
                raise ValueError(f"Cached VACF CSV {water_vacf_csv} is empty after applying max_time_fs")
            water_cache = (times_fs, vacf, vacf_norm)
        elif compute_vacf and reuse_vacf_csv and compute_water_self and not compute_msd:
            print(f"[WARN] Recomputing VACF for {water_species_label}: cached file not found at {water_vacf_csv}")

        need_species_compute = any(sp not in cached_vacf for sp in analysis_species)
        need_water_compute = compute_water_self and compute_vacf and water_cache is None
        lammps_cache: Dict[str, np.ndarray] = {}
        water_self_vel = None
        dt_lammps = None
        if is_lammps and (need_species_compute or need_water_compute):
            lammps_cache, dt_lammps, water_self_vel, _ = load_velocities_from_lammpstrj(
                input_path,
                species=species,
                user_type_map=lammps_type_map,
                compute_water_self=compute_water_self,
                water_o_species=water_o_species,
                water_h_species=water_h_species,
                water_self_target=water_self_target,
                mass_o=mass_o,
                mass_h=mass_h,
                max_frames=max_frames,
            )
            if dt_fs is None:
                raise ValueError(
                    "Per input LAMMPS devi passare --dt-fs dal wrapper/script bash. "
                    f"Il dump contiene solo delta TIMESTEP={dt_lammps}, non il tempo fisico in fs."
                )

        onefile = input_path / "solute.csv" if input_path.is_dir() else None
        per_species = {sp: (input_path / f"solute_{sp}.csv") for sp in analysis_species} if input_path.is_dir() else {}
        have_one = bool(onefile and onefile.exists())
        have_per = bool(per_species) and all(p.exists() for p in per_species.values())
        if not is_lammps and analysis_species and not (have_one or have_per):
            raise FileNotFoundError(f"Could not find {onefile} nor all of {list(per_species.values())}")

        times_fs_common = None
        for sp in analysis_species:
            vel = None
            if sp in cached_vacf:
                times_fs, vacf, vacf_norm = cached_vacf[sp]
                dt_local = float(times_fs[1] - times_fs[0]) if len(times_fs) > 1 else (dt_fs if dt_fs is not None else 1.0)
            else:
                if reuse_vacf_csv:
                    print(f"[WARN] Recomputing VACF for species '{sp}': using input velocities instead of cached CSV.")
                if is_lammps:
                    vel = lammps_cache[sp] * V_AFS_TO_AU
                    dt_local = dt_fs if dt_fs is not None else dt_lammps
                elif have_one:
                    vel, dt_infer = load_velocities_from_long_csv(onefile, species_col=species_col, id_col=id_col, t_col=t_col, vx_col=vx_col, vy_col=vy_col, vz_col=vz_col, species_filter=sp)
                    dt_local = dt_fs if dt_fs is not None else dt_infer
                else:
                    vel, dt_local = load_velocities_from_long_csv(per_species[sp], species_col=species_col, id_col=id_col, t_col=t_col, vx_col=vx_col, vy_col=vy_col, vz_col=vz_col, species_filter=None)
                    if dt_fs is not None:
                        dt_local = dt_fs
                vacf = vacf_from_velocities(vel)
                times_fs = np.arange(len(vacf)) * dt_local
                if vacf_max_time_fs is not None:
                    n_keep = int(vacf_max_time_fs / dt_local)
                    vacf = vacf[:n_keep]
                    times_fs = times_fs[:n_keep]
                vacf_norm = vacf / (vacf[0] if vacf[0] != 0 else 1.0)
                vacf_csv = _vacf_output_csv(out_dir, sp, vacf_max_time_fs)
                _save_vacf_csv(vacf_csv, times_fs, vacf, vacf_norm)
            if times_fs_common is None:
                times_fs_common = times_fs
            vacf_norm_by_species[sp] = vacf_norm
            vacf_raw_by_species[sp] = vacf
            D = compute_diffusion_from_vacf(vacf, dt_local, dim=3)
            if vel is not None:
                D.update(estimate_diffusion_error_blocks(vel, dt_local, n_blocks_err, max_time_fs=vacf_max_time_fs, dim=3))
            diff_rows.append({"species": sp, **D})

        if is_lammps and compute_water_self:
            D_w = {}
            if compute_vacf:
                if water_cache is not None:
                    times_fs_w, vacf_w, vacf_norm_w = water_cache
                    dt_local = float(times_fs_w[1] - times_fs_w[0]) if len(times_fs_w) > 1 else (dt_fs if dt_fs is not None else 1.0)
                else:
                    if reuse_vacf_csv:
                        print(f"[WARN] Recomputing VACF for {water_species_label}: using trajectory data instead of cached CSV.")
                    if water_self_vel is None:
                        raise RuntimeError("Internal error: water self velocities were not computed")
                    dt_local = dt_fs if dt_fs is not None else dt_lammps
                    vacf_w = vacf_from_velocities(water_self_vel * V_AFS_TO_AU)
                    times_fs_w = np.arange(len(vacf_w)) * dt_local
                    if vacf_max_time_fs is not None:
                        n_keep = int(vacf_max_time_fs / dt_local)
                        vacf_w = vacf_w[:n_keep]
                        times_fs_w = times_fs_w[:n_keep]
                    vacf_norm_w = vacf_w / (vacf_w[0] if vacf_w[0] != 0 else 1.0)
                    _save_vacf_csv(_vacf_output_csv(out_dir, water_species_label, vacf_max_time_fs), times_fs_w, vacf_w, vacf_norm_w)
                vacf_norm_by_species[water_species_label] = vacf_norm_w
                vacf_raw_by_species[water_species_label] = vacf_w
                D_w.update(compute_diffusion_from_vacf(vacf_w, dt_local, dim=3))
                if water_cache is None:
                    D_w.update(estimate_diffusion_error_blocks(water_self_vel * V_AFS_TO_AU, dt_local, n_blocks_err, max_time_fs=vacf_max_time_fs, dim=3))
                if times_fs_common is None:
                    times_fs_common = times_fs_w
            if compute_msd:
                water_pos, dt_msd, box_l_msd = load_water_com_positions_from_lammpstrj(
                    input_path,
                    species=species,
                    user_type_map=lammps_type_map,
                    water_o_species=water_o_species,
                    water_h_species=water_h_species,
                    water_self_target=water_self_target,
                    mass_o=mass_o,
                    mass_h=mass_h,
                    max_frames=max_frames,
                )
                if dt_fs is not None and abs(float(dt_fs) - float(dt_msd)) > 1e-12:
                    print(
                        f"[INFO] Delta TIMESTEP nel dump per MSD = {dt_msd}; "
                        f"uso --dt-fs={dt_fs} come tempo fisico tra frame."
                    )
                dt_msd_use = dt_fs if dt_fs is not None else dt_msd
                max_lag = None
                if msd_max_time_fs is not None:
                    max_lag = int(msd_max_time_fs / dt_msd_use)
                rng = np.random.default_rng(msd_random_seed) if msd_origin_mode == "random" else None
                msd = msd_multi_origin(
                    water_pos,
                    max_lag=max_lag,
                    n_origins=msd_n_origins,
                    origin_mode=msd_origin_mode,
                    rng=rng,
                )
                t_msd = np.arange(len(msd), dtype=float) * dt_msd_use
                D_msd = diffusion_from_msd(
                    msd,
                    dt_msd_use,
                    dim=3,
                    fit_start_fs=msd_fit_start_fs,
                    fit_end_fs=msd_fit_end_fs,
                )
                D_msd.update(
                    estimate_msd_error_blocks(
                        water_pos,
                        dt_msd_use,
                        n_blocks=n_blocks_err,
                        max_time_fs=msd_max_time_fs,
                        fit_start_fs=msd_fit_start_fs,
                        fit_end_fs=msd_fit_end_fs,
                        n_origins=msd_n_origins,
                        origin_mode=msd_origin_mode,
                        rng_seed=msd_random_seed,
                    )
                )
                D_w.update(D_msd)
                msd_df = pd.DataFrame({"time_fs": t_msd, "msd_A2": msd})
                msd_df.to_csv(_msd_output_csv(out_dir, water_species_label, msd_max_time_fs), index=False)
            diff_rows.append({"species": water_species_label, **D_w})
        all_times = times_fs_common
    else:
        raise ValueError("mode must be 'from-vel' or 'from-vacf'")

    diff_df = pd.DataFrame(diff_rows)
    if n_blocks_err >= 2 and compute_vacf and vacf_max_time_fs is None:
        print("[WARN] Block error estimate VACF richiesta ma --max-time-fs non impostato: stima errori VACF non eseguita.")
    if n_blocks_err >= 2 and compute_msd and msd_max_time_fs is None:
        print("[WARN] Block error estimate MSD richiesta ma --msd-max-time-fs non impostato: stima errori MSD non eseguita.")
    if n_blocks_err >= 2:
        for c in ["D_1e-5_cm^2/s_std_blocks", "D_1e-5_cm^2/s_sem_blocks", "D_n_blocks_used"]:
            if c not in diff_df.columns:
                diff_df[c] = np.nan
        if float(diff_df["D_n_blocks_used"].fillna(0).max()) < 2:
            print("[WARN] Block error estimate richiesta, ma traiettoria troppo corta per usare almeno 2 blocchi validi.")
        for c in ["D_msd_1e-5_cm^2/s_sem_blocks", "D_msd_n_blocks_used"]:
            if c not in diff_df.columns:
                diff_df[c] = np.nan
    diff_compact = build_compact_diffusion_table(diff_df)
    diff_base = (
        diffusion_output_name
        if diffusion_output_name is not None
        else ("diffusion_coefficients.csv" if not compute_water_self else f"diffusion_coefficients_{water_species_label}.csv")
    )
    diff_csv = _with_backup_suffix(out_dir / diff_base)
    diff_compact.to_csv(diff_csv, index=False)
    diff_full_base = (
        diffusion_full_output_name
        if diffusion_full_output_name is not None
        else ("diffusion_coefficients_full.csv" if not compute_water_self else f"diffusion_coefficients_full_{water_species_label}.csv")
    )
    diff_csv_full = _with_backup_suffix(out_dir / diff_full_base)
    if diff_csv_full.exists():
        diff_csv_full.unlink()

    if compute_vacf and all_times is not None and vacf_norm_by_species:
        plot_base = "vacf_with_inset.pdf" if not compute_water_self else f"vacf_with_inset_{water_species_label}.pdf"
        plot_path = _with_backup_suffix(out_dir / plot_base)
        plot_vacf_normalized2_with_inset(
            plot_path,
            all_times,
            vacf_norm_by_species,
            vacf_raw_by_species,
            inset_xlim=inset_xlim,
            title=title,
            water_species_label=water_species_label if compute_water_self else None,
        )
        print(f"Saved VACF to {out_dir}/vacf_<species>[_max-time-fs].csv")
        print(f"Saved plot to {plot_path}")
    else:
        plot_path = None
    print(f"Saved diffusion coefficients (compact) to {diff_csv}")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Unified VACF + diffusion + plot tool")
    p.add_argument("input_path", help="Folder containing CSV data, or LAMMPS trajectory file (.lammpstrj/.dump)")
    p.add_argument("--mode", choices=["from-vel", "from-vacf"], default="from-vel", help="Data source mode")
    p.add_argument("--input-format", choices=["auto", "csv", "lammps"], default="auto", help="Input format for --mode from-vel")
    p.add_argument("--dt-fs", type=float, default=None, help="Time step in fs (needed if not inferable)")
    p.add_argument("--output-dir", default=None, help="Where to write outputs (default: input_path)")
    p.add_argument("--max-time-fs", type=float, default=None, help="Truncate VACF/plot at this time (fs)")
    p.add_argument("--msd-max-time-fs", type=float, default=None, help="Maximum time window for MSD calculation (fs)")
    p.add_argument("--inset-xlim", type=float, nargs=2, default=(0.0, 250.0), help="Inset xlim in fs, e.g. --inset-xlim 0 250")
    p.add_argument("--title", default=None, help="Optional plot title")

    p.add_argument("--n_species", type=int, default=None, help="Number of species (e.g., 2). If provided without --species, you'll be prompted.")
    p.add_argument("--species", nargs="*", default=None, help="List of species names, e.g., --species Na Cl")

    p.add_argument("--species-col", default="species")
    p.add_argument("--id-col", default="id")
    p.add_argument("--t-col", default="iter")
    p.add_argument("--vx-col", default="vx")
    p.add_argument("--vy-col", default="vy")
    p.add_argument("--vz-col", default="vz")
    p.add_argument("--lammps-type-map", default=None, help="Map LAMMPS atom types to species, e.g. '1:O,2:H'")
    p.add_argument("--compute-water-self", action="store_true", help="Also compute H2O self-diffusion from molecular COM (LAMMPS mode)")
    p.add_argument("--water-o-species", default="O", help="Species label to use as oxygen for H2O COM")
    p.add_argument("--water-h-species", default="H", help="Species label to use as hydrogen for H2O COM")
    p.add_argument("--water-self-target", choices=["com", "o"], default="com", help="Target used for water self-diffusion: molecular COM or oxygen atoms")
    p.add_argument("--mass-o", type=float, default=15.999, help="O mass used for H2O COM velocity")
    p.add_argument("--mass-h", type=float, default=1.008, help="H mass used for H2O COM velocity")
    p.add_argument("--n-blocks-err", type=int, default=0, help="Enable block error estimate when >=2; block length is max-time-fs")
    p.add_argument("--compute-water-self-msd", action="store_true", help="Also compute H2O self-diffusion from MSD/Einstein (LAMMPS)")
    p.add_argument("--compute-msd", action="store_true", help="Compute MSD/Einstein diffusion estimate")
    p.add_argument("--skip-msd", action="store_true", help="Disable MSD/Einstein diffusion estimate")
    p.add_argument("--skip-vacf", action="store_true", help="Disable VACF calculation/output")
    p.add_argument("--msd-n-origins", type=int, default=50, help="Number of time origins for MSD averaging")
    p.add_argument("--msd-origin-mode", choices=["uniform", "all", "random"], default="uniform", help="How to choose MSD time origins")
    p.add_argument("--msd-random-seed", type=int, default=None, help="Random seed used when --msd-origin-mode=random")
    p.add_argument("--msd-fit-start-fs", type=float, default=None, help="Linear-fit start time for MSD diffusion estimate (fs)")
    p.add_argument("--msd-fit-end-fs", type=float, default=None, help="Linear-fit end time for MSD diffusion estimate (fs)")
    p.add_argument("--max-frames", type=int, default=None, help="Read at most this many trajectory frames (speed-up for large dumps)")
    p.add_argument("--omp-threads", type=int, default=None, help="Set OMP/BLAS threads for faster FFT/linear algebra")
    p.add_argument("--reuse-vacf-csv", action="store_true", help="Reuse existing vacf_<species>.csv outputs instead of recomputing VACF when possible")

    args = p.parse_args(argv)

    if args.species is None or len(args.species) == 0:
        if args.n_species is None:
            try:
                n = int(input("Numero di specie: ").strip())
            except Exception:
                n = 2
        else:
            n = args.n_species
        species = []
        for i in range(n):
            name = input(f"Nome specie {i + 1}: ").strip()
            while not name:
                name = input(f"Nome specie {i + 1} (non vuoto): ").strip()
            species.append(name)
        args.species = species
    else:
        if args.n_species is not None and args.n_species != len(args.species):
            print(f"[WARN] --n_species={args.n_species} ma hai passato {len(args.species)} nomi in --species. Adeguo a {len(args.species)}.")
            args.n_species = len(args.species)

    return args


def main(argv=None):
    t0 = time.perf_counter()
    if argv is None and len(sys.argv) == 1:
        help_message = """
Usage: python msd_vacf.py <input_path> [options]
------------------------------------------------
Main modes:
  --mode from-vel     Compute VACF and diffusion coefficients from velocities (CSV or LAMMPS dump)
  --mode from-vacf    Use precomputed VACF files (<species>/AvVCT_slice.csv)

Key parameters:
  --input-format <f>  auto | csv | lammps (default: auto)
  --dt-fs <value>     Time step in femtoseconds (required for from-vacf)
  --n-species <N>     Number of species (e.g., 2)
  --species <names>   List of species names separated by spaces (e.g., --species Na Cl)
  --output-dir <path> Output folder (default: input_path)
  --max-time-fs <val> Truncate VACF and plots at this maximum time
  --inset-xlim a b    x-axis range for inset zoom (default: 0 250)
  --title <string>    Optional title for the plot

CSV column options:
  --species-col <name>  Column name for species (default: species)
  --id-col <name>       Column name for particle ID (default: id)
  --t-col <name>        Column name for time/iteration (default: iter)
  --vx-col, --vy-col, --vz-col  Column names for velocity components
  --lammps-type-map   Type->species map for LAMMPS dump, e.g. "1:O,2:H"
  --compute-water-self Compute also H2O COM VACF/diffusion (LAMMPS)

Units:
  - LAMMPS input velocities (A/fs) are converted to a.u. velocity internally.
  - Diffusion output is legacy-compatible and includes D_1e-5_cm^2/s (paper units).
  - Optional Yeh-Hummer corrected output: D_1e-5_cm^2/s_corrected.
  - Optional statistical uncertainty from block averaging: --n-blocks-err.
  - Optional Einstein estimate for H2O self-diffusion: --compute-water-self-msd.
  - MSD origins can be sampled uniformly, randomly, or using all available origins.

Typical usage:
  python msd_vacf.py data/ --mode from-vel --dt-fs 0.25 --n-species 2 --species Na Cl
  python msd_vacf.py traj_prod.lammpstrj --mode from-vel --input-format lammps --dt-fs 2.0 --species O H --lammps-type-map 1:O,2:H
"""
        print(help_message)
        sys.exit(0)

    args = parse_args(argv)
    configure_numpy_threads(args.omp_threads)
    compute_msd = (args.compute_water_self_msd or args.compute_msd) and not args.skip_msd
    run_unified(
        input_path=args.input_path,
        species=args.species,
        mode=args.mode,
        input_format=args.input_format,
        dt_fs=args.dt_fs,
        output_dir=args.output_dir,
        max_time_fs=args.max_time_fs,
        msd_max_time_fs=args.msd_max_time_fs,
        inset_xlim=tuple(args.inset_xlim),
        title=args.title,
        species_col=args.species_col,
        id_col=args.id_col,
        t_col=args.t_col,
        vx_col=args.vx_col,
        vy_col=args.vy_col,
        vz_col=args.vz_col,
        lammps_type_map=parse_lammps_type_map(args.lammps_type_map),
        compute_water_self=args.compute_water_self,
        water_o_species=args.water_o_species,
        water_h_species=args.water_h_species,
        water_self_target=args.water_self_target,
        mass_o=args.mass_o,
        mass_h=args.mass_h,
        n_blocks_err=args.n_blocks_err,
        compute_water_self_msd=compute_msd,
        msd_n_origins=args.msd_n_origins,
        msd_origin_mode=args.msd_origin_mode,
        msd_random_seed=args.msd_random_seed,
        msd_fit_start_fs=args.msd_fit_start_fs,
        msd_fit_end_fs=args.msd_fit_end_fs,
        max_frames=args.max_frames,
        reuse_vacf_csv=args.reuse_vacf_csv,
        compute_vacf=not args.skip_vacf,
        compute_msd=compute_msd,
    )
    elapsed = time.perf_counter() - t0
    print(f"[INFO] Total execution time: {elapsed:.2f} s")


if __name__ == "__main__":
    main()
