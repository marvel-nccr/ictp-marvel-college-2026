#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


def configure_numpy_threads(n_threads):
    if n_threads is None:
        return
    s = str(int(n_threads))
    os.environ["OMP_NUM_THREADS"] = s
    os.environ["OPENBLAS_NUM_THREADS"] = s
    os.environ["MKL_NUM_THREADS"] = s
    os.environ["VECLIB_MAXIMUM_THREADS"] = s
    os.environ["NUMEXPR_NUM_THREADS"] = s


def parse_lammps_type_map(text):
    if text is None:
        return None
    out = {}
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        k, v = item.split(":")
        out[int(k.strip())] = v.strip()
    return out


def format_time_tag(fs_value):
    if fs_value is None:
        return ""
    fs_value = float(fs_value)
    if fs_value >= 1e6:
        return f"_{fs_value/1e6:.3f}ns".replace(".", "p")
    if fs_value >= 1e3:
        return f"_{fs_value/1e3:.3f}ps".replace(".", "p")
    return f"_{fs_value:.1f}fs".replace(".", "p")


def p1(x):
    return x


def p2(x):
    return 0.5 * (3.0 * x * x - 1.0)


def minimum_image(dr, box):
    return dr - box * np.round(dr / box)


def parse_frame_arrays(fh):
    line = fh.readline()
    if not line:
        return None
    if not line.startswith("ITEM: TIMESTEP"):
        raise ValueError("Unexpected dump format: expected 'ITEM: TIMESTEP'")

    timestep = int(fh.readline().strip())

    line = fh.readline().strip()
    if line != "ITEM: NUMBER OF ATOMS":
        raise ValueError("Unexpected dump format: expected 'ITEM: NUMBER OF ATOMS'")
    n_atoms = int(fh.readline().strip())

    line = fh.readline().strip()
    if not line.startswith("ITEM: BOX BOUNDS"):
        raise ValueError("Unexpected dump format: expected 'ITEM: BOX BOUNDS'")

    bounds = []
    for _ in range(3):
        lo, hi, *_ = fh.readline().split()
        bounds.append((float(lo), float(hi)))
    origin = np.array([lo for lo, hi in bounds], dtype=float)
    box = np.array([hi - lo for lo, hi in bounds], dtype=float)

    line = fh.readline().strip()
    if not line.startswith("ITEM: ATOMS"):
        raise ValueError("Unexpected dump format: expected 'ITEM: ATOMS ...'")
    columns = line.split()[2:]
    col_idx = {c: i for i, c in enumerate(columns)}

    required_any = (
        {"x", "y", "z"} <= set(columns) or
        {"xu", "yu", "zu"} <= set(columns)
    )
    if not required_any:
        raise ValueError("Trajectory must contain x/y/z or xu/yu/zu")
    for c in ("id", "mol", "type"):
        if c not in col_idx:
            raise ValueError(f"Trajectory must contain '{c}' column")

    data = np.empty((n_atoms, len(columns)), dtype=float)
    for i in range(n_atoms):
        data[i] = np.fromstring(fh.readline(), sep=" ", dtype=float)

    return {
        "timestep": timestep,
        "n_atoms": n_atoms,
        "box": box,
        "origin": origin,
        "columns": columns,
        "col_idx": col_idx,
        "data": data,
    }


def extract_wrapped_xyz_from_frame(frame):
    data = frame["data"]
    col = frame["col_idx"]
    box = frame["box"]
    origin = frame["origin"]

    if {"x", "y", "z"} <= set(col):
        return data[:, [col["x"], col["y"], col["z"]]]

    xyz = data[:, [col["xu"], col["yu"], col["zu"]]]
    return origin[None, :] + np.mod(xyz - origin[None, :], box[None, :])


def build_water_topology_from_first_frame(frame, type_map, water_o_species="O", water_h_species="H"):
    data = frame["data"]
    col = frame["col_idx"]

    atom_id = data[:, col["id"]].astype(int)
    mol = data[:, col["mol"]].astype(int)
    atype = data[:, col["type"]].astype(int)

    species = np.array([type_map.get(int(t), None) for t in atype], dtype=object)
    if np.any(species == None):  # noqa: E711
        bad = sorted(set(int(t) for t, s in zip(atype, species) if s is None))
        raise ValueError(f"Unmapped LAMMPS types: {bad}")

    order = np.argsort(atom_id)
    atom_id = atom_id[order]
    mol = mol[order]
    species = species[order]
    xyz_idx = order

    unique_mol = np.unique(mol)
    o_idx = []
    h1_idx = []
    h2_idx = []
    mol_ids = []

    for m in unique_mol:
        mask = mol == m
        idx_local = np.where(mask)[0]
        idx_O = idx_local[species[idx_local] == water_o_species]
        idx_H = idx_local[species[idx_local] == water_h_species]
        if len(idx_O) != 1 or len(idx_H) != 2:
            continue
        o_idx.append(xyz_idx[idx_O[0]])
        h1_idx.append(xyz_idx[idx_H[0]])
        h2_idx.append(xyz_idx[idx_H[1]])
        mol_ids.append(m)

    if not mol_ids:
        raise ValueError("No valid water molecules found in first frame")

    return {
        "mol_ids": np.asarray(mol_ids, dtype=int),
        "o_idx": np.asarray(o_idx, dtype=int),
        "h1_idx": np.asarray(h1_idx, dtype=int),
        "h2_idx": np.asarray(h2_idx, dtype=int),
        "n_mol": len(mol_ids),
    }


def frame_to_orientation_vectors(frame, topo, target="dipole"):
    xyz = extract_wrapped_xyz_from_frame(frame)
    box = frame["box"]

    O = xyz[topo["o_idx"]]
    H1w = xyz[topo["h1_idx"]]
    H2w = xyz[topo["h2_idx"]]

    r1 = minimum_image(H1w - O, box[None, :])
    r2 = minimum_image(H2w - O, box[None, :])

    if target == "dipole":
        v = r1 + r2
    elif target == "bisector":
        u1 = r1 / np.linalg.norm(r1, axis=1)[:, None]
        u2 = r2 / np.linalg.norm(r2, axis=1)[:, None]
        v = u1 + u2
    elif target == "oh":
        v = r1
    elif target == "hh":
        v = minimum_image(H2w - H1w, box[None, :])
    else:
        raise ValueError(f"Unknown target: {target}")

    nrm = np.linalg.norm(v, axis=1)
    bad = nrm <= 0.0
    if np.any(bad):
        raise ValueError("Found zero-norm orientation vector")
    return v / nrm[:, None]


def build_orientation_timeseries(
    input_path,
    type_map,
    water_o_species="O",
    water_h_species="H",
    target="dipole",
    max_frames=None,
):
    with open(input_path, "r", encoding="utf-8") as fh:
        first = parse_frame_arrays(fh)
        if first is None:
            raise ValueError("No frames found")

        topo = build_water_topology_from_first_frame(
            first, type_map, water_o_species=water_o_species, water_h_species=water_h_species
        )

        U_list = []
        timesteps = []

        U0 = frame_to_orientation_vectors(first, topo, target=target)
        U_list.append(U0)
        timesteps.append(first["timestep"])

        n_read = 1
        while True:
            if max_frames is not None and n_read >= max_frames:
                break
            frame = parse_frame_arrays(fh)
            if frame is None:
                break
            U = frame_to_orientation_vectors(frame, topo, target=target)
            U_list.append(U)
            timesteps.append(frame["timestep"])
            n_read += 1

    U = np.stack(U_list, axis=0)
    timesteps = np.asarray(timesteps, dtype=float)
    return U, timesteps, topo["mol_ids"]


def choose_origins(n_frames, max_lag_frames, mode="uniform", n_origins=50, seed=None):
    n_valid = n_frames - max_lag_frames
    if n_valid <= 0:
        raise ValueError("Not enough frames for requested max lag")
    all_origins = np.arange(n_valid, dtype=int)

    if mode == "all":
        return all_origins
    if n_origins is None or n_origins <= 0:
        raise ValueError("n_origins must be positive for uniform/random")
    if n_origins >= n_valid:
        return all_origins
    if mode == "uniform":
        return np.linspace(0, n_valid - 1, n_origins, dtype=int)
    if mode == "random":
        rng = np.random.default_rng(seed)
        return np.sort(rng.choice(all_origins, size=n_origins, replace=False))
    raise ValueError(f"Unknown origin mode: {mode}")


def rotational_correlation(U, dt_fs, max_time_fs=None, origin_mode="uniform", n_origins=50, random_seed=None):
    n_frames, _, _ = U.shape

    if dt_fs is None:
        raise ValueError("--dt-fs is required")

    dt_fs = float(dt_fs)
    if max_time_fs is None:
        max_lag_frames = n_frames - 1
    else:
        max_lag_frames = min(int(math.floor(float(max_time_fs) / dt_fs)), n_frames - 1)

    origins = choose_origins(
        n_frames=n_frames,
        max_lag_frames=max_lag_frames,
        mode=origin_mode,
        n_origins=n_origins,
        seed=random_seed,
    )

    n_lags = max_lag_frames + 1
    c1 = np.zeros(n_lags, dtype=float)
    c2 = np.zeros(n_lags, dtype=float)

    for o in origins:
        dots = np.einsum("mci,mci->mc", U[o:o + n_lags], U[o][None, :, :])
        dots = np.clip(dots, -1.0, 1.0)
        c1 += np.mean(p1(dots), axis=1)
        c2 += np.mean(p2(dots), axis=1)

    c1 /= float(len(origins))
    c2 /= float(len(origins))

    time_fs = np.arange(n_lags, dtype=float) * dt_fs
    return time_fs, c1, c2, origins


def estimate_tau_integral(time_fs, corr):
    integrate = getattr(np, "trapezoid", None)
    if integrate is None:
        integrate = np.trapz
    return float(integrate(corr, x=time_fs))


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Rotational relaxation tool for water from LAMMPS trajectories")
    p.add_argument("input_path", help="LAMMPS trajectory file (.lammpstrj/.dump)")
    p.add_argument("--input-format", choices=["auto", "lammps"], default="auto")
    p.add_argument("--dt-fs", type=float, default=None, help="Physical time step between saved frames in fs")
    p.add_argument("--output-dir", default=None, help="Where to write outputs (default: input file directory)")
    p.add_argument("--max-time-fs", type=float, default=None, help="Maximum lag time for rotational correlation (fs)")
    p.add_argument("--title", default=None)
    p.add_argument("--species", nargs="*", default=None, help="Species labels, e.g. --species O H")
    p.add_argument("--n_species", type=int, default=None)
    p.add_argument("--lammps-type-map", default=None, help="Map LAMMPS atom types to species, e.g. 1:O,2:H")
    p.add_argument("--water-o-species", default="O")
    p.add_argument("--water-h-species", default="H")
    p.add_argument("--rot-target", nargs="+", choices=["dipole", "bisector", "oh", "hh", "all"], default=["dipole"])
    p.add_argument("--rot-n-origins", type=int, default=50)
    p.add_argument("--rot-origin-mode", choices=["uniform", "all", "random"], default="uniform")
    p.add_argument("--rot-random-seed", type=int, default=None)
    p.add_argument("--max-frames", type=int, default=None)
    p.add_argument("--omp-threads", type=int, default=None)

    args = p.parse_args(argv)

    if args.species is None or len(args.species) == 0:
        if args.n_species is None:
            args.species = ["O", "H"]
        else:
            args.species = [f"S{i+1}" for i in range(args.n_species)]

    targets = args.rot_target
    if "all" in targets:
        args.rot_target = ["dipole", "oh", "hh"]
    else:
        # remove duplicates while preserving order
        seen = set()
        cleaned = []
        for t in targets:
            if t not in seen:
                cleaned.append(t)
                seen.add(t)
        args.rot_target = cleaned

    return args


def run_one_target(args, target):
    output_dir = Path(args.output_dir) if args.output_dir else Path(args.input_path).resolve().parent
    output_dir.mkdir(parents=True, exist_ok=True)

    type_map = parse_lammps_type_map(args.lammps_type_map)

    U, timesteps, mol_ids = build_orientation_timeseries(
        input_path=args.input_path,
        type_map=type_map,
        water_o_species=args.water_o_species,
        water_h_species=args.water_h_species,
        target=target,
        max_frames=args.max_frames,
    )

    time_fs, c1, c2, origins = rotational_correlation(
        U=U,
        dt_fs=args.dt_fs,
        max_time_fs=args.max_time_fs,
        origin_mode=args.rot_origin_mode,
        n_origins=args.rot_n_origins,
        random_seed=args.rot_random_seed,
    )

    tau1_int_fs = estimate_tau_integral(time_fs, c1)
    tau2_int_fs = estimate_tau_integral(time_fs, c2)

    tag = format_time_tag(args.max_time_fs)
    csv_path = output_dir / f"rot_relax_{target}{tag}.csv"
    summary_path = output_dir / f"rot_relax_summary_{target}{tag}.csv"

    df = pd.DataFrame({
        "time_fs": time_fs,
        "time_ps": time_fs / 1e3,
        "C1": c1,
        "C2": c2,
    })
    df.to_csv(csv_path, index=False)

    summary = {
        "input_path": str(args.input_path),
        "target": target,
        "n_frames": U.shape[0],
        "n_molecules": U.shape[1],
        "dt_fs": float(args.dt_fs),
        "max_time_fs": float(time_fs[-1]),
        "origin_mode": args.rot_origin_mode,
        "n_origins_used": int(len(origins)),
        "tau1_integral_fs": tau1_int_fs,
        "tau1_integral_ps": tau1_int_fs / 1e3,
        "tau2_integral_fs": tau2_int_fs,
        "tau2_integral_ps": tau2_int_fs / 1e3,
    }
    pd.DataFrame([summary]).to_csv(summary_path, index=False)

    print(f"[INFO] Target           : {target}")
    print(f"[INFO] Wrote           : {csv_path}")
    print(f"[INFO] Wrote           : {summary_path}")
    print(f"[INFO] tau1_integral   = {tau1_int_fs/1e3:.4f} ps")
    print(f"[INFO] tau2_integral   = {tau2_int_fs/1e3:.4f} ps")


def main(argv=None):
    t0 = time.perf_counter()

    if argv is None and len(sys.argv) == 1:
        print(
            """
Usage: python rot_relax.py <traj.lammpstrj> [options]
-----------------------------------------------------
Compute rotational correlation functions C1(t), C2(t)
and associated relaxation times for water.

Examples:
  python rot_relax.py traj_prod_long_recursive.lammpstrj \\
    --input-format lammps \\
    --dt-fs 10 \\
    --species O H \\
    --lammps-type-map 1:O,2:H \\
    --rot-target dipole

  python rot_relax.py traj_prod_long_recursive.lammpstrj \\
    --input-format lammps \\
    --dt-fs 10 \\
    --species O H \\
    --lammps-type-map 1:O,2:H \\
    --rot-target dipole oh hh

  python rot_relax.py traj_prod_long_recursive.lammpstrj \\
    --input-format lammps \\
    --dt-fs 10 \\
    --species O H \\
    --lammps-type-map 1:O,2:H \\
    --rot-target all
"""
        )
        sys.exit(0)

    args = parse_args(argv)
    configure_numpy_threads(args.omp_threads)

    for target in args.rot_target:
        run_one_target(args, target)

    elapsed = time.perf_counter() - t0
    print(f"[INFO] Total execution time: {elapsed:.2f} s")


if __name__ == "__main__":
    main()
