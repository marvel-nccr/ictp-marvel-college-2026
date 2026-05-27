#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


def parse_lammps_type_map(text: str | None) -> dict[int, str]:
    if text is None:
        return {}
    out: dict[int, str] = {}
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        atom_type, species = item.split(":")
        out[int(atom_type.strip())] = species.strip()
    return out


def parse_frame_arrays(handle):
    line = handle.readline()
    if not line:
        return None
    if not line.startswith("ITEM: TIMESTEP"):
        raise ValueError("Unexpected dump format: expected 'ITEM: TIMESTEP'")

    timestep = int(handle.readline().strip())

    if handle.readline().strip() != "ITEM: NUMBER OF ATOMS":
        raise ValueError("Unexpected dump format: expected 'ITEM: NUMBER OF ATOMS'")
    n_atoms = int(handle.readline().strip())

    if not handle.readline().strip().startswith("ITEM: BOX BOUNDS"):
        raise ValueError("Unexpected dump format: expected 'ITEM: BOX BOUNDS'")

    bounds = []
    for _ in range(3):
        lo, hi, *_ = handle.readline().split()
        bounds.append((float(lo), float(hi)))
    origin = np.array([lo for lo, _ in bounds], dtype=float)
    box = np.array([hi - lo for lo, hi in bounds], dtype=float)

    atom_header = handle.readline().strip()
    if not atom_header.startswith("ITEM: ATOMS"):
        raise ValueError("Unexpected dump format: expected 'ITEM: ATOMS ...'")

    columns = atom_header.split()[2:]
    col_idx = {name: i for i, name in enumerate(columns)}
    if "type" not in col_idx:
        raise ValueError("Trajectory must contain the 'type' column")
    if {"x", "y", "z"} <= set(columns):
        xyz_cols = ("x", "y", "z")
    elif {"xu", "yu", "zu"} <= set(columns):
        xyz_cols = ("xu", "yu", "zu")
    else:
        raise ValueError("Trajectory must contain x/y/z or xu/yu/zu")

    data = np.empty((n_atoms, len(columns)), dtype=float)
    for i in range(n_atoms):
        data[i] = np.fromstring(handle.readline(), sep=" ", dtype=float)

    xyz = data[:, [col_idx[c] for c in xyz_cols]]
    if xyz_cols == ("xu", "yu", "zu"):
        xyz = origin[None, :] + np.mod(xyz - origin[None, :], box[None, :])

    return {
        "timestep": timestep,
        "box": box,
        "columns": columns,
        "col_idx": col_idx,
        "data": data,
        "xyz": xyz,
    }


def minimum_image(delta: np.ndarray, box: np.ndarray) -> np.ndarray:
    return delta - box * np.round(delta / box)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compute O-O, O-H, and H-H radial distribution functions from a LAMMPS dump.")
    parser.add_argument("trajectory", type=Path, help="LAMMPS trajectory in dump format")
    parser.add_argument("--output-dir", type=Path, default=Path("."), help="Directory for rdf_oo.csv")
    parser.add_argument("--lammps-type-map", default="1:O,2:H", help="Map atom types to species labels")
    parser.add_argument("--bins", type=int, default=250, help="Number of histogram bins")
    parser.add_argument("--r-max", type=float, default=None, help="Maximum distance in angstrom")
    parser.add_argument("--max-frames", type=int, default=None, help="Use at most this many frames")
    return parser


def main() -> int:
    args = build_arg_parser().parse_args()
    type_map = parse_lammps_type_map(args.lammps_type_map)
    if not type_map:
        raise ValueError("A valid --lammps-type-map is required")

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    histograms: dict[str, np.ndarray] = {}
    edges = None
    n_frames = 0
    volume_sum = 0.0
    counts_by_species: dict[str, int] | None = None
    masks_by_species: dict[str, np.ndarray] | None = None

    with args.trajectory.open(encoding="utf-8") as handle:
        while True:
            frame = parse_frame_arrays(handle)
            if frame is None:
                break

            atom_types = frame["data"][:, frame["col_idx"]["type"]].astype(int)
            if masks_by_species is None:
                species = np.array([type_map.get(atom_type, "") for atom_type in atom_types], dtype=object)
                masks_by_species = {
                    "O": species == "O",
                    "H": species == "H",
                }
                counts_by_species = {name: int(np.count_nonzero(mask)) for name, mask in masks_by_species.items()}
                if counts_by_species["O"] < 2 or counts_by_species["H"] < 2:
                    raise ValueError("Need at least two oxygen atoms and two hydrogen atoms to compute g(r)")
                r_max = args.r_max if args.r_max is not None else 0.5 * float(np.min(frame["box"])) * 0.98
                edges = np.linspace(0.0, r_max, args.bins + 1)
                histograms = {
                    "oo": np.zeros(args.bins, dtype=float),
                    "oh": np.zeros(args.bins, dtype=float),
                    "hh": np.zeros(args.bins, dtype=float),
                }

            oxygen_xyz = frame["xyz"][masks_by_species["O"]]
            hydrogen_xyz = frame["xyz"][masks_by_species["H"]]
            box = frame["box"]
            delta_oo = oxygen_xyz[:, None, :] - oxygen_xyz[None, :, :]
            delta_oo = minimum_image(delta_oo, box[None, None, :])
            distances_oo = np.linalg.norm(delta_oo[np.triu_indices(counts_by_species["O"], k=1)], axis=1)
            histograms["oo"] += np.histogram(distances_oo, bins=edges)[0]

            delta_hh = hydrogen_xyz[:, None, :] - hydrogen_xyz[None, :, :]
            delta_hh = minimum_image(delta_hh, box[None, None, :])
            distances_hh = np.linalg.norm(delta_hh[np.triu_indices(counts_by_species["H"], k=1)], axis=1)
            histograms["hh"] += np.histogram(distances_hh, bins=edges)[0]

            delta_oh = oxygen_xyz[:, None, :] - hydrogen_xyz[None, :, :]
            delta_oh = minimum_image(delta_oh, box[None, None, :])
            distances_oh = np.linalg.norm(delta_oh.reshape(-1, 3), axis=1)
            histograms["oh"] += np.histogram(distances_oh, bins=edges)[0]
            volume_sum += float(np.prod(box))
            n_frames += 1

            if args.max_frames is not None and n_frames >= args.max_frames:
                break

    if n_frames == 0 or not histograms or edges is None or counts_by_species is None:
        raise FileNotFoundError(f"No readable frames found in {args.trajectory}")

    shell_volumes = (4.0 / 3.0) * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    r_centers = 0.5 * (edges[1:] + edges[:-1])
    mean_volume = volume_sum / n_frames
    densities = {
        "O": counts_by_species["O"] / mean_volume,
        "H": counts_by_species["H"] / mean_volume,
    }

    pair_specs = {
        "oo": ("O", "O", 0.5, "rdf_oo.csv"),
        "oh": ("O", "H", 1.0, "rdf_oh.csv"),
        "hh": ("H", "H", 0.5, "rdf_hh.csv"),
    }
    print(f"[INFO] Frames used for RDF: {n_frames}")
    for pair_name, (species_a, species_b, prefactor, filename) in pair_specs.items():
        n_a = counts_by_species[species_a]
        rho_b = densities[species_b]
        ideal_counts = prefactor * n_a * rho_b * shell_volumes
        g_r = histograms[pair_name] / (n_frames * ideal_counts)
        out_path = output_dir / filename
        with out_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["r_A", "g_r", "counts", "n_frames", "n_a", "n_b", "rho_b_A-3", "pair"])
            for r_value, g_value, count in zip(r_centers, g_r, histograms[pair_name].astype(int)):
                writer.writerow([r_value, g_value, count, n_frames, counts_by_species[species_a], counts_by_species[species_b], rho_b, pair_name])
        print(f"[INFO] Wrote {out_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
