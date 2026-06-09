"""Exercise 3: NMR-driven structure determination with ShiftML3.

Run from the exercise_3 directory with the correct environment:

    workon nmr
    python solutions/solution.py

Produces three lollipop plots (PDF):
    lollipop_default.pdf       — per-structure intercept calibration (all protons)
    lollipop_no_acidic.pdf     — per-structure intercept, acidic proton excluded
    lollipop_global_cal.pdf    — global calibration (all protons included)
"""

import shutil
from collections import OrderedDict

import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from platformdirs import user_cache_path
from shiftml.ase import ShiftML
from sklearn.metrics import root_mean_squared_error

# ── settings ───────────────────────────────────────────────────────────────────
NUM_H_PER_MOLECULE = 31

# Acidic proton: key "5" sits at 15.37 ppm — well above the 10–13 ppm typical
# range, indicating a strongly deshielded exchangeable proton.
ACIDIC_PROTON_KEY = "5"

# Global calibration parameters (PBE/GIPAW fitted against experimental 1H shifts
# for a diverse benchmark of molecular crystals).
GLOBAL_SLOPE     = -0.9024
GLOBAL_INTERCEPT =  28.05

# ── copy model to cache and initialise ShiftML3 ────────────────────────────────
print("Setting up ShiftML3...")
shutil.copytree(
    "/home/max/cosmo_models/shiftml",
    user_cache_path() / "shiftml",
    dirs_exist_ok=True,
)
calculator = ShiftML("ShiftML3")

# ── load candidate structures ──────────────────────────────────────────────────
frames = read("azd_molecular_crystals.xyz", ":")
print(f"Loaded {len(frames)} candidate structures  "
      f"({len(frames[0])} atoms per unit cell)")

# ── experimental 1H chemical shifts (assigned, ppm) ───────────────────────────
assigned_experimental_shifts = OrderedDict({
    "1":                          6.92,
    "2":                          8.69,
    "3":                          9.01,
    "4":                          8.47,
    "5":                         15.37,   # ← acidic proton
    "6":                          7.73,
    "7":                          9.64,
    "8":                          2.90,
    "9":                          1.78,
    "10,11":                      1.88,
    "12":                         1.80,
    "13":                         1.60,
    "14":                         0.44,
    "15":                         1.54,
    "16,17":                      1.88,
    "18,19":                      0.80,
    "20":                         1.00,
    "21,22":                      1.74,
    "23,24,25,26,27,28,29,30,31": 0.73,
})
experimental_shifts = np.array(list(assigned_experimental_shifts.values()))

# ── helper functions ───────────────────────────────────────────────────────────

def assign_shieldings(per_h_shieldings, assigned_experimental_shifts):
    """Pick one molecule's H shieldings from the unit cell and average over
    symmetry-equivalent groups listed in assigned_experimental_shifts."""
    per_mol = per_h_shieldings.reshape(NUM_H_PER_MOLECULE, -1)[:, 0]
    out = []
    for atom_string in assigned_experimental_shifts.keys():
        idx = [int(s) - 1 for s in atom_string.split(",")]
        out.append(per_mol[idx].mean())
    return np.array(out)


def calibrated_rmse(shieldings_per_candidate, experimental_shifts,
                    slope=-1.0, intercept=None):
    """RMSE for each candidate after linear calibration σ → δ."""
    rmses = []
    for sigmas in shieldings_per_candidate:
        b = intercept
        if b is None:
            b = np.mean(experimental_shifts) - slope * np.mean(sigmas)
        predicted_shifts = slope * sigmas + b
        rmses.append(root_mean_squared_error(predicted_shifts, experimental_shifts))
    return np.array(rmses)


def make_lollipop_plot(frames, rmse_gipaw, rmse_sml, title=None):
    """Lollipop plot comparing GIPAW and ShiftML3 RMSE vs experiment."""
    candidate_idx = np.arange(len(frames))
    dx = 0.18
    fig, ax = plt.subplots(figsize=(8.5, 5.2), constrained_layout=True, dpi=120)

    band_lo, band_hi = 0.33 - 0.16, 0.33 + 0.16
    ax.axhspan(band_lo, band_hi, color="0.85", alpha=0.7, zorder=0,
               label=r"DFT vs experiment noise floor ($0.33 \pm 0.16$ ppm)")

    ax.vlines(candidate_idx - dx, 0, rmse_gipaw, color="C0", lw=2.5, alpha=0.85, zorder=2)
    ax.vlines(candidate_idx + dx, 0, rmse_sml,   color="C1", lw=2.5, alpha=0.85, zorder=2)
    ax.scatter(candidate_idx - dx, rmse_gipaw, color="C0", label="GIPAW (DFT reference)",
               s=60, edgecolor="white", lw=0.9, zorder=3)
    ax.scatter(candidate_idx + dx, rmse_sml,   color="C1", label="ShiftML3",
               s=60, edgecolor="white", lw=0.9, zorder=3)

    ax.set_xlabel("Candidate structure index", fontsize=13)
    ax.set_ylabel(r"$^1$H shift RMSE / ppm", fontsize=13)
    ax.set_xticks(candidate_idx[::2])
    ax.set_xlim(-0.7, len(frames) - 0.3)
    ax.set_ylim(0, max(rmse_sml.max(), rmse_gipaw.max()) * 1.15)
    ax.grid(axis="y", color="0.92", lw=0.6, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper left", frameon=False, fontsize=12)
    if title:
        ax.set_title(title, fontsize=13)
    return fig


# ── run ShiftML3 inference ─────────────────────────────────────────────────────
print("\nRunning ShiftML3 inference (this may take ~30 s)...")
shieldings_sml   = []
shieldings_gipaw = []

for frame in frames:
    is_h = frame.get_atomic_numbers() == 1
    sml   = calculator.get_cs_iso(frame).ravel()[is_h]
    gipaw = frame.arrays["CS"][is_h]
    shieldings_sml.append(assign_shieldings(sml,   assigned_experimental_shifts))
    shieldings_gipaw.append(assign_shieldings(gipaw, assigned_experimental_shifts))

shieldings_sml   = np.array(shieldings_sml)
shieldings_gipaw = np.array(shieldings_gipaw)
print("Inference complete.")

# ── plot 1: default per-structure calibration ──────────────────────────────────
print("\n── Plot 1: per-structure intercept calibration (all protons) ──")
rmse_sml   = calibrated_rmse(shieldings_sml,   experimental_shifts)
rmse_gipaw = calibrated_rmse(shieldings_gipaw, experimental_shifts)

best = int(np.argmin(rmse_sml))
print(f"  Best ShiftML3 : candidate #{best}  (RMSE = {rmse_sml[best]:.3f} ppm)")
print(f"  Best GIPAW    : candidate #{int(np.argmin(rmse_gipaw))}  "
      f"(RMSE = {rmse_gipaw.min():.3f} ppm)")

fig = make_lollipop_plot(
    frames, rmse_gipaw, rmse_sml,
    title="Per-structure calibration (all protons)",
)
fig.savefig("lollipop_default.pdf")
print("  Saved → lollipop_default.pdf")

# ── plot 2: exclude acidic proton ──────────────────────────────────────────────
print(f"\n── Plot 2: acidic proton '{ACIDIC_PROTON_KEY}' "
      f"({assigned_experimental_shifts[ACIDIC_PROTON_KEY]:.2f} ppm) excluded ──")

slice_idxs = [
    i for i, key in enumerate(assigned_experimental_shifts.keys())
    if key != ACIDIC_PROTON_KEY
]
print(f"  Using {len(slice_idxs)} / {len(assigned_experimental_shifts)} peaks")

rmse_sml   = calibrated_rmse(shieldings_sml[:,   slice_idxs],
                              experimental_shifts[slice_idxs])
rmse_gipaw = calibrated_rmse(shieldings_gipaw[:, slice_idxs],
                              experimental_shifts[slice_idxs])

best = int(np.argmin(rmse_sml))
print(f"  Best ShiftML3 : candidate #{best}  (RMSE = {rmse_sml[best]:.3f} ppm)")
print(f"  Best GIPAW    : candidate #{int(np.argmin(rmse_gipaw))}  "
      f"(RMSE = {rmse_gipaw.min():.3f} ppm)")

fig = make_lollipop_plot(
    frames, rmse_gipaw, rmse_sml,
    title=f"Per-structure calibration (key '{ACIDIC_PROTON_KEY}' excluded)",
)
fig.savefig("lollipop_no_acidic.pdf")
print("  Saved → lollipop_no_acidic.pdf")

# ── plot 3: global calibration ─────────────────────────────────────────────────
print(f"\n── Plot 3: global calibration  "
      f"(a = {GLOBAL_SLOPE}, b = {GLOBAL_INTERCEPT} ppm, all protons) ──")

rmse_sml   = calibrated_rmse(shieldings_sml,   experimental_shifts,
                              slope=GLOBAL_SLOPE, intercept=GLOBAL_INTERCEPT)
rmse_gipaw = calibrated_rmse(shieldings_gipaw, experimental_shifts,
                              slope=GLOBAL_SLOPE, intercept=GLOBAL_INTERCEPT)

best = int(np.argmin(rmse_sml))
print(f"  Best ShiftML3 : candidate #{best}  (RMSE = {rmse_sml[best]:.3f} ppm)")
print(f"  Best GIPAW    : candidate #{int(np.argmin(rmse_gipaw))}  "
      f"(RMSE = {rmse_gipaw.min():.3f} ppm)")

fig = make_lollipop_plot(
    frames, rmse_gipaw, rmse_sml,
    title=fr"Global calibration ($a={GLOBAL_SLOPE}$, $b={GLOBAL_INTERCEPT}$ ppm)",
)
fig.savefig("lollipop_global_cal.pdf")
print("  Saved → lollipop_global_cal.pdf")
