"""Exercise 3a: Accelerating SCF convergence with ML initial guesses.

Run from the exercise_2 directory with the correct environment:

    workon density
    python solution_3a.py

Produces: scf_acceleration.pdf
"""

import time

import ase.io
import matplotlib.pyplot as plt
import numpy as np
from metatomic.torch import ModelOutput, load_atomistic_model
from metatomic.torch.ase_calculator import MetatomicCalculator
from pyscf import dft

from rho_utils import atoms_to_pyscf, dm_from_ri_coefficients, run_scf

# ── settings ───────────────────────────────────────────────────────────────────
BASIS       = "def2-svp"
XC          = "pbe"
AUXBASIS    = "def2-universal-jfit"
MODEL       = "/home/max/cosmo_models/pet-scfbench-rho-coeffs-jfit.pt"
TARGET_NAME = "mtt::rho_c_jfit_overlap"

# ── load model ─────────────────────────────────────────────────────────────────
print("Loading model...")
model      = load_atomistic_model(MODEL)
calculator = MetatomicCalculator(model)

# ── load molecules ─────────────────────────────────────────────────────────────
frames = ase.io.read("molecules.xyz", index=":")
print(f"\nLoaded {len(frames)} molecules:")
for i, f in enumerate(frames):
    print(f"  {i}: {f.get_chemical_formula():8s}  ({len(f)} atoms)")

# ── SAD vs ML comparison for each molecule ─────────────────────────────────────
names  = [f.get_chemical_formula() for f in frames]
n_sads, n_mls = [], []

for atoms_i in frames:
    formula = atoms_i.get_chemical_formula()
    print(f"\n── {formula} ──")

    mol_i     = atoms_to_pyscf(atoms_i, BASIS)
    mf_i      = dft.RKS(mol_i)
    mf_i.xc   = XC
    dm_sad_i  = mf_i.get_init_guess()

    t0 = time.time()
    _, ns = run_scf(atoms_i, XC, BASIS, dm0=dm_sad_i)
    t_sad = time.time() - t0
    print(f"  SAD  → {ns:2d} cycles  ({t_sad:.1f} s)")

    t0 = time.time()
    ml_coeff_i = calculator.run_model(
        atoms_i, {TARGET_NAME: ModelOutput(per_atom=True)}
    )[TARGET_NAME]
    dm_ml_i = dm_from_ri_coefficients(atoms_i, ml_coeff_i, XC, BASIS, AUXBASIS)
    _, nm = run_scf(atoms_i, XC, BASIS, dm0=dm_ml_i)
    t_ml = time.time() - t0
    print(f"  ML   → {nm:2d} cycles  ({t_ml:.1f} s)")

    n_sads.append(ns)
    n_mls.append(nm)

# ── plot ───────────────────────────────────────────────────────────────────────
x   = np.arange(len(names))
fig, ax = plt.subplots(figsize=(9, 4), constrained_layout=True, dpi=120)
ax.bar(x - 0.2, n_sads, 0.4, label="SAD",      color="C7", edgecolor="white")
ax.bar(x + 0.2, n_mls,  0.4, label="ML guess", color="C1", edgecolor="white")
ax.set_xticks(x)
ax.set_xticklabels(names, rotation=30, ha="right")
ax.set_ylabel("SCF iterations to convergence")
ax.set_title("SAD vs ML initial guess")
ax.legend(frameon=False)
ax.spines[["top", "right"]].set_visible(False)

out = "scf_acceleration.pdf"
fig.savefig(out)
print(f"\nSaved → {out}")
