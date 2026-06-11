"""Exercise 3b: Investigating approximate equivariance of the ML density model.

Run from the exercise_2 directory with the correct environment:

    workon density
    python solution_3b.py

Produces: equivariance_test.pdf
"""

import ase
import matplotlib.pyplot as plt
import numpy as np
from metatomic.torch import ModelOutput, load_atomistic_model
from metatomic.torch.ase_calculator import MetatomicCalculator
from pyscf import dft
from scipy.spatial.transform import Rotation

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

# ── build water molecule ───────────────────────────────────────────────────────
d_OH  = 0.9584                    # O–H bond length, Å
theta = np.deg2rad(104.45 / 2)   # half the H–O–H angle

water = ase.Atoms(
    symbols   = ["O", "H", "H"],
    positions = [
        [ 0.0,  0.0,                    0.0                   ],
        [ 0.0,  d_OH * np.sin(theta),  -d_OH * np.cos(theta)  ],
        [ 0.0, -d_OH * np.sin(theta),  -d_OH * np.cos(theta)  ],
    ],
)
water.positions -= water.get_center_of_mass()


# ── helper ─────────────────────────────────────────────────────────────────────
def energy_from_dm(atoms, dm, xc, basis):
    """Total DFT energy E[D] without any SCF iteration."""
    mol        = atoms_to_pyscf(atoms, basis)
    mf         = dft.RKS(mol)
    mf.xc      = xc
    mf.verbose = 0
    return mf.energy_tot(dm=dm)


# ── zero-shot energy comparison ────────────────────────────────────────────────
print("\n── Zero-shot energy comparison ──")

mol_ref        = atoms_to_pyscf(water, BASIS)
mf_ref         = dft.RKS(mol_ref)
mf_ref.xc      = XC
mf_ref.verbose = 0
dm_sad_ref     = mf_ref.get_init_guess()

ml_coeff_ref = calculator.run_model(
    water, {TARGET_NAME: ModelOutput(per_atom=True)}
)[TARGET_NAME]
dm_ml_ref = dm_from_ri_coefficients(water, ml_coeff_ref, XC, BASIS, AUXBASIS)

mf_conv, _ = run_scf(water, XC, BASIS, dm0=dm_sad_ref)
e_conv     = mf_conv.e_tot
e_sad_ref  = energy_from_dm(water, dm_sad_ref, XC, BASIS)
e_ml_ref   = energy_from_dm(water, dm_ml_ref,  XC, BASIS)

print(f"Converged SCF energy : {e_conv:.6f} Ha  (reference)")
print(f"SAD zero-shot energy : {e_sad_ref:.6f} Ha  (error = {(e_sad_ref - e_conv)*1e3:+.2f} mHa)")
print(f"ML  zero-shot energy : {e_ml_ref:.6f} Ha  (error = {(e_ml_ref  - e_conv)*1e3:+.2f} mHa)")

# ── rotation matrices ──────────────────────────────────────────────────────────
angles_deg = np.linspace(0, 360, 11)[:-1]   # 0°, 36°, 72°, …, 324°

axis = np.array([1.0, 1.0, 0.5])
axis /= np.linalg.norm(axis)

rotation_matrices = [
    Rotation.from_rotvec(np.deg2rad(a) * axis).as_matrix()
    for a in angles_deg
]

# ── rotation loop ──────────────────────────────────────────────────────────────
print("\n── Rotation test ──")
energies_sad = []
energies_ml  = []

for R, angle in zip(rotation_matrices, angles_deg):
    atoms_rot = water.copy()
    atoms_rot.set_positions(atoms_rot.get_positions() @ R.T)

    # SAD initial guess for this rotation
    mol_rot        = atoms_to_pyscf(atoms_rot, BASIS)
    mf_rot         = dft.RKS(mol_rot)
    mf_rot.xc      = XC
    mf_rot.verbose = 0
    dm_sad_rot     = mf_rot.get_init_guess()

    # ML prediction for this rotation
    ml_coefficients = calculator.run_model(
        atoms_rot, {TARGET_NAME: ModelOutput(per_atom=True)}
    )[TARGET_NAME]
    dm_ml_rot = dm_from_ri_coefficients(atoms_rot, ml_coefficients, XC, BASIS, AUXBASIS)

    e_sad = energy_from_dm(atoms_rot, dm_sad_rot, XC, BASIS)
    e_ml  = energy_from_dm(atoms_rot, dm_ml_rot,  XC, BASIS)
    energies_sad.append(e_sad)
    energies_ml.append(e_ml)
    print(f"  θ = {angle:6.1f}°  →  E_SAD = {e_sad:.6f}  E_ML = {e_ml:.6f} Ha")

energies_sad = np.array(energies_sad)
energies_ml  = np.array(energies_ml)
print(f"\nSAD energy std : {energies_sad.std() * 1e3:.4f} mHa")
print(f"ML  energy std : {energies_ml.std()  * 1e3:.4f} mHa")

# ── plot ───────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True, dpi=120)
ax.plot(angles_deg, (energies_sad - energies_sad.mean()) * 1e3,
        "s-", color="C7", lw=2, ms=7, label="SAD")
ax.plot(angles_deg, (energies_ml  - energies_ml.mean())  * 1e3,
        "o-", color="C1", lw=2, ms=7, label="ML")
ax.axhline(0, color="0.5", lw=0.8, ls="--")
ax.set_xlabel("Rotation angle / degrees")
ax.set_ylabel(r"$E - \langle E \rangle$ / mHa")
ax.set_title("Equivariance test: no-SCF energy vs rotation")
ax.legend(frameon=False)
ax.spines[["top", "right"]].set_visible(False)

out = "equivariance_test.pdf"
fig.savefig(out)
print(f"\nSaved → {out}")
