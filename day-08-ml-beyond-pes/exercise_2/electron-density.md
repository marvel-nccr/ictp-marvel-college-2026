# Accelerating SCF convergence with an ML surrogate for the electron density

Activate the environment for this exercise:

```bash
cdd 08
cd exercise_2
workon metatomic
```

## 1. Background and context


Kohn-Sham DFT finds the ground-state electron density by iteratively diagonalising the
Fock matrix $F[D] = h + V_J[D] + V_{xc}[D]$ until self-consistency. Each iteration is
expensive, so the fewer cycles needed the better.

The default starting point in most codes — including PySCF — is the **superposition of
atomic densities** (SAD): pre-computed atomic density matrices are stacked
block-diagonally. This ignores all bonding and molecular charge redistribution, so it
can be a poor starting point for molecules with complex bonding effects.

A pretrained ML model predicts the RI expansion coefficients $\{c_P\}$ of the full
molecular electron density directly from the atomic positions. Converting these to a
density matrix and using it as the initial guess starts the SCF much closer to the
solution, reducing the number of iterations to convergence.

The architecture of choice that forms the ML surrogate model for the electron density in
this example is PET, an unconstrained transformer-based GNN. PET does not enforce
roto-equivariance at the level of its architecture, instead learning it via data
augmentation. Predictions of equivariant quantities -- in this case the electron density
coefficients -- are thus only approximately equivariant, though with sufficient training
time PET can learn equivariance to a small error.

## 2. Reading material

The recipe [**"ML surrogate for the electron density and derived properties"**](https://atomistic-cookbook.org/examples/ml-density/ml-density.html)
computes an SCF initial guess for an organic molecule with an ML surrogate model and
uses it to speed up SCF. The key ideas:

- **Resolution-of-identity (RI) approximation.** Rather than working with the full
  density matrix $D_{\mu\nu}$ (an $N_\text{ao}\times N_\text{ao}$ object), the density
  is re-expanded on a smaller, atom-centred *auxiliary* basis $\{\chi_P\}$:

  $$\tilde{\rho}(\mathbf{r}) = \sum_P c_P\,\chi_P(\mathbf{r})$$

  The expansion coefficients $\{c_P\}$ are found by minimising the pointwise squared
  error. The result is a compact, atom-centred, and systematically improvable
  representation — attractive properties for a machine-learning target.

- **ML model for the density.** A pretrained PET model maps atomic positions directly
  to the RI coefficients $\{c_P\}$, predicting the full molecular electron density
  without any self-consistent field calculation. Because the basis functions are
  atom-centred, the model is naturally transferable to molecules not seen during
  training.

- **Accelerated SCF.** Given predicted coefficients, a single Fock-matrix build and
  diagonalisation yields a density matrix $D_0$ already close to self-consistency.
  The SCF then needs far fewer iterations than when starting from SAD.

- **Downstream properties.** Electronic observables that are linear functionals of the
  density — such as the dipole moment — can be evaluated directly from $\tilde{\rho}$
  without running any SCF at all. The quality depends on how closely $\tilde{\rho}$
  approximates the true ground-state density.

👉 https://atomistic-cookbook.org/examples/ml-density/ml-density.html

By the end you should have some ideas about these questions: *What is the RI
approximation and why are the coefficients $\{c_P\}$ a natural ML target? Why does the
SAD initial guess work at all, and in which situations is it worst? How does the ML
density compare to SAD?*

This exercise comes in two (optional) parts: a) accelerating SCF convergence with an ML
initial guess, and b) investigating the approximate equivariance of derived energy
predictions.

## 3a. Exercise: accelerating SCF convergence with ML initial guesses

In this exercise you will apply an ML density model to a small set of organic molecules
and measure the reduction in SCF iterations compared to the SAD baseline. The workflow
mirrors the cookbook recipe but runs on a batch of molecules so you can see how the
speedup varies with the system.

1. Load a set of 5 small organic molecules from the SCFBench dataset.
2. Visualise them with chemiscope and choose your favourite.
3. For your chosen molecule, compare the number of SCF iterations needed with:
   - the SAD initial guess (default)
   - the ML-predicted density as the initial guess
4. Repeat for the rest and plot the results.
5. [Bonus] If your virtual machine is running smoothly and is fast, you can repeat the
   exercise for some larger organic molecules (stored in `molecules_larger.xyz`). While
   this should run relatively fast on a personal laptop, significant slowdowns in
   running DFT with PySCF on the virtual machines may be exhibited. Proceed with
   caution!

### Provided code

**(i) Get the model.** 

The electron density model `pet-scfbench-rho-coeffs-jfit.pt` is stored elsewhere on your
virtual machine (assuming you have run `update` before starting).

**(ii) Imports.**

The imports for the rest of the exercise:

```python
import ase.io
import chemiscope
import matplotlib.pyplot as plt
import numpy as np

from metatomic.torch import ModelOutput, load_atomistic_model
from metatomic.torch.ase_calculator import MetatomicCalculator
from pyscf import dft

from rho_utils import atoms_to_pyscf, ri_coeffs_mts_to_pyscf, dm_from_ri_coefficients, run_scf

BASIS = "def2-svp"
XC = "pbe"
AUXBASIS = "def2-universal-jfit"
MODEL = "/home/max/cosmo_models/pet-scfbench-rho-coeffs-jfit.pt"
TARGET_NAME = "mtt::rho_c_jfit_overlap"
```

**(iii) Load the model.** and create a metatomic calculator:

```python
model = load_atomistic_model(MODEL)
calculator = MetatomicCalculator(model)
```

**(iv) Load molecules and visualize.**

Load with `ASE`:

```python
frames = ase.io.read("molecules.xyz", index=":")

for i, f in enumerate(frames):
    print(f"  {i}: {f.get_chemical_formula():8s}  ({len(f)} atoms)")
```

Visualise with `chemiscope`:

```python
chemiscope.show(frames, mode="structure")
```

**(v) Choose one molecule and run SCF**

Pick the index of your favourite molecule (0-indexed).

```python
i = ...  # TODO
atoms = frames[i]
print(f"Chosen: {atoms.get_chemical_formula()} ({len(atoms)} atoms)")

chemiscope.show([atoms], mode="structure")
```

**(vi) Compare SAD baseline versus ML initial guesses**

SAD baseline:
```python
# Run SCF using the SAD initial guess
mol = atoms_to_pyscf(atoms, BASIS)
mf_sad = dft.RKS(mol)
mf_sad.xc = XC
dm_sad = mf_sad.get_init_guess()   # superposition of atomic densities

t0 = time.time()
_, n_sad = run_scf(atoms, XC, BASIS, dm0=dm_sad)
print(f"SAD initial guess → {n_sad} SCF cycles in {time.time() - t0} s")
```

ML initial guess:
```python
# Run SCF using the ML initial guess
ml_coefficients = calculator.run_model(
    atoms, {TARGET_NAME: ModelOutput(per_atom=True)}
)[TARGET_NAME]

dm_ml  = dm_from_ri_coefficients(atoms, ml_coefficients, XC, BASIS, AUXBASIS)

# Run SCF
t0 = time.time()
_, n_ml = run_scf(atoms, XC, BASIS, dm0=dm_ml)
print(f"ML initial guess  → {n_ml} SCF cycles in {time.time() - t0} s")
print(f"Speedup: {(n_ml / n_sad) * 100:.1f} of the SAD iterations")
```
The function `dm_from_ri_coefficients` is required to convert the
RI coefficients to the density matrix. If you are curious, open file `rho_utils.py` and
inspect this function to see how it works.


**(vii) Repeat for all 5 molecules**


Repeat the comparison for all 5 molecules and plot the result. What do you notice about
the reported wall times?


```python
n_sads, n_mls = [], []
for atoms_i in frames:

    print("Molecule:", atoms_i.get_chemical_formula())
    mol_i  = atoms_to_pyscf(atoms_i, BASIS)
    mf_i   = dft.RKS(mol_i); mf_i.xc = XC
    dm_sad_i = mf_i.get_init_guess()

    t0 = time.time()
    _, ns = run_scf(atoms_i, XC, BASIS, dm0=dm_sad_i)
    print(f"   SAD initial guess → {n_sad} SCF cycles in {time.time() - t0} s")


    t0 = time.time()
    ml_coeff_i = calculator.run_model(
        atoms_i, {TARGET_NAME: ModelOutput(per_atom=True)}
    )[TARGET_NAME]
    dm_ml_i = dm_from_ri_coefficients(atoms_i, ml_coeff_i, XC, BASIS, AUXBASIS)
    _, nm = run_scf(atoms_i, XC, BASIS, dm0=dm_ml_i)
    print(f"   ML initial guess  → {n_ml} SCF cycles in {time.time() - t0} s")

    n_sads.append(ns); n_mls.append(nm)
```

and plot the results:
```python
# Plot results
names  = [f.get_chemical_formula() for f in frames]
x = np.arange(len(names))
fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)
ax.bar(
    x - 0.2,
    ...,  # TODO: replace with your list of SAD iteration counts
    0.4,
    label="SAD",
    color="C7",
    edgecolor="white"
)
ax.bar(
    x + 0.2,
    ...,  # TODO: replace with your list of ML iteration counts
    0.4,
    label="ML guess",
    color="C1",
    edgecolor="white"
)
ax.set_xticks(x)
ax.set_xticklabels(names, rotation=30, ha="right")
ax.set_ylabel("SCF iterations to convergence")
ax.set_title("SAD vs ML initial guess")
ax.legend()
ax.spines[["top","right"]].set_visible(False)
```

## 3b. Exercise: investigating approximate equivariance

In this part you will use the pretrained ML model to predict the density, and from
the resulting density matrix compute the zero-shot total energy. This will be performed
for 10 rotated copies of a water molecule.

A physical density model must be **rotationally equivariant**: rotating the input
coordinates should rotate the predicted density in the same way, so that
rotation-invariant quantities such as the total energy come out identical regardless of
orientation.

The ML surrogate model you have been using in this exercise (PET) is unconstrained to
roto-equivariance on the architectural level, instead learning it implicitly through
data augmentation. Therefore, any equivariant quantities predicted by the model will
only be approximately equivariant. This also affects derived quantities such as the
total energy: physically speaking, the total energy of an isolated system is invariant
under rotations. However, as the density coefficients predicted by PET are only
approximately equivariant, the resulting total energy will be only approximately
invariant.

In this exercise 

Here you will measure the residual numerical error by computing the no-SCF DFT energy —
the energy functional $E[\tilde\rho]$ evaluated at the ML-predicted density — for 10
different orientations of a water molecule and checking how much it varies.

### Provided code

**(i) Additional import.**

```python
import ase
from scipy.spatial.transform import Rotation
```

**(ii) Build a water molecule with `ase.Atoms`.**

```python
d_OH  = 0.9584                         # O–H bond length, Å
theta = np.deg2rad(104.45 / 2)         # half the H–O–H angle

water = ase.Atoms(
    symbols  = ["O", "H", "H"],
    positions= [
        [ 0.0,  0.0,                    0.0                   ],
        [ 0.0,  d_OH * np.sin(theta),  -d_OH * np.cos(theta)  ],
        [ 0.0, -d_OH * np.sin(theta),  -d_OH * np.cos(theta)  ],
    ],
)
water.positions -= water.get_center_of_mass()   # centre at origin
```

**(iii) Helper: evaluate total DFT energy at a given density matrix (no SCF).**

```python
def energy_from_dm(atoms, dm, xc, basis):
    """Total DFT energy E[D] without any SCF iteration."""
    mol        = atoms_to_pyscf(atoms, basis)
    mf         = dft.RKS(mol)
    mf.xc      = xc
    mf.verbose = 0
    return mf.energy_tot(dm=dm)
```

**(iv) Zero-shot energy comparison.**

Before looking at rotations, compute the no-SCF energy from the SAD and ML initial
density matrices for the original molecule and compare both to the converged SCF
energy. This tells you how close each initial guess sits to the true ground state in
terms of energy prediction.

After this, for the rotational analysis, we will only look at the deviations of the
rotation-dependent energies against the mean for each initial guess separately. This
will help us decouple the accuracy of the energy from its rotational dependence.

```python
# SAD initial guess for the original molecule
mol_ref  = atoms_to_pyscf(water, BASIS)
mf_ref   = dft.RKS(mol_ref)
mf_ref.xc = XC
mf_ref.verbose = 0
dm_sad_ref = mf_ref.get_init_guess()

# ML prediction for the original molecule
ml_coeff_ref = calculator.run_model(
    water, {TARGET_NAME: ModelOutput(per_atom=True)}
)[TARGET_NAME]
dm_ml_ref = dm_from_ri_coefficients(water, ml_coeff_ref, XC, BASIS, AUXBASIS)

# Converged SCF energy (reference)
mf_conv, _ = run_scf(water, XC, BASIS, dm0=dm_sad_ref)
e_conv = mf_conv.e_tot

# No-SCF energies
e_sad_ref = energy_from_dm(water, dm_sad_ref, XC, BASIS)
e_ml_ref  = energy_from_dm(water, dm_ml_ref,  XC, BASIS)

print(f"Converged SCF energy : {e_conv:.6f} Ha  (reference)")
print(f"SAD zero-shot energy : {e_sad_ref:.6f} Ha  (error = {(e_sad_ref - e_conv)*1e3:+.2f} mHa)")
print(f"ML  zero-shot energy : {e_ml_ref:.6f} Ha  (error = {(e_ml_ref  - e_conv)*1e3:+.2f} mHa)")
```

**(v) Define 10 rotation matrices.**

We generate 10 Cartesian rotation matrices evenly spaced from 0° to 360° around a
general axis (not aligned with any molecular symmetry element).

```python
angles_deg = np.linspace(0, 360, 11)[:-1]   # 0°, 36°, 72°, …, 324°

axis = np.array([1.0, 1.0, 0.5])
axis /= np.linalg.norm(axis)                 # normalise

rotation_matrices = [
    Rotation.from_rotvec(np.deg2rad(a) * axis).as_matrix()
    for a in angles_deg
]
```

**(vi) Loop over rotations — compute no-SCF energies for SAD and ML.**

For each rotated geometry we evaluate the no-SCF energy for both the SAD and the ML
initial density matrix.

```python
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
print(f"\nSAD energy std : {(energies_sad.std()) * 1e3:.4f} mHa")
print(f"ML  energy std : {(energies_ml.std())  * 1e3:.4f} mHa")
```

**(vii) Plot.**

Both curves are shown relative to their own mean so they share the same axis scale.

```python
fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True, dpi=120)
ax.plot(
    angles_deg,
    (energies_sad - energies_sad.mean()) * 1e3,
    "s-",
    color="C7",
    lw=2,
    ms=7,
    label="SAD"
)
ax.plot(
    angles_deg,
    (energies_ml  - energies_ml.mean())  * 1e3,
    "o-", 
    color="C1", 
    lw=2, 
    ms=7, 
    label="ML"
)
ax.axhline(0, color="0.5", lw=0.8, ls="--")
ax.set_xlabel("Rotation angle / degrees")
ax.set_ylabel(r"$E - \langle E \rangle$ / mHa")
ax.set_title("Equivariance test: no-SCF energy vs rotation")
ax.legend(frameon=False)
ax.spines[["top", "right"]].set_visible(False)
```


## 4. Discussion points

Once you've finished the exercise, or while you go, think about these prompts and
discuss with your neighbor.

### SCF acceleration

<details>
<summary>💭 Why does the ML guess help?</summary>

The ML model has learned to predict the full molecular density, including the effects
of bonding and charge polarisation, directly from the atomic positions. Its prediction
is much closer to the converged density than the SAD superposition, so the SCF needs
fewer iterations to reach self-consistency. 

The speedup tends to be larger for molecules where SAD is worst (strong polarisation,
conjugation), and due to the ~ O(3) complexity of DFT may be the systems where the
benefit of an ML initial guess in most pronounced.

</details>

<details>
<summary>💭 Beyond the number of SCF cycles saved, what impacts whether the actual wall time of the SCF calculation is decreased?</summary>

The raw iteration count is only one factor. The actual wall-time speedup also depends on:

- **Model inference cost.** Evaluating the ML model on the given geometry takes time
  that grows with model size and auxiliary basis size. For small molecules the saving of
  cycles may not be benficial.
- **Per-iteration cost.** Each SCF step involves Fock matrix construction and a
  diagonalisation (scaling as $\mathcal{O}(N^3)$ in the basis size). Saving even one
  or two iterations is worth more for large, expensive systems.
- **Software integration overhead.** Passing data between the ML framework and the DFT
  code introduces latency that does not appear in the raw iteration count.

In practice the ML initial guess pays off most clearly for larger molecules, where
each SCF cycle is expensive and the model inference cost is a small fraction of the
total.

</details>


<details>
<summary>💭 Instead of going through the density matrix, could you compute, for example, the dipole moment directly from the predicted coefficients? </summary>

Given the predicted RI coefficients, the electric dipole moment can be computed
**directly** — without constructing a density matrix or running any SCF diagonalisation:

$$\boldsymbol{\mu} = \underbrace{\sum_i Z_i\, \mathbf{R}_i}_{\text{nuclear}} \;-\; \underbrace{\sum_P c_P \int \mathrm{d}^3r\; \chi_P(\mathbf{r})\, \mathbf{r}}_{\text{electronic}}$$

The nuclear term is a simple weighted sum of atomic positions. The electronic term
requires the first spatial moment $\int \mathrm{d}^3r\; \chi_P(\mathbf{r})\, \mathbf{r}$
for each auxiliary function — a three-component vector evaluated numerically on a
quadrature grid.

</details>

<details>
<summary>💭 What lack of mathematical contraints might make this not as effective as going through the density matrix? </summary>

However, the predicted density $\tilde{\rho}$ satisfies no strict mathematical
constraints: it is not guaranteed to integrate to the correct number of electrons
(particle-number conservation), nor to be $N$-representable. Errors in $\tilde{\rho}$
propagate directly into the dipole without any self-correcting mechanism. Going through
a density matrix — via one Fock diagonalisation, as `dm_from_ri_coefficients` does —
enforces the correct electron count and a physical density, typically improving the
accuracy of derived properties at the additional cost of a diagonalization.

### Approximate equivariance

</details>

<details>
<summary>💭 Why do you expect the energy deviations as a function of rotation curve for the SAD guess to be flat? Why is it not? </summary>

They should be by construction, as the atomic densities are speherically symmetric and
thus invariant to rotations. However in practice this isn't exactly the case, due to
numerics and perhaps the "egg box" problem, where the discretization of the DFT grids
plays a role.

</details>


<details>
<summary>💭 What would you expect the deviations curves to look like if you used an equivariant architecture, such as MACE? </summary>

One would expect that they should be a lot flatter, as the predicted coefficients would
transform equivariantly under rotation of the system, by design.

</details>


<details>
<summary>💭 Computing the standard deviation of the derived energies as a function of rotation is fairly straightforward. How could you assess the 'equivariance error' of the coefficients themselves? </summary>

<TODO>: One would need to back-transform the predictions on the rotated systems to the
original reference frame. As the density coefficients are defined on a basis of
spherical harmonics, one would need to use Wigner-D matrices for these rotations. For
accurate metrics, one would need to evaluate the model on a sufficiently converged
quadrature that samples the O(3) group, too.

For more details on equivariance and how unconstrained models learn them, read our
paper:

"How unconstrained machine learning models learn physical symmetries", M. Domina,
J.W. Abbott, P. Pegolo, F. Bigi, M. Ceriotti, https://arxiv.org/pdf/2603.24638

</details>