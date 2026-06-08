# Accelerating SCF convergence with an ML density guess

Activate the environment for this exercise:

```bash
workon density
```

## 1. Background and context


Kohn-Sham DFT finds the ground-state electron density by iteratively diagonalising the
Fock matrix $F[D] = h + V_J[D] + V_{xc}[D]$ until self-consistency. Each iteration
is expensive, so the fewer cycles needed the better.

The default starting point in most codes — including PySCF — is the **superposition of
atomic densities** (SAD): pre-computed atomic density matrices are stacked
block-diagonally. This ignores all bonding and molecular charge redistribution, so it
can be a poor starting point for molecules with complex bonding effects.

A pretrained ML model predicts the RI expansion coefficients $\{c_P\}$ of the full
molecular electron density directly from the atomic positions. Converting these to a
density matrix and using it as the initial guess starts the SCF much closer to the
solution, reducing the number of iterations to convergence.

## 2. Reading material

The recipe [**"Machine-learned dipoles and infrared spectroscopy of liquid water"**](https://atomistic-cookbook.org/examples/ml-density/ml-density.html)
computes an SCF initial guess for an organic molecule with an ML surrogate model and
uses it to speed up SCF. The key ideas:

<TODO>


👉 https://atomistic-cookbook.org/examples/ml-density/ml-density.html

By the end you should have some ideas about these questions: <TODO>

## 3. Exercise: accelerating SCF convergence with ML intial guesses

<TODO>

1. Load a set of 10 small organic molecules.
2. Visualise them with chemiscope and choose your favourite.
3. For your chosen molecule, compare the number of SCF iterations needed with:
   - the SAD initial guess (default)
   - the ML-predicted density as the initial guess
4. Repeat for 9 more organic molecules and plot the results.

### Provided code

**(i) Get the model.** The electron density model
`pet-scfbench-rho-coeffs-jfit.pt` is stored elsewhere on your virtual machine
(assuming you have run `update` before starting).

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
chemiscope.show(
    frames,
    properties={"n_atoms": [len(f) for f in frames]},
    settings=chemiscope.quick_settings(
        trajectory=True,
        map_settings={"x": {"property": "n_atoms"}, "y": {"property": "n_atoms"}},
    ),
    mode="default",
)
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

_, n_sad = run_scf(atoms, XC, BASIS, dm0=dm_sad)
print(f"SAD initial guess → {n_sad} SCF cycles")
```

ML initial guess:
```python
# Run SCF using the ML initial guess
ml_coefficients = calculator.run_model(
    atoms, {TARGET_NAME: ModelOutput(per_atom=True)}
)[TARGET_NAME]

dm_ml  = dm_from_ri_coefficients(atoms, ml_coefficients, XC, BASIS, AUXBASIS)

# Run SCF
_, n_ml = run_scf(atoms, XC, BASIS, dm0=dm_ml)
print(f"ML initial guess  → {n_ml} SCF cycles")
print(f"Speedup: {n_sad / n_ml:.1f}x fewer iterations")
```
The function `dm_from_ri_coefficients` is required to convert the
RI coefficients to the density matrix. If you are curious, open file `rho_utils.py` and
inspect this function to see how it works.


**(vii) Repeat for all 10 molecules**


If you have time, repeat the comparison for all 10 molecules and plot the result.


```python
names  = [f.get_chemical_formula() for f in frames]
n_sads, n_mls = [], []

for atoms_i in frames:
    mol_i  = atoms_to_pyscf(atoms_i, BASIS)
    mf_i   = dft.RKS(mol_i); mf_i.xc = XC
    dm_sad_i = mf_i.get_init_guess()
    _, ns = run_scf(atoms_i, XC, BASIS, dm0=dm_sad_i)

    ml_coeff_i = calculator.run_model(
        atoms_i, {TARGET_NAME: ModelOutput(per_atom=True)}
    )[TARGET_NAME]
    dm_ml_i = dm_from_ri_coefficients(atoms_i, ml_coeff_i, XC, BASIS, AUXBASIS)
    _, nm = run_scf(atoms_i, XC, BASIS, dm0=dm_ml_i)

    n_sads.append(ns); n_mls.append(nm)

# Plot results
x = np.arange(len(names))
fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)
ax.bar(
    x - 0.2,
    ...,  # TODO
    0.4,
    label="SAD",
    color="C7",
    edgecolor="white"
)
ax.bar(
    x + 0.2,
    ...,  # TODO
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

## 4. Discussion points

Once you've finished the exercise, or while you go, think about these prompts and discuss with your neighbor.

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

<TODO> Mention: model inference (model size and auxiliary basis size dependent), Fock
construction and diagonalization, integrating Ml codes with DFT codes!

</details>


<details>
<summary>💭 Instead of going through the denisty matrix, could you compute, for example, the dipole moment directly from the predicted coefficients? </summary>

Given the predicted RI coefficients, the electric dipole moment can be computed
**directly** — without constructing a density matrix or running any SCF diagonalisation:

$$\boldsymbol{\mu} = \underbrace{\sum_i Z_i\, \mathbf{R}_i}_{\text{nuclear}} \;-\; \underbrace{\sum_P c_P \int \mathrm{d}^3r\; \chi_P(\mathbf{r})\, \mathbf{r}}_{\text{electronic}}$$

The nuclear term is a simple weighted sum of atomic positions. The electronic term
requires the first spatial moment $\int \mathrm{d}^3r\; \chi_P(\mathbf{r})\, \mathbf{r}$
for each auxiliary function — a three-component vector that is evaluated numerically
on a quadrature grid.

<TODO> Finish: what lack of mathematical constraints might make this not as effective as
going through the density matrix?

</details>