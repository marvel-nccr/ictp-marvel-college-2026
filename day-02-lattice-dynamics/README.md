# Density functional perturbation theory: Hands-on session

**Tutors**: Alberto Carta, Alfredo Fiorentino, Libor Vojacek

## Environment setup

Activate the working environment before running any exercise:

```bash
workon qe
```

## Exercises

### Example 1: Phonons in silicon bulk (`01_Si_QE_a`, `01_Si_QE_b`)

**a) Phonons at the Gamma point**
Compute the vibrational frequencies of bulk Si at the zone centre using Quantum ESPRESSO (`ph.x`). This potentially gives direct access to the Raman/IR-active modes and is the simplest starting point for a phonon calculation.

**b) Phonon dispersion**
Extend the calculation to a uniform **q**-point grid, Fourier-interpolate the interatomic force constants with `q2r.x` and `matdyn.x`, and plot the full phonon band structure along a high-symmetry path.

---

### Example 2: Phonons in a polar material — AlAs (`02_AlAs_QE_a`, `02_AlAs_QE_b`)

**a) Phonons at the Gamma point**
Repeat the zone-centre calculation for AlAs. Because AlAs is a polar semiconductor, non-analytic corrections (Born effective charges and the macroscopic dielectric tensor) are essential to reproduce the LO-TO splitting.

**b) Phonon dispersion**
Compute and plot the full dispersion for AlAs, including the non-analytic corrections along the high-symmetry path to correctly capture the LO-TO splitting near Γ.

---

### Example 3: Relaxation and phonons with a MLIP (`03_phonons_mlip`)

> **Environment**: this part uses a different Python environment. Deactivate the QE environment first, then activate the MLIP one:
> ```bash
> conda deactivate qe
> ```

Open the notebook in JupyterLab, either launch it from the desktop shortcut, or from the terminal:

```bash
jupyter-lab 03_phonons_mlip/day2_mlip_ase_eos_relax_phonons.ipynb
```

Use the **UPET PET-MAD** machine-learning interatomic potential (via ASE) to:

1. Scan the equation of state (EOS) and fit a Birch–Murnaghan curve for Si and AlAs.
2. Relax the unit cell with `FrechetCellFilter` + BFGS and compare the relaxed volume against the EOS minimum.
3. Compute phonon band structures from finite displacements for both materials along the Γ → X → W → X → Γ → L path and compare with the DFPT results.


# Example 4: Anharmonicity and thermal lattice expansion

Here we will compute the thermal expansion of Gold using the stochastic self-consistent harmonic approximation.
We need to activate the `scha` environment:

```bash
conda activate sscha
```

Inside the `04_thermal_lattice_expansion` directory we have the results of a DFPT calculation with quantum espresso
of gold, with the dynamical matrices saved as `gold_harmonic_dyn_X` (X goes from 1 to 13, indicating the index of the 
irreducible q points in a `4x4x4` mesh of q-points).

1. Learn how to generate a random ensemble distributed according to the harmonic dynamical matrix
2. Run the first SSCHA free energy minimization, to quantify anharmonicity at fixed volume.
3. Compare the anharmonic and harmonic phonon dispersion of gold at 300 K
4. Run the variable cell relaxation as a function of temperature
5. Fit the data and compute the lattice thermal expansion coefficient of gold.

For any of this task, we provide the commented scripts that can inspected. 
