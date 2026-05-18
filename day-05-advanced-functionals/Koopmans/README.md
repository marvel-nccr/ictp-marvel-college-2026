# Koopmans functionals: Hands-on session

**Tutors**: Nicola Colonna and Edward Linscott

In this session we will use the [`koopmans`](https://koopmans-functionals.org) package to run Koopmans spectral-functional calculations on two systems: the ozone molecule and bulk zinc oxide.

## Background

Standard (semi-local) density-functional theory describes ground-state properties remarkably well, but it is poorly suited to *spectral* properties — ionisation potentials, electron affinities, band gaps, photoemission spectra. The underlying problem is that the Kohn–Sham eigenvalues do not have a rigorous physical meaning (except for the HOMO in exact KS-DFT, via Janak's theorem), and approximate functionals violate the piecewise-linearity of the total energy with respect to fractional electron number.

Koopmans functionals impose piecewise linearity at the level of each individual variational orbital. The result is an orbital-density-dependent functional in which the eigenvalues become approximations of charged-excitation energies — *i.e.* they predict the IP, EA, and band structure directly.

Any given Koopmans calculation requires several ingredients:

1. **A choice of variational orbitals.** For molecules, the Kohn–Sham orbitals; for periodic solids, **maximally-localised Wannier functions (MLWFs)**
2. **A set of screening parameters $\alpha_i$**, one per variational orbital, that account for orbital relaxation. These can be computed in two ways:
   - the **ΔSCF** method, which performs explicit constrained $N\pm1$ calculations in a supercell;
   - the **DFPT** method, which uses density-functional perturbation theory in the primitive cell.

## Exercises

- [`01-ozone/`](01-ozone/) — Compute the ionisation potential and electron affinity of ozone using KI with the **ΔSCF** screening method.
- [`02-zno-dfpt/`](02-zno-dfpt/) — Compute the band structure of bulk ZnO using KI with the **DFPT** screening method.

> **Note**
>
> Before starting, check that `koopmans` is installed by running `koopmans --help` in your terminal. If you don't see a help message, ask a tutor.

## Further reading

- The `koopmans` documentation: <https://koopmans-functionals.org>
- E. Linscott *et al.*, *koopmans: An Open-Source Package for Accurately and Efficiently Predicting Spectral Properties with Koopmans Functionals*, J. Chem. Theory Comput. **19**, 7097 (2023)
- N. Colonna *et al.*, *Koopmans Spectral Functionals in Periodic Boundary Conditions*, J. Chem. Theory Comput. **18**, 5435 (2022)
- R. De Gennaro *et al.*, *Bloch's theorem in orbital-density-dependent functionals: Band structures from Koopmans spectral functionals*, Phys. Rev. B **106**, 035106 (2022)
