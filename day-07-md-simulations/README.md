# MD simulations: Hands-on session

**Tutors:** Federica Troni, Edrick Solis Gonzalez

---

This folder contains three exercises on molecular dynamics simulations of water with LAMMPS and Jupyter notebooks.

## Setup

Activate the Python environment before running anything:

```bash
workon metatomic
```

Then launch JupyterLab from the terminal inside this folder:

```bash
jupyter lab
```

Open the notebook for the exercise you want to work on from the JupyterLab file browser.

Select the kernel: Python (metatomic) in the notebook.

---

## Exercise 1 — Thermostat Tuning in NVT Equilibration

**Goal.** Study how the thermostat damping time `tdamp` and simulation lenght affects NVT equilibration for SPC/Fw and TIP3P water.

---

## Exercise 2 — Diffusion from VACF and MSD

**Goal.** Compute the self-diffusion coefficient of water from simulations using both the VACF and the MSD in NVE.

---

## Exercise 3 — Thermostat Effects on Water Observables

**Goal.** Compare how different thermostats affect structural and dynamical observables of water, including temperature, radial distribution functions, diffusion, and rotational relaxation.
