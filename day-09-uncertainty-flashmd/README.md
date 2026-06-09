# Uncertainty quantification & FlashMD: Hands-on session

**Tutors**: Filippo and Matthias

---

This session has two independent parts, each in its own subfolder:

- **`uncertainty/`** — quantifying the prediction error of a machine-learning
  interatomic potential (PET-MAD) and propagating it from single-point energies all
  the way to thermodynamic averages.
- **`flashmd/`** — long-stride molecular dynamics with the universal **FlashMD** model.


---

## Part 1 — Uncertainty quantification with PET-MAD

Folder: [`uncertainty/`](uncertainty/) (see its `README.md` for details).

To get started, open a terminal and run the following commands:

```bash
cd uncertainty  # cd into the uncertainty folder
workon metatomic  # activate the `metatomic` environment
```

This moves into this directory and activates the `metatomic` Python environment. Then
launch JupyterLab from the same terminal:

```bash
jupyter lab
```

Open the notebook for the part you want to work on from the JupyterLab file browser and
select the **Python (metatomic)** kernel.

**Goal.** Learn three increasingly demanding ways of attaching error bars to the
predictions of an ML potential or simulations with the potential.
All using PET-MAD's built-in uncertainty quantification
(last-layer prediction rigidity and committee/shallow-ensemble models):

1. **Dataset uncertainties** — estimate per-structure uncertainties on a validation set
   and check them against the true errors with a parity plot.
2. **Vacancy formation energies** — propagate ensemble uncertainties through a simple
   function of energies (a vacancy formation energy in aluminum).
3. **Error propagation in MD** — propagate committee uncertainties to a thermodynamic
   average (the radial distribution function of liquid water) via reweighting in i-PI.

---

## Part 2 — Long-stride MD with FlashMD

Folder: [`flashmd/`](flashmd/).

_Content to be added by the FlashMD tutors._

---

## References

- **PET-MAD**: Mazitov *et al.*, 2025 — <https://arxiv.org/abs/2503.14118>
- **Last-layer prediction rigidity (LLPR) uncertainties**: Bigi *et al.*, 2024 —
  <https://arxiv.org/abs/2403.02251>
- **Shallow-ensemble UQ**: Kellner & Ceriotti, *Uncertainty quantification by direct
  propagation of shallow ensembles*, <https://arxiv.org/abs/2402.16621>
- **Error propagation / calibration**: Imbalzano *et al.*, 2021 —
  <https://arxiv.org/abs/2011.08828>
