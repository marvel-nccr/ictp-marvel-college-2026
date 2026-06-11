# Uncertainty quantification & FlashMD: Hands-on session

**Tutors**: Filippo and Matthias

---

This session has two independent parts, each in its own subfolder:

- **`uncertainty/`** — quantifying the prediction error of a machine-learning
  interatomic potential (PET-MAD) and propagating it from single-point energies all
  the way to thermodynamic averages.
- **`nc-and-flashmd/`** — using non-conservative models and long-stride molecular
  dynamics with the universal **FlashMD** model.


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

## Part 2 — Machine learning in molecular dynamics beyond potentials

Folder: [`nc-and-flashmd/`](nc-and-flashmd/).

To get started, open a terminal and run the following commands:

```bash
cd nc-and-flashmd/  # cd into the uncertainty folder
workon metatomic  # activate the `metatomic` environment
```

This moves into this directory and activates the `metatomic` Python environment. Then
launch JupyterLab from the same terminal:

```bash
jupyter lab
```

Open the notebook for the part you want to work on from the JupyterLab file browser and
select the **Python (metatomic)** kernel.

**Goal.** Learn how to run accelerated MD with machine-learning methods, first using
the PET-MAD universal interatomic potential with non-conservative forces and multiple
time stepping, then using the universal FlashMD model for molecular dynamics prediction
using long strides.

1. **Non-conservative forces** — Accelerate MD compared to MLIPs by leveraging faster
   force predictions.
2. **FlashMD** — Accelerate MD even further with models predicting simulations in long
   strides.
3. **Interactive viewers** — Interact with MD systems in your browser and explore!

---

## References

- **PET-MAD**: Mazitov *et al.*, 2025 — <https://arxiv.org/abs/2503.14118>
- **Last-layer prediction rigidity (LLPR) uncertainties**: Bigi *et al.*, 2024 —
  <https://arxiv.org/abs/2403.02251>
- **Shallow-ensemble UQ**: Kellner & Ceriotti, *Uncertainty quantification by direct
  propagation of shallow ensembles*, <https://arxiv.org/abs/2402.16621>
- **Error propagation / calibration**: Imbalzano *et al.*, 2021 —
  <https://arxiv.org/abs/2011.08828>
- **Correct use of non-conservative forces in simulations**: Bigi *et al.*, 2024 —
  <https://proceedings.mlr.press/v267/bigi25a.html>
- **FlashMD**: Bigi *et al.*, 2025 — <https://papers.nips.cc/paper_files/paper/2025/hash/7ffade093764c9a59a777c3cbe346b97-Abstract-Conference.html>
