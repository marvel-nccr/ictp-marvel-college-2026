# Uncertainty quantification with PET-MAD

**Tutors**: List here

---

This part shows three ways of computing and propagating errors on the outputs of an
ML-potential-driven simulation, using **PET-MAD** and its built-in uncertainty
quantification (UQ). PET-MAD v1.5.0 ships with both *last-layer prediction rigidity*
(LLPR) uncertainties and a *committee / shallow ensemble* of last-layer models, so no
extra training is needed.

## Getting started

From a terminal:

```bash
cdd 9          # cd into day-09 and activate the `metatomic` environment
cd uncertainty  # cd into the uncertainty folder
jupyter lab    # launch JupyterLab
```

Open the notebook in this folder and select the **Python (metatomic)** kernel.

The example uses the PET-MAD model (extra-small variant, for speed),
given the 


---

## Exercise 1 — Uncertainties on a dataset

**Goal.** Estimate per-structure uncertainties for single-point energy predictions on a
validation set and check whether they are trustworthy.

We run PET-MAD over a subset of the MAD-1.5 r2SCAN validation set, request the
`energy_uncertainty` output alongside the energy, and compare the predicted
uncertainties to the *true* errors (predicted minus reference energy). The comparison is
visualized as a **parity plot of uncertainties vs. errors** on log–log axes — a
well-calibrated model has points scattered around the diagonal.

Questions to think about:
- Why are larger structures more uncertain? - is this useful to know?
- How does the uncertainty correlate with the error? - is this useful to know?
- How can uncertainties be used to build better datasets?


---

## Exercise 2 — Uncertainties in vacancy formation energies

**Goal.** Propagate ensemble uncertainties through a simple *function* of energies.

Using an aluminum crystal, we compute the **vacancy formation energy** by tracking the
ensemble energies (`energy_ensemble`) at three stages: the pristine bulk, immediately
after removing one atom, and after relaxing the defected cell. Because every ensemble
member gives its own value, we obtain a mean **and** a standard deviation for the
vacancy formation energy directly.

Questions to think about:
- Is it justified to only relax with the mean of the ensemble, and then compute ensemble dissociation-energies on the mean-relaxed structure? What would be the alternative?
- How cames that the predicted total potential energies of the bulk vary significantly across the ensemble, but the vacancy formation energy from different ensemble members quantitatively disagree much less? - think about error cancellation and the fact that the vacancy formation energy is a *difference* of energies.


---

## Exercise 3 — Error propagation in MD

**Goal.** Propagate committee uncertainties to a *thermodynamic average* computed over an
MD trajectory.

We run a short i-PI simulation of a box of 32 water molecules with a committee model,
compute the **radial distribution functions** (H–H and O–O) frame by frame with ASE,
and then use i-PI's `i-pi-committee-reweight` tool to reweight the RDFs by the committee
energies. This yields a mean RDF and envelopping bands showing how
potential-energy uncertainty maps onto a structural observable. The driver's
`uncertainty_threshold` option also flags frames whose atomic-energy uncertainty exceeds
a chosen value (in eV/atom).

> The in-session run uses only a few hundred MD steps for speed, so the RDFs and their
> error bands will look noisy. Increasing the number of steps (e.g. to 10000) gives
> much smoother results.

Questions to think about:
- The uncertainties of the O-O RDFs are quite elevated, even when moving to the precomputed values that have been run over extended trajectories. 
How could the uncertainty be reduced? - is this a problem of the sampling, or of the model? How could you test this?

---

## References

- PET-MAD: Mazitov *et al.*, 2025 — <https://arxiv.org/abs/2503.14118>
- LLPR uncertainties: Bigi *et al.*, 2024   <https://arxiv.org/abs/2403.02251>
- Shallow-ensemble UQ: Kellner & Ceriotti, 
  <https://arxiv.org/abs/2402.16621>
- Calibration & error propagation: Imbalzano *et al.*, 2021 —
  <https://arxiv.org/abs/2011.08828>
- Full recipe: *Uncertainty Quantification with PET-MAD* in the
  [atomistic cookbook](https://atomistic-cookbook.org).
