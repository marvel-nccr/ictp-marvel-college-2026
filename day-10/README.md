# Embedded Cluster Expansion with pyeCE: Hands-on session

**Tutors**: Yann Lorris Mueller & Anirudh Raju Natarajan

In this session you will build, train, and use an **embedded cluster expansion** (eCE) model — a machine-learned, on-lattice atomistic model with a neural network backend — to study a six-component (senary) BCC alloy containing Cr, Mo, Nb, Ta, V, and W.

The session is divided into three parts:

1. **Construct an eCE model** — define the crystal structure and chemical degrees of freedom via a PRIM file, then specify the interaction cutoffs, basis, and neural network architecture.
2. **Train the eCE model** — explore a database of DFT formation energies, fit the model using ladder training, and validate it by comparing DFT and eCE convex hulls.
3. **Compute finite-temperature thermodynamic properties** — run canonical Monte Carlo simulations to study short-range order, and evaluate the free energy via thermodynamic integration.

## Getting started

Navigate to this directory and launch Jupyter:

```bash
cdd 10
jupyter lab tutorial.ipynb
```

Work through `tutorial.ipynb` from top to bottom. Each section contains worked examples followed by tasks for you to complete. Solutions are available as collapsible blocks beneath each task.

> [!NOTE]
>
> If your `model_BCC.pth` did not build or train correctly, a pre-trained model is provided in `Files/model_BCC.pth`. You can use it for all remaining exercises.

> [!NOTE]
>
> If you are unable to run the canonical Monte Carlo simulations, pre-computed results are provided in `Files/results_canonical_Mo-Ta.json` and `Files/results_canonical.json` and can be used directly for the analysis tasks.

## Files provided

Pre-computed data and models are in the `Files/` directory:

| File                                 | Description                                                                                     |
|--------------------------------------|-------------------------------------------------------------------------------------------------|
| `Files/model_BCC.pth`                | Pre-trained BCC senary eCE model                                                                |
| `Files/database.pkl`                 | DFT training database — over 3600 symmetrically unique BCC arrangements of Cr, Mo, Nb, Ta, V, W |
| `Files/results_canonical_Mo-Ta.json` | Canonical MC results for the equimolar Mo–Ta binary at 1000 K                                   |
| `Files/results_canonical.json`       | Canonical MC results for the equiatomic senary alloy at 2000 K and 1000 K                       |
| `Files/cooldown.dat`                 | Enthalpy vs. temperature from a cooling simulation of the equiatomic senary alloy               |
