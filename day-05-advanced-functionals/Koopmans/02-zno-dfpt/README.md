# Exercise 2: Band structure of bulk ZnO with KI and DFPT

**Tutors**: Nicola Colonna and Edward Linscott

In this exercise you will compute the KI band structure of wurtzite ZnO. There are two important differences with the ozone calculation:

1. **The variational orbitals are now maximally-localised Wannier functions (MLWFs)** instead of Kohn–Sham orbitals.
2. **The screening parameters are computed via DFPT** instead of via constrained ΔSCF calculations in a supercell.

This exercise is adapted from [tutorial 3](https://koopmans-functionals.org/en/latest/tutorials/tutorial_3.html) of the `koopmans` documentation.

> **Note**
>
> A full DFPT screening calculation for ZnO takes longer than the time available in this session. We have therefore pre-computed the screening parameters $\alpha_i$ and provided them in `zno.json` via the `alpha_guess` field, with `calculate_alpha` set to `false`. The KI band structure is computed using these pre-computed values.

## Files provided

- `zno.json` — the input file for the calculation
- `plot_bandstructures.py` — a `python` script to plot the LDA and KI band structures side by side

---

## Problem 1: Understanding the input file

Open `zno.json` and locate the `workflow` block:

```json
{
  "workflow": {
    "task": "singlepoint",
    "functional": "ki",
    "base_functional": "lda",
    "method": "dfpt",
    "init_orbitals": "mlwfs",
    "calculate_alpha": false,
    "alpha_guess": [[0.358, 0.364, ...]],
    "pseudo_library": "PseudoDojo/0.4/LDA/SR/standard/upf",
    "gb_correction": true,
    "eps_inf": 5.3
  }
}
```

### Part A

Compare this block against the one from Exercise 1. Identify and explain each difference:

- `task: singlepoint` (we'll change this in Problems 2 and 3)
- `base_functional: lda` instead of the default PBE
- `method: dfpt` instead of `dscf`
- `init_orbitals: mlwfs` instead of `kohn-sham`
- `calculate_alpha: false` with an explicit `alpha_guess`
- `gb_correction` and `eps_inf`

### Part B

The `gb_correction` is the *Gygi–Baldereschi* correction, which accounts for the spurious interaction of an orbital with its periodic images. It requires the high-frequency dielectric constant `eps_inf` of the material. Why is this finite-size correction necessary for solids but not for ozone?

### Part C

Inspect the `calculator_parameters.w90` block. This sets up the Wannierization. Why are the projections split into multiple sub-lists (so-called "blocks")?

---

## Problem 2: The Wannierization

Before doing any Koopmans physics, we need a good set of Wannier functions. Whether they are "good" depends on whether the band structure they interpolate matches the explicit DFT band structure.

### Part A

Change `"task": "singlepoint"` to `"task": "wannierize"` in `zno.json`, then run

```bash
koopmans run zno.json | tee zno_wannierize.md
```

This task performs an SCF calculation, an NSCF calculation, a Wannierization of each block of bands, and then an explicit DFT band-structure calculation on the same _**k**_-path for comparison. The workflow generates a plot `zno_bandstructure.png` that overlays the interpolated and the explicit band structures.

### Part B

Inspect `zno_bandstructure.png`. Do the interpolated bands lie on top of the explicit ones? If not, where do the largest discrepancies appear, and in which energy window?

### Part C

The Wannierization of the empty (Zn 4s) manifold uses two **disentanglement** keywords, which you can find in the `w90` block of `zno.json`:

```json
"dis_froz_max": 14.5,
"dis_win_max": 17.0
```

`dis_win_max` is the upper bound of the disentanglement window (within which Bloch states are mixed to construct the Wannier functions); `dis_froz_max` is the upper bound of the frozen window (within which the Bloch states are *exactly* preserved). Try varying these values by ±1 eV and see how the interpolated band structure changes.

> **Tip**
>
> When you run the same workflow more than once, `koopmans` will reuse intermediate results from previous runs by default. If you want to start from scratch, set `"from_scratch": true` in the `engine` block of `zno.json`.

---

## Problem 3: Running the KI calculation

Now change `"task"` back to `"singlepoint"` and re-run

```bash
koopmans run zno.json | tee zno_ki.md
```

This time the workflow proceeds beyond the Wannierization to actually compute the KI band structure. Open `zno_ki.md` and identify the following steps:

- a **Wannierization** block (the same one you just ran),
- a **`wann2kc`** step, which converts the Wannier90 files into a format readable by `kcw.x`,
- a **`ham`** step, which constructs the Koopmans Hamiltonian and computes the band structure.

> **Note**
>
> If we had set `"calculate_alpha": true` we would also see a **`screen`** step between `wann2kc` and `ham`, in which the $\alpha_i$ values are computed via DFPT. With the pre-computed values from `alpha_guess`, this step is skipped.

### Part A

Contrast this workflow with the one you saw for ozone. In the DFPT workflow, where does the cost of computing the screening parameters go (when we *do* compute them)? Why is this expected to scale much better with system size than ΔSCF?

### Part B

What is the role of `wann2kc`? Why is this conversion step needed when going from a Wannier90 output to a Koopmans calculation?

---

## Problem 4: Plotting and analysing the band structure

The `singlepoint` task generates `zno_bandstructure.png` automatically. For a more comprehensive plot, run

```bash
python plot_bandstructures.py
```

This produces `zno_bandstructures.png`, comparing the LDA band structure against the KI@LDA one and labelling the band gap.

### Part A

Inspect `zno_bandstructures.png`. What is the KI@LDA band gap? Compare it against

- the LDA band gap (read off the same plot),
- the experimental value of 3.4 eV,
- and the converged KI@LDA result of 3.62 eV from Colonna *et al.* (2022)[^Colonna2022].

### Part B

The KI calculation does not change the *shape* of the LDA band structure dramatically — it primarily acts as a *scissor* that opens up the gap. Looking at the bands away from the gap, can you identify any features that *are* modified beyond a rigid shift?

### Part C

By construction, the Koopmans correction depends on which orbital you are looking at — orbitals with $\alpha_i$ close to 1 are corrected strongly, while orbitals with $\alpha_i$ close to 0 are barely shifted. Look at the `alpha_guess` in `zno.json` and predict which manifolds of bands (Zn 3d–O 2p valence vs. Zn 4s conduction) get the largest correction. Verify your prediction against the plot.

---

## Problem 5: Take-aways

### Part A

Summarise the differences between Exercise 1 (ozone, ΔSCF, KS orbitals) and Exercise 2 (ZnO, DFPT, Wannier orbitals). For each difference, explain *why* the choice is appropriate to the system.

### Part B

In both exercises we obtained excellent agreement with experiment for a charged-excitation property (IP/EA for ozone; band gap for ZnO), starting from a semi-local DFT functional that on its own gives the wrong answer by several eV. What is the physical content of the screening parameters $\alpha_i$, and why does enforcing the Koopmans condition fix the eigenvalues so dramatically?

---

[^Colonna2022]: N. Colonna, R. De Gennaro, E. Linscott, and N. Marzari, *Koopmans Spectral Functionals in Periodic Boundary Conditions*, J. Chem. Theory Comput. **18**, 5435 (2022). [doi:10.1021/acs.jctc.2c00161](https://doi.org/10.1021/acs.jctc.2c00161).
