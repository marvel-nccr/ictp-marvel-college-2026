# NMR-driven structure determination with ShiftML3

Activate the environment for this exercise:

```bash
cdd 8
cd exercise_3
workon nmr
```

## 1. Background and context

**NMR-driven structure determination** combines solid-state NMR measurements with computational
crystal-structure prediction (CSP). Given a pool of candidate crystal structures, one
identifies the experimentally observed polymorph by scoring each candidate against
measured NMR data. This is especially valuable for organic pharmaceuticals: they are
often obtained only as powders (ruling out single-crystal X-ray diffraction), and
different polymorphs can have very different solubility and bioavailability.

A solid-state NMR experiment on a powdered organic crystal reports a set of *chemical
shifts* $\delta_i$ (ppm), one for each magnetically distinct nuclear site. Computationally,
what is naturally accessible is the *chemical shielding* $\sigma_i$ — the response of the
local electronic structure to the applied magnetic field. The two are related by an
approximately linear *calibration*:

$$\delta_i = a\,\sigma_i + b$$

The slope $a \approx -1$ (shielding *reduces* the resonance frequency relative to a
reference compound, hence the sign), and the intercept $b$ encodes the reference compound
shielding. Both constants depend on the level of theory.

Given a candidate crystal structure $X$, we predict shieldings $\{\sigma_i(X)\}$,
for $N$ atoms in the structure, calibrate them to shifts, 
and compute the root-mean-square error against the experimental
spectrum for each structure:

$$\mathrm{RMSE}(X) = \sqrt{\frac{1}{N}\sum_{i=1}^{N}\left(\delta_i^{\mathrm{pred}}(X)\
-\delta_i^{\mathrm{exp}}\right)^{2}}$$

The structure with the **lowest RMSE** is the best match. This works because the isotropic
chemical shielding at each nucleus is highly sensitive to its local crystal environment,
so different polymorphs produce measurably different <sup>1</sup>H spectra.

### The noise floor

Even for the *correct* structure, a finite RMSE is expected due to:
- residual DFT functional errors,
- the ShiftML3 model error,
- vibrational averaging not captured by static DFT.

### The dataset

The dataset contains 10 candidate crystal structures of a pharmaceutical compound, each
with pre-computed GIPAW (DFT) reference shieldings stored as per-atom arrays. Your goal
is to identify the correct structure using ShiftML3 predictions compared against a
measured solid-state <sup>1</sup>H spectrum.


## 2. Reading material

The recipe [**"NMR-shielding-driven structure determination with
ShiftML3"**](https://atomistic-cookbook.org/examples/shiftml-structure-match/shiftml-structure-match.html)
works through the cocaine benchmark, identifying the correct polymorph from a pool of
29 candidates using both GIPAW-DFT and ShiftML3 shieldings. The key ideas:

- **ShiftML3.** An ensemble of 8 PET models trained on GIPAW-level shieldings for
  ~14 000 molecular crystals from the Cambridge Structural Database. It delivers
  DFT-level accuracy for <sup>1</sup>H structure matching in seconds per candidate — versus
  hundreds of CPU-hours per GIPAW calculation.

- **Shielding → shift calibration.** The model outputs shieldings $\sigma_i$; the
  experiment reports shifts $\delta_i$. A linear calibration $\delta_i = a\sigma_i + b$
  is needed. The slope is conventionally fixed at $a = -1$; the intercept can be fitted
  per-structure (aligning predicted and experimental means, thereby comparing spectral
  *patterns*) or globally (a single pair from a broad benchmark, preserving absolute
  shift information).

- **Structure matching by RMSE.** Each candidate is scored against the experimental
  <sup>1</sup>H spectrum; the one with the lowest RMSE is taken as the best match. The method
  works because <sup>1</sup>H shieldings are highly sensitive to local environments — different
  polymorphs produce detectably different spectra.

- **Noise floor.** Even for the correct structure a finite irreducible RMSE is expected
  from DFT errors and the ML model's finite accuracy. Candidates whose RMSE falls
  within this band ($0.33 \pm 0.16$ ppm for PBE/GIPAW) are statistically
  indistinguishable from experiment.


👉 https://atomistic-cookbook.org/examples/shiftml-structure-match/shiftml-structure-match.html

By the end you should have some ideas about these questions: *Why are <sup>1</sup>H shieldings
sensitive enough to distinguish crystal polymorphs? What are the trade-offs between
per-structure and global calibration? What is the noise floor, and what does it mean
practically if two candidates both fall within it?*

## 3. Exercise: comparing DFT-to-experimental calibration methods

In this exercise you will apply the NMR-crystallography workflow to **AZD8329**, a
pharmaceutical compound. You have 10 candidate crystal structures with pre-computed
GIPAW shieldings. You will run ShiftML3 inference on all candidates and compare three
calibration strategies — per-structure intercept fitting, exclusion of the acidic
proton, and global calibration — examining how each affects the ability to identify the
correct polymorph.

The following are the assigned shifts (in ppm) for the molecular crystal experimentally
observed, for which we want to match the structure. You'll need these later.
```python
{
    "1":             6.92,
    "2":             8.69,
    "3":             9.01,
    "4":             8.47,
    "5":            15.37,
    "6":             7.73,
    "7":             9.64,
    "8":             2.90,
    "9":             1.78,
    "10,11":         1.88,
    "12":            1.80,
    "13":            1.60,
    "14":            0.44,
    "15":            1.54,
    "16,17":         1.88,
    "18,19":         0.80,
    "20":            1.00,
    "21,22":         1.74,
    "23,24,25,26,27,28,29,30,31": 0.73,
}
```

Later in the exercise, you'll need to identify the acidic proton in order to remove it
from the shift calibration. In normal organic molecules, this sits at 10-13 ppm, but can
be higher at 15+ ppm if the environment is particularly deshielded (i.e. in the case of
aromatic functional groups that withdraw electron density.)

### Provided code

**(i) Get the model.**

The ShiftML3 model is an ensemble of 8 members, and are
stored elsewhere on your virtual machine (assuming you have run `update` before
starting). So that the code can access them, they should be copied to the cache directory:

```python
import shutil
from platformdirs import user_cache_path

shutil.copytree("/home/max/cosmo_models/shiftml", user_cache_path() / "shiftml", dirs_exist_ok=True)
```

Initialize the ShiftML3 calculator:

```python
from shiftml.ase import ShiftML

calculator = ShiftML("ShiftML3")
```

**(ii) Imports.**

The imports for the rest of the exercise:

```python
from collections import OrderedDict

import chemiscope
import matplotlib.pyplot as plt
import numpy as np
from ase.io import read

from sklearn.metrics import root_mean_squared_error

NUM_H_PER_MOLECULE = 31
```

**(iii) Load molecular crystal candidates and visualise.**

Load with `ASE`:

```python
frames = read("azd_molecular_crystals.xyz", ":")

print(f"Loaded {len(frames)} candidate structures")
print(f"Atoms per unit cell: {len(frames[0])}")
```

Visualise with `chemiscope`:

```python
chemiscope.show(
    frames,
    mode="structure",
    settings=chemiscope.quick_settings(
        trajectory=True, structure_settings={"unitCell": True}
    ),
)
```

**(iv) Helper functions**

`assign_shieldings` picks one molecule's hydrogen shieldings from the unit cell (the
unit cell may contain several symmetry-equivalent copies) and averages over
symmetry-equivalent proton groups listed in `assigned_experimental_shifts`. For
example, a freely-rotating methyl group contributes three H atoms that all map to a
single observed NMR peak, so their predicted shieldings are averaged before comparison.

```python

def assign_shieldings(per_h_shieldings, assigned_experimental_shifts):
    """
    Pick one molecule's H shieldings from the unit cell and average over
    symmetry-equivalent groups listed in assigned_experimental_shifts.
    """
    per_mol = per_h_shieldings.reshape(NUM_H_PER_MOLECULE, -1)[:, 0]
    out = []
    for atom_string in assigned_experimental_shifts.keys():
        idx = [int(s) - 1 for s in atom_string.split(",")]
        out.append(per_mol[idx].mean())
    return np.array(out)

```

`calibrated_rmse` applies the linear calibration $\delta = a\sigma + b$ for each
candidate and returns one RMSE against experiment per candidate. If no intercept is
provided, it is fitted per-structure by matching the predicted and experimental shift
*means* — this corrects for the unknown reference-compound offset while preserving
the spectral pattern information.

```python

def calibrated_rmse(shieldings_per_candidate, experimental_shifts,
                    slope=-1.0, intercept=None):
    """RMSE for each candidate after linear calibration σ → δ."""
    rmses = []
    for sigmas in shieldings_per_candidate:
        b = intercept
        if b is None:
            b = np.mean(experimental_shifts) - slope * np.mean(sigmas)
        predicted_shifts = slope * sigmas + b
        rmses.append(root_mean_squared_error(predicted_shifts, experimental_shifts))
    return np.array(rmses)
```

`make_lollipop_plot` produces the standard visualisation for NMR crystallography:
each candidate appears as a vertical lollipop whose height is its RMSE against
experiment. GIPAW (blue) and ShiftML3 (orange) are shown side by side; the grey band
marks the noise floor — candidates inside it are statistically indistinguishable from
the correct structure.

```python

def make_lollipop_plot(frames, rmse_gipaw, rmse_sml):
    """Lollipop plot comparing GIPAW and ShiftML3 RMSE vs experiment."""
    candidate_idx = np.arange(len(frames))
    dx = 0.18
    fig, ax = plt.subplots(figsize=(8.5, 5.2), constrained_layout=True, dpi=120)

    band_lo, band_hi = 0.33 - 0.16, 0.33 + 0.16
    ax.axhspan(band_lo, band_hi, color="0.85", alpha=0.7, zorder=0,
               label=r"DFT vs experiment noise floor ($0.33 \pm 0.16$ ppm)")

    ax.vlines(candidate_idx - dx, 0, rmse_gipaw, color="C0", lw=2.5, alpha=0.85, zorder=2)
    ax.vlines(candidate_idx + dx, 0, rmse_sml,   color="C1", lw=2.5, alpha=0.85, zorder=2)
    ax.scatter(candidate_idx - dx, rmse_gipaw, color="C0", label="GIPAW (DFT reference)",
               s=60, edgecolor="white", lw=0.9, zorder=3)
    ax.scatter(candidate_idx + dx, rmse_sml,   color="C1", label="ShiftML3",
               s=60, edgecolor="white", lw=0.9, zorder=3)

    ax.set_xlabel("Candidate structure index", fontsize=13)
    ax.set_ylabel(r"$^1$H shift RMSE / ppm", fontsize=13)
    ax.set_xticks(candidate_idx[::2])
    ax.set_xlim(-0.7, len(frames) - 0.3)
    ax.set_ylim(0, max(rmse_sml.max(), rmse_gipaw.max()) * 1.15)
    ax.grid(axis="y", color="0.92", lw=0.6, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper left", frameon=False, fontsize=12)
```

**(v) Defining the experimental shifts**

Complete the dict with the assigned experimental shifts provided earlier.

```python
assigned_experimental_shifts = OrderedDict({
    "1":             6.92,
    "2":             8.69,
    ...,  # TODO: copy the rest of the experimental shieldings here
})
```

**(vi) Running ShiftML3 inference**

This may take around 30 seconds. Remember this is many orders of magnitude faster than DFT!

```python
shieldings_sml   = []
shieldings_gipaw = []

for frame in frames:
    is_h = frame.get_atomic_numbers() == 1

    sml   = calculator.get_cs_iso(frame).ravel()[is_h]
    gipaw = frame.arrays["CS"][is_h]

    shieldings_sml.append(assign_shieldings(sml,   assigned_experimental_shifts))
    shieldings_gipaw.append(assign_shieldings(gipaw, assigned_experimental_shifts))

shieldings_sml   = np.array(shieldings_sml)
shieldings_gipaw = np.array(shieldings_gipaw)
print("Inference complete.")
```

**(vii) Default calibration**

The following performs standard per-structure intercept fitting with the slope fixed at
$-1$, and computes the RMSE for each candidate. The resulting lollipop plot is generated.


```python
# Compute RMSEs
experimental_shifts = np.array(list(assigned_experimental_shifts.values()))
rmse_sml = calibrated_rmse(shieldings_sml, experimental_shifts)
rmse_gipaw = calibrated_rmse(shieldings_gipaw, experimental_shifts)

best = int(np.argmin(rmse_sml))
print(f"Best ShiftML3 candidate: #{best}  (RMSE = {rmse_sml[best]:.3f} ppm)")
print(f"Best GIPAW    candidate: #{int(np.argmin(rmse_gipaw))}  (RMSE = {rmse_gipaw.min():.3f} ppm)")

# Plot figure
make_lollipop_plot(frames, rmse_gipaw, rmse_sml)
```

**(viii) Excluding the acidic proton**


In absence of a proper global
calibration and especially without proper finite temperature sampling
it may be useful to exclude acidic protons from the RMSE computation. After all due to 
broadining effects it might even be difficult to determine the shift value experimentally
with sufficient accuracy.

The variable `acidic_proton_key` below gives the key in `assigned_experimental_shifts`
corresponding to the acidic proton. Using the information given to you above, identify
the key of this proton (a string number, like `"X"`).

```python
acidic_proton_key = ...  # TODO

slice_idxs = [
    i
    for i, key in enumerate(assigned_experimental_shifts.keys())
    if key != acidic_proton_key   # YOUR CONDITION HERE
]

print(f"Using {len(slice_idxs)} peaks (excluded key '{acidic_proton_key}')")
```

Then recompute and plot again the lollipop.

```python
# Re-compute RMSEs
rmse_sml = calibrated_rmse(shieldings_sml[:,   slice_idxs], experimental_shifts[slice_idxs])
rmse_gipaw = calibrated_rmse(shieldings_gipaw[:, slice_idxs], experimental_shifts[slice_idxs])

# Plot figure
make_lollipop_plot(frames, rmse_gipaw, rmse_sml)
```

**(ix) Using an experimentally-calibrated global calibration**

The per-structure intercept fit above adjusts $b$ separately for each candidate to
align the predicted mean shift with the experimental mean. This is convenient but
discards information carried by the absolute offset, which can differ between
polymorphs due to differences in crystal packing.

A better approach is a **global calibration** — a single $(a, b)$ pair fitted once by
linear regression of computed shieldings against experimental shifts across many
molecular crystals. This corrects systematic DFT functional errors (slope $a \neq -1$)
and removes the per-structure degree of freedom. The values below were obtained by
fitting PBE/GIPAW shieldings against experimental <sup>1</sup>H shifts for a diverse benchmark
set:

$$a = -0.9024 \qquad b = 28.05 \text{ ppm}$$

With a global calibration the acidic proton can also be *included* again, given
that in practice the slope corrects strongly for finite temperature effects in
the computed shielding values of hydrogen-bonded and acidic protons.

```python
# Define slope and intercept from given values
slope = ...  # TODO
intercept = ...  # TODO

# Compute RMSEs
rmse_sml = calibrated_rmse(
    shieldings_sml, experimental_shifts, slope=slope, intercept=intercept
)
rmse_gipaw = calibrated_rmse(
    shieldings_gipaw, experimental_shifts, slope=slope, intercept=intercept
)

# Plot figure
make_lollipop_plot(frames, rmse_gipaw, rmse_sml)
```


## 4. Discussion points

Once you've finished the exercise, or while you go, think about these questions and discuss with your neighbor.

<details>
<summary> Why does the grey band in the plot represent?</summary>

The grey band is the **noise floor**: the irreducible RMSE expected even for the
correct crystal structure, arising from DFT functional errors, finite ShiftML3 model
error, and the neglect of vibrational averaging in the static calculations. Candidates
whose RMSE falls inside this band (here $0.33 \pm 0.16$ ppm, estimated from a broad
benchmark of organic solids) are statistically indistinguishable from experiment —
the measured spectrum alone cannot rule them out.

</details>


<details>
<summary> What would it mean if more than one candidate was within the noise floor?</summary>

In an ambiguous case, the <sup>1</sup>H isotropic RMSE alone cannot resolve the degeneracy — the
lollipop plot can confidently *exclude* high-RMSE candidates but cannot choose between
those inside the band. In practice one could then: include shieldings from additional
nuclei (e.g. <sup>13</sup>C, <sup>15</sup>N); incorporate the full chemical shielding tensor
(anisotropy); or collect additional experimental data. The key figure of merit for
structure determination is an unambiguous unique minimum, not merely a low absolute
RMSE.

</details>

<details>
<summary> What changed when we excluded the acidic proton and why?</summary>

After excluding the acidic proton the calibration fits the bulk of the spectrum much
more accurately, pulling the RMSE values down toward (or into) the noise floor.
The correct candidate should now stand out clearly as the lowest-RMSE structure —
and ideally both GIPAW and ShiftML3 agree on the same winner.

Notice also that the overall *ranking* of candidates often improves even if the absolute
RMSE values are still slightly above the floor: the key figure of merit for structure
determination is whether the correct polymorph is uniquely identifiable, not whether it
sits precisely inside the grey band.

</details>


<details>
<summary> Why does the global calibration help?</summary>

The per-structure intercept fit assumes $a = -1$ exactly. In practice, PBE/GIPAW
systematically underestimates electron density in certain chemical environments,
giving a slope slightly different from $-1$. Furthermore, the slopes have been determined
on static 0 K structures, neglecting quantum-nuclear effects (delocalizations of the protons)
and generally finite temperature effects.
The global calibration corrects this with a slope of $-0.9024$, 
which better captures the true shielding-to-shift relationship
across diverse chemical environments. 

</details>


<details>
<summary> What methods could be used to model better the acidic protons?</summary>

Acidic (exchangeable) protons are particularly sensitive to nuclear quantum effects:
tunnelling and zero-point delocalization shift the proton away from its classical
equilibrium position, systematically moving the predicted resonance. Better treatments
include:

- **Path-integral MD.** Sampling the quantum nuclear distribution via ring-polymer
  simulations captures vibrational averaging and proton delocalization, bringing
  predicted shifts closer to experiment.
- **Empirical corrections.** Fitting a separate calibration for exchangeable protons
  using a benchmark of known acidic-proton shifts can reduce systematic errors at low
  cost.
- **Electronic structure improvements.** GGA DFT functionals tend to overstabilize 
  the localization of protons, why in reality fractional protonation states (smeared out
  proton positions) are observed. Hybrid DFT might be required to accurately describe
  bond breaking events and hence structural ensembles generated at this level of theory
  and the average shieldings predicted from it, might improve overall prediction accuracy.

</details>
