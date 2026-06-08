# NMR-driven structure determination with ShiftML3

Activate the environment for this exercise:

```bash
workon nmr
```

## 1. Background and context

<TODO>: clean up this section.

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
calibrate them to shifts, and compute the root-mean-square error against the experimental
spectrum:

$$\mathrm{RMSE}(X) = \sqrt{\frac{1}{N}\sum_{i=1}^{N}\!\left(\delta_i^{\mathrm{pred}}(X)
- \delta_i^{\mathrm{exp}}\right)^{2}}$$

The structure with the **lowest RMSE** is the best match. This works because the isotropic
chemical shielding at each nucleus is highly sensitive to its local crystal environment,
so different polymorphs produce measurably different $^1$H spectra.

### The noise floor

Even for the *correct* structure, a finite RMSE is expected due to:
- residual DFT functional errors,
- the ShiftML3 model error,
- vibrational averaging not captured by static DFT.

### The dataset

The dataset contains 10 candidate crystal structures of a pharmaceutical compound, each
with pre-computed GIPAW (DFT) reference shieldings stored as per-atom arrays. Your goal
is to identify the correct structure using ShiftML3 predictions compared against a
measured solid-state $^1$H spectrum.


## 2. Reading material

The recipe [**"NMR-shielding-driven structure determination with
ShiftML3"**](https://atomistic-cookbook.org/examples/shiftml-structure-match/shiftml-structure-match.html)
<TODO>. The key ideas:

<TODO>


👉 https://atomistic-cookbook.org/examples/shiftml-structure-match/shiftml-structure-match.html

By the end you should have some ideas about these questions: <TODO>

## 3. Exercise: comparing DFT-to-experimental calibration methods

<TODO>

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

Later in the exercise, you'll need to identify the acidic proton in oredr to remove it
from the shift calibration. In normal organic molecules, this sits at 10-13 ppm, but can
be higher at 15+ ppm if the environment is particularly deshielded (i.e. in the case of
aromatic functional groups that withdraw electron density.)

<TODO>: clean up the following section

Eventually, you use a global calibration scheme with the following parameters, computed
against experiment.

The per-structure intercept fitting above adjusts $b$ separately for each candidate to
align the mean predicted shift with the experimental mean. This is convenient but has
a drawback: it *removes any information carried by the absolute offset*, which can vary
between polymorphs due to differences in crystal packing.

A better approach is to use a **global calibration** — a single $(a, b)$ pair fitted by
linear regression of computed shieldings against experimental shifts across many molecular
crystals. This accounts for systematic DFT functional errors (reflected in the slope
$a \neq -1$) and reference-compound offset without discarding polymorph-specific
information.

The values below were obtained by fitting PBE/GIPAW shieldings against experimental
$^1$H shifts for a diverse benchmark set:

$$a = -0.9024 \qquad b = 28.05 \text{ ppm}$$

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

<TODO>: explain briefly the function

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

<TODO>: explain briefly the function

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

<TODO>: explain briefly the function

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
    ...,  # TODO
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


As acidic protons are well captured at this DFT level, in absence of a proper global
calibration it may be useful to exclude acidic protons too.

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

**(ix) Using an experimentally-calculated global calibration**

<TODO>: explain

```python
# Define slope and intercept from given values
slope     = ...
intercept =  ...

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

Once you've finished the exercise, or while you go, think about these prompts and discuss with your neighbor.

<details>
<summary>💭 Why does the grey band in the plot represent?</summary>

<TODO>: briefly explain.

</details>


<details>
<summary>💭 What would it mean if more than one candidate was within the noise floor?</summary>

<TODO> In an ambiguous case, the lollipop plot only really helps *exclude* structures
with high RMSE relatve to experiment. Further screening would be needed to distinguish
between equally viable candidates.

</details>

<details>
<summary>💭 What changed when we excluded the acidic proton and why?</summary>

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
<summary>💭 Why does the global calibration help?</summary>

The per-structure intercept fit assumes $a = -1$ exactly. In practice, PBE/GIPAW
systematically underestimates electron density in certain chemical environments,
giving a slope slightly different from $-1$. The global calibration corrects this
with a slope of $-0.9024$, which better captures the true shielding-to-shift relationship
across diverse chemical environments. The fixed intercept also removes the freedom
that was inadvertently masking polymorph-specific offsets.

In addition, with the global calibration the acidic proton peak can be *included* again
(it is no longer a special case because the calibration was not fitted to these
particular structures): the result uses all available experimental information.

</details>


<details>
<summary>💭 What methods could be used to model better the acidic protons?</summary>

<TODO>

</details>



## 5. Further reading

The following papers, in chronological order, are good reads for density learning:

1. <TODO>