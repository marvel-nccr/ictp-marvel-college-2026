# Infrared spectra from molecular dynamics: reading and exercise

Activate the environment for this exercise:

```bash
cdd 08
cd exercise_1
workon metatomic
```

## 1. Background and context

Infrared (IR) spectroscopy measures how a sample absorbs infrared light. Those photons
carry the energy of **molecular vibrations**, so an IR spectrum is a direct window onto
how atoms are bonded and how they move: bond stiffness, hydrogen bonding, local structure,
dynamics. It is one of the most widely used probes in chemistry, physics, materials
science and biology.

The observable is a frequency-resolved absorption coefficient, **n(ω)·α(ω)**. Each band
carries three pieces of information:

- **position**: the vibrational frequency (bond stiffness),
- **width**: the vibrational lifetime and environment,
- **intensity**: how strongly the dipole responds to the motion.

## 2. Reading material

The recipe [**"Machine-learned dipoles and infrared spectroscopy of liquid water"**](https://atomistic-cookbook.org/examples/water-ir-spectrum/water-ir-spectrum.html)
computes the IR spectrum of *liquid water* from molecular dynamics (MD). The key ideas:

- **Linear response theory and Green-Kubo formula.** The IR lineshape is proportional to the
  Fourier transform of the equilibrium **dipole–dipole autocorrelation function**,

  $n(\omega)\alpha(\omega)\ \propto\ \omega^2 \int \langle \boldsymbol{\mu}(0)\cdot \boldsymbol{\mu}(t) \rangle e^{-i \omega t} dt ,$

  where $n(\omega)$ In practice you run dynamics collecting the total dipole
  **μ(t)**, and take its power spectrum. Anharmonicity, mode coupling,
  temperature and the liquid environment come along for free with the dynamics.

- **Two ingredients from the simulation.**
  1. *Forces* to drive the MD, here from a **machine-learning interatomic potential (MLIP)**
     (PET-MAD fine-tuned to a dedicated SCAN dataset).
  2. *The dipole* **μ(t)** along the trajectory. We will compare a crude model that uses **fixed point
     charges** with a **machine-learned dipole** model and see the limits of the former.

- **Workflow.** Point-charge baseline → fine-tune a joint energy+dipole model → run an MD
  trajectory → evaluate the ML dipole along it → compute and compare the spectra against
  the point-charge baseline and against experiment. The recipe also discusses
  **equivariance** (how a dipole, being a vector, must rotate with the molecule).

**Read the recipe and make sure you understand it before starting the exercise:**

👉 https://atomistic-cookbook.org/examples/water-ir-spectrum/water-ir-spectrum.html

By the end you should have some ideas about these questions: *How do you train an ML model for dipoles? Why
does the point-charge model fail? What does the autocorrelation/power-spectrum step
actually compute? What do the band positions, widths and intensities tell you?*

## 3. Exercise: the IR spectrum of a single gas-phase water molecule

You will now compute the IR spectrum of one water molecule in the gas phase (no periodic
box). It is small enough to run on the virtual machine or on your laptop in a couple of minutes.

**Goal.** Get the vibrational IR spectrum of a single H₂O molecule from an MD trajectory,
and identify its vibrational bands.

The parts that differ from the cookbook recipe are provided below: running MD with ASE, a simplified function to compute the
power spectrum, and one to remove the molecular rotation (the *Eckart frame*). What is left for
you is building the molecule and setting up the model so it returns
both the forces (to run the dynamics) and the dipole (to build the spectrum), and finally putting the pieces together and plotting.

### Provided code

**(i) Get the model.** The fine-tuned energy+dipole checkpoint
`pet-mad-xs-v1.5.0_SCAN_dipole.ckpt` is stored elsewhere on your virtual machine
(assuming you have run `update` before starting).

Export it to a TorchScript model the calculator can load. (In a Jupyter notebook, prefix
the line with `!` to run it as a shell command.)

```bash
# export the checkpoint to a TorchScript model (.pt) that the calculator can load
mtt export /home/max/cosmo_models/pet-mad-xs-v1.5.0_SCAN_dipole.ckpt -o ./pet-mad-xs-v1.5.0_SCAN_dipole.pt
```

**(ii) Imports.** The imports for the rest of the exercise:

```python
import warnings

import ase.build
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal
from ase import units
from ase.md.bussi import Bussi
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from metatomic.torch import ModelOutput
from metatomic.torch.ase_calculator import MetatomicCalculator
from scipy.spatial.transform import Rotation

# isolated molecules have no stress -> silence two harmless warnings
warnings.filterwarnings("ignore", message=".*compute_requested_neighbors_from_options.*")
warnings.filterwarnings("ignore", message=".*invalid value encountered in scalar add.*")

MODEL = "pet-mad-xs-v1.5.0_SCAN_dipole.pt"   # the model you exported in step (i)
```

**(iii) The molecular dynamics.** A short constant-temperature (CSVR) run that
records the molecular dipole at every saved frame. It uses the two objects *you*
will define, `atoms` and `calc`.

```python
TIMESTEP = 0.5 * units.fs
TEMPERATURE = 300.0          # K
N_EQUIL = 2000               # equilibration steps (discarded)
N_FRAMES = 2000              # recorded production frames
STRIDE = 4                   # MD steps between saved frames -> save every 2 fs
DT_FS = STRIDE * TIMESTEP / units.fs
TAUT = 100 * units.fs        # thermostat coupling time


def run_single_trajectory(rng):
    """Run one constant-temperature trajectory; return (frames, dipoles)."""
    mol = atoms.copy()
    mol.calc = calc
    MaxwellBoltzmannDistribution(mol, temperature_K=TEMPERATURE, rng=rng)
    Stationary(mol)                                # remove drift, keep rotation
    dyn = Bussi(mol, TIMESTEP, temperature_K=TEMPERATURE, taut=TAUT, rng=rng)
    dyn.run(N_EQUIL)                               # equilibration

    frames = []
    dipoles = np.zeros((N_FRAMES, 3))
    for i in range(N_FRAMES):
        mol.get_potential_energy()              # run the model at the current geometry
        block = calc.additional_outputs[
            "mtt::dipole"
        ].block(0) # extra outputs are stored here
        mu = block.values.detach().cpu().numpy()[0, :, 0]
        frame = mol.copy()
        frame.info["dipole"] = mu
        frames.append(frame)
        dipoles[i] = mu
        dyn.run(STRIDE)
    return frames, dipoles
```

**(iv) The spectrum.** Turns a dipole time series into an IR-like spectrum (an FFT power
spectrum weighted by ω²).

```python
def dipole_spectrum(dipole_series, dt_fs):
    """IR-like spectrum (arbitrary units). Returns (frequencies_cm, intensity)."""
    c_cms = 2.99792458e10  # speed of light, cm/s
    freq_hz, psd = scipy.signal.periodogram(
        dipole_series, fs=1.0, detrend="constant", axis=0
    )
    psd = psd.sum(axis=1)                          # the FT of mu . mu
    freqs_cm = freq_hz / (dt_fs * 1e-15) / c_cms
    return freqs_cm, freqs_cm**2 * psd
```

**(v) The Eckart frame.** The motion of a molecule splits into three parts: overall
translation, overall rotation, and internal vibrations. Rotations of the molecule, in particular, make the overall dipole sweep around and produce a large low-frequency signal that has nothing to do with vibration. To isolate the vibrations only, we have to get rid of translational and rotational effects.

The Eckart frame is a molecule-fixed frame that strips those two rigid motions away. For
each frame the function below:

1. subtracts the center of mass, which removes the overall translation;
2. finds the rigid rotation that best superimposes the current geometry onto a fixed
   reference geometry, in a *mass-weighted* least-squares sense (the [Kabsch algorithm](https://en.wikipedia.org/wiki/Kabsch_algorithm)).
   The mass weighting is what makes rotation and vibration separate cleanly. It then
   rotates the atoms by that rotation;
3. applies the *same* rotation to the dipole, since a dipole is a vector and has to turn
   with the molecule.

What survives in the dipole is only its change due to the internal vibrations, which is
exactly what the spectrum should see. You can use the function as a black box, but it is
worth understanding what it removes and why.

```python
def apply_eckart_frame(frames, ref_frame=None, dipole_key="dipole"):
    """Remove overall translation and rotation (leaving only vibrations)."""
    if not frames:
        return []
    if ref_frame is None:
        ref_frame = frames[0]
    masses = ref_frame.get_masses()
    coords_ref = ref_frame.get_positions() - ref_frame.get_center_of_mass()
    out = []
    for frame in frames:
        new = frame.copy()
        coords = new.get_positions() - new.get_center_of_mass()
        # mass-weighted best-fit rotation onto the reference...
        R, _ = Rotation.align_vectors(coords_ref, coords, weights=masses)
        new.set_positions(R.apply(coords))
        if dipole_key in new.info:                 # ...applied to the dipole too
            new.info[dipole_key] = R.apply(new.info[dipole_key])
        out.append(new)
    return out
```

### Your task

Three short pieces are left for you.

**(a) Create the molecule**, a single non-periodic water molecule. *Hint: see the `ase.build.molecule` function.*

**(b) Set up the calculator so it returns both the energy and the dipole.** With
`MetatomicCalculator` (see how the recipe builds it), you need to:

- select the fine-tuned **SCAN energy head**, so the dynamics use the right forces
  (otherwise the calculator falls back to the base model and the molecule misbehaves), and
- request the **`mtt::dipole`** output in addition to the energy, so every step also
  returns the dipole.

Fill in the two arguments marked `...`:

```python
dipole_request = {
    "mtt::dipole": ModelOutput(
        quantity="", unit="", per_atom=False, explicit_gradients=[]
    )
}

calc = MetatomicCalculator(
    MODEL,
    additional_outputs=...,   # <- request the dipole (use dipole_request)
    variants=...,             # <- select the "scan" energy head
    device="cpu",
)
atoms.calc = calc
```

**(c) Run the dynamics, build the spectra, and plot.** Use the provided functions to get
the spectrum from the raw (lab-frame) dipole and from the Eckart-projected dipole, then
make the plot:

```python
rng = np.random.default_rng(1)
frames, dipoles = run_single_trajectory(rng)

# raw spectrum (includes the molecule's rotation/translation)
freqs, spec_lab = dipole_spectrum(dipoles, DT_FS)

# vibration-only spectrum (rotation/translation projected out)
eckart = apply_eckart_frame(frames, ref_frame=ase.build.molecule("H2O"))
dipoles_eckart = np.array([f.info["dipole"] for f in eckart])
freqs, spec_eckart = dipole_spectrum(dipoles_eckart, DT_FS)

# TODO: plot spec_lab and spec_eckart against freqs (e.g. over 0-4000 cm^-1) and compare
```

Which curve shows clean vibrational bands, and which one is dominated by the molecular
rotation? Why?

### Experimental gas-phase spectrum of H₂O
- bend ν₂ ≈ **1595 cm⁻¹**
- symmetric stretch ν₁ ≈ **3657 cm⁻¹**
- asymmetric stretch ν₃ ≈ **3756 cm⁻¹**
- Check actual experimental data at [NIST](https://webbook.nist.gov/cgi/inchi?ID=C7732185&Type=IR-SPEC&Index=0). What do you see? How does it compare with your simulated spectrum?

A non-linear triatomic has exactly three vibrations (a bend and two O–H stretches), so make
sure you can spot all three. Your classical MD peaks will sit a bit below the experimental
values, which is expected for classical dynamics (anharmonicity, no zero-point energy).
Don't try to match them exactly; aim for the right pattern.

### Optional extensions
- **Average** several trajectories (different `rng` seeds) for smoother, less noisy peaks.
- Replace the ML dipole with **fixed point charges** and compare: how do the band
  *intensities* change?
- (Advanced) Compare your dynamical peaks with the **harmonic** frequencies (normal-mode analysis at
  the relaxed geometry, e.g. with `ase.vibrations.Vibrations`). How do the harmonic lines compare with the
  peak centers from MD? Can you guess why?
