# MD with MLIPs

**Tutors:** Paolo Pegolo, Joseph W. Abbott, Matthias Kellner, Filippo Bigi

---

## 1. Background and context

A machine-learning interatomic potential (MLIP) predicts energies and forces
(and stress), and you can use it to run molecular dynamics in place of expensive
*ab initio* calculations. There are (at least) two ways to get one. The classic
way is to train it **from scratch** on a dataset of reference structures:
accurate, but it might need a lot of data.

The other way starts from a *universal* MLIP such as **PET-MAD**, trained across
a broad chemical space. For most systems a universal model already works well
out of the box, with no extra training. When you need more, e.g. better accuracy
on a specific system you care about, or a different (more accurate, more
expensive) level of theory you only have a little data for, you can
**fine-tune** it: start from the pretrained checkpoint and keep training on your
small dataset. Because the model already knows a lot of chemistry, a few hundred
structures can be enough, whereas training from scratch might need orders of
magnitude more.

In this exercise you will do both on the same tiny ethanol dataset: train a PET
model from scratch, fine-tune PET-MAD-XS, compare how well the two reproduce the
reference data, and finally run molecular dynamics of a single ethanol molecule
with the fine-tuned model.

## 2. Reading material

The recipe [**"Introduction to foundational models for molecular
dynamics"**](https://atomistic-cookbook.org/examples/water-md/water-md.html)
takes a *single* universal model, PET-MAD-XS, and uses it unchanged to run
molecular dynamics on four qualitatively different aqueous systems (liquid
water, a salt solution, an ethanol-water mixture, superionic ice). Read it to
get a feel for what a foundational MLIP can do straight out of the box, and for
how an MD simulation with one is set up and run.

The exercise then goes a step further and trains models of your own. The
fine-tuning part has three things worth keeping in mind:

- **Start from a checkpoint, not from scratch.** Fine-tuning continues training
  from the pretrained model checkpoint (a `.ckpt` file). You point `metatrain`
  at it through the `finetune.read_from` field of the training options.

- **A new energy head for your reference.** In general, the new reference
  energies will not have the same absolute scale as the data the universal model
  was trained on (different code, functional, pseudopotentials, even a different
  zero of energy). This means that the new model will need a new *head*,
  including a different compositional baseline. Instead of overwriting the
  original `energy` head, you train a *new* one, a target named, e.g.,
  `energy/finetune`. Note that any other "variant" name could be chosen instead
  of `finetune`.

- **The `metatrain` workflow.** Write an `options.yaml` → `mtt train
  options.yaml -o model.pt` to train → `mtt eval model.pt eval.yaml` to evaluate
  the performance, then use the exported model from ASE or LAMMPS. At inference
  you pick the fine-tuned head by passing `variants={"energy": "finetune"}` to
  the calculator.

**Read the recipe before starting the exercise:**

👉 https://atomistic-cookbook.org/examples/water-md/water-md.html

## 3. Exercise: train, fine-tune, compare, then run MD

**Goal.** Build two ethanol models, one trained from scratch, one fine-tuned
from PET-MAD-XS, and compare how well each reproduces the reference energies and
forces, then run a short, stable MD of a single ethanol molecule with the
fine-tuned model. To solve the exercise, we suggest using a new Jupyter notebook
to analyze the training and evaluation outputs and plot what's needed.

### Provided data

The dataset, `ethanol_ev.xyz`, contains 100 ethanol configurations (9 atoms
each), each labelled with a reference `energy` (in eV) and `forces` (in eV/Å).
It is deliberately small and narrow (one molecule, just a few geometries), so it
is quick to train on.

Take a look at the second line of the file: the `Properties=...` and
`energy=...` fields tell you under which keys the energy and forces are stored.
You will need those keys in the options files below.

### (1) Train a model from scratch

First, train a PET model on the ethanol data from scratch, where the model
weights are randomly initialized. Write `options.yaml` and fill the two blanks
marked `# <-- ...`. If unsure, feel free to explore the [`metatrain`
docs](https://docs.metatensor.org/metatrain/latest/index.html) or ask questions
to the tutors.

Write the following text to `options.yaml`:

```yaml
device: cpu
seed: 0

architecture:
  name: pet
  training:
    num_epochs: 200           # from scratch the model needs many passes over the data
    learning_rate: 1e-3       # a sensible starting point for training from scratch
    batch_size: 16
    loss:
      energy:
        type: mse
        weight: 10            # weigh the energy 10x the forces in the loss
        gradients:
          positions:          # the "positions gradient" of the energy is the forces
            type: mse
            weight: 1

training_set:
  systems:
    read_from: ethanol_ev.xyz
    length_unit: angstrom
  targets:
    energy:                   # the standard energy head
      key: ...                # <-- the energy key in the xyz file
      unit: eV
      forces:
        key: ...              # <-- the forces key in the xyz file

validation_set: 0.1           # hold out 10% for validation, 40% for testing
test_set: 0.4
```

You can run the training with:

```bash
mtt train options.yaml -o model-scratch.pt
```

This trains the model and writes the exported potential `model-scratch.pt`. On a
CPU the run takes a while (at least a few minutes for a few hundred epochs).
Since this is probably your first `metatrain` run, it is worth knowing what it
prints and what it leaves behind.

**The log.** `metatrain` prints one line per epoch, of the form:

```text
Epoch: 99 | learning rate: ... | training loss: ... | training energy RMSE (per atom): ... meV | training forces RMSE: ... meV/A | validation energy RMSE (per atom): ... meV | validation forces RMSE: ... meV/A
```

Each line gives the current learning rate, the loss, and the energy and force
errors (RMSE and MAE) on the training and validation sets. Watch the
*validation* errors in particular: they should come down and stay close to the
training errors. A widening gap between training and validation means the model
is overfitting; if none of the errors move, training is too slow (often the
learning rate is too low).

**The `outputs/` folder.** Each run also creates a timestamped directory
`outputs/<date>/<time>/` containing:

- `train.log`: the full log, in case it scrolled past;
- `train.csv`: the same per-epoch metrics as a table, handy for plotting
  learning curves;
- `model_0.ckpt`, `model_100.ckpt`, …: checkpoints saved along the way, plus a
  final `model-scratch.ckpt` (a checkpoint is what you restart or fine-tune
  from);
- `options_restart.yaml`: the exact options used, with every default filled in;
- `indices/`: which structures went into the train / validation / test splits.

The `model-scratch.pt` in your working directory is the file you load for
evaluation and the MD below.

### (2) Evaluate it and make parity plots

`mtt eval` runs a model over a dataset and reports its errors. Write a small
`eval.yaml` that points at the same data and asks for the `energy` head:

```yaml
systems:
  read_from: ethanol_ev.xyz
  length_unit: angstrom
targets:
  energy:
    key: energy
    unit: eV
    forces:
      key: forces
```

and run the evaluation with:

```bash
mtt eval model-scratch.pt eval.yaml -o pred-scratch.xyz
```

You should see the energy and force RMSE/MAE printed in the terminal. At the end
of the evaluation, you should find the `pred-scratch.xyz` file with the
predictions.

A **parity plot** puts the predicted value against the reference value, point by
point. It is the quickest way to *see* whether a model is any good: the closer
the cloud sits to the diagonal, the better.

One subtlety: `metatrain` split the data into training, validation and test
sets, and the model has already seen the first two. To judge it fairly we look
only at the **test** structures, the ones it never trained on. `metatrain`
recorded that split under the run's `outputs/` folder (`indices/test.txt`), so
we load those indices and plot only those configurations.

The helpers below read the reference from `ethanol_ev.xyz` and the predictions
from the `mtt eval` output, keep only the test configurations, and save the
figure to a file.

```python
from glob import glob

import ase.io
import matplotlib.pyplot as plt
import numpy as np

frames = ase.io.read("ethanol_ev.xyz", ":")


def reference(frames):
    """Per-atom reference energies (eV) and flattened reference forces (eV/Å)."""
    e = np.array([a.get_potential_energy() / len(a) for a in frames])
    f = np.concatenate([a.get_forces().ravel() for a in frames])
    return e, f


def read_predictions(pred_file):
    """Per-atom predicted energies (eV) and flattened forces (eV/Å) from an mtt eval file.

    `mtt eval` names its output columns after the model's target: the plain `energy` head
    lands on a calculator (`get_potential_energy()` / `get_forces()`), while the
    `energy/finetune` head lands in `info` / `arrays` under keys carrying the head name.
    This reads whichever layout is present.
    """
    frames = ase.io.read(pred_file, ":")
    e, f = [], []
    for a in frames:
        energy = next((v for k, v in a.info.items() if k.startswith("energy")), None)
        if energy is None:
            energy = a.get_potential_energy()
        fkey = next((k for k in a.arrays if "force" in k.lower()), None)
        forces = a.arrays[fkey] if fkey is not None else a.get_forces()
        e.append(energy / len(a))
        f.append(forces.ravel())
    return np.array(e), np.concatenate(f)


def parity(e_ref, f_ref, predictions, title, filename, indices=None):
    """Energy and force parity plots; `predictions` is a list of (label, e_pred, f_pred).

    With `indices` (here the test-set structures) only those configurations are plotted. The
    forces are regrouped per structure first, so the selection picks whole configurations.
    """
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    if indices is not None:
        n = len(e_ref)
        e_ref = e_ref[indices]
        f_ref = f_ref.reshape(n, -1)[indices].ravel()
        predictions = [
            (label, e_pred[indices], f_pred.reshape(n, -1)[indices].ravel())
            for label, e_pred, f_pred in predictions
        ]
    for label, e_pred, f_pred in predictions:
        ax[0].scatter(e_ref, e_pred, s=10, alpha=0.5, label=label)
        ax[1].scatter(f_ref, f_pred, s=10, alpha=0.5, label=label)
    for axis, x, name in [
        (ax[0], e_ref, "energy / eV/atom"),
        (ax[1], f_ref, "force / eV/Å"),
    ]:
        axis.axline((x.min(), x.min()), slope=1, ls="--", c="k", lw=1)
        axis.set_xlabel(f"reference {name}")
        axis.set_ylabel(f"predicted {name}")
        axis.legend()
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(filename, dpi=250)
    plt.show()


# which structures metatrain held out for testing (latest run)
test_idx = np.loadtxt(sorted(glob("outputs/*/*/indices/test.txt"))[-1], dtype=int)

e_ref, f_ref = reference(frames)
e_scratch, f_scratch = read_predictions("pred-scratch.xyz")
parity(
    e_ref,
    f_ref,
    [("from scratch", e_scratch, f_scratch)],
    "from scratch",
    "parity-scratch.png",
    indices=test_idx,
)
```

### (3) Fine-tune PET-MAD-XS

Now build the second model by fine-tuning the universal PET-MAD-XS instead of
starting from random weights. The pre-trained checkpoints are distributed on
HuggingFace, but to avoid too many competing requests to the server you will
find them already downloaded in the Virtual Machine. Copy the PET-MAD-XS
checkpoint in the current exercise's folder with

```bash
cdd 07
cd exercise_4
cp /home/max/cosmo_models/pet-mad-xs-v1.5.0.ckpt .
```

Fine-tuning needs the *checkpoint* rather than an exported `.pt` model, because
the checkpoint still carries everything needed to keep training (the full
architecture and possibly the training state), while the `.pt` is a frozen model
meant only for running simulations.

Write `options-ft.yaml` and fill the three blanks. They are the fields that
connect the run to *your* checkpoint and your data.

```yaml
device: cpu
seed: 0

architecture:
  name: pet
  training:
    num_epochs: 50           # far fewer than from scratch: you start from a trained model
    learning_rate: 1e-4       # lower too, since you are adjusting a model, not building one
    batch_size: 16
    finetune:
      method: full            # update all of the model's weights
      read_from: ...          # <-- the base checkpoint you downloaded
      inherit_heads:
        # {new head: source head}: initialise the new `energy/finetune` head from
        # the pretrained `energy` head instead of from scratch
        energy/finetune: energy
    loss:
      energy/finetune:
        type: mse
        weight: 10            # weigh the energy 10x the forces in the loss
        gradients:
          positions:          # the "positions gradient" of the energy is the forces
            type: mse
            weight: 1

training_set:
  systems:
    read_from: ethanol_ev.xyz
    length_unit: angstrom
  targets:
    energy/finetune:          # a NEW head, separate from the pretrained "energy"
      key: ...                # <-- the energy key in the xyz file
      unit: eV
      forces:
        key: ...              # <-- the forces key in the xyz file

validation_set: 0.1
test_set: 0.4
```

```bash
mtt train options-ft.yaml -o model-ft.pt
```

Note that 50 epochs of fine-tuning take a fraction of the 200 or more you needed
from scratch. That is the whole point: the model already knows most of the
chemistry, so it only has to adjust.

### (4) Evaluate the fine-tuned model and compare

Evaluate it the same way. The only change in `eval-ft.yaml` is the target name:
the fine-tuned model predicts through the `energy/finetune` head, while the
reference values still live under the `energy` key in the file.

```yaml
systems:
  read_from: ethanol_ev.xyz
  length_unit: angstrom
targets:
  energy/finetune:
    key: energy
    unit: eV
    forces:
      key: forces
```

```bash
mtt eval model-ft.pt eval-ft.yaml -o pred-ft.xyz
```

Read its predictions from `pred-ft.xyz`, make its parity plot, and overlay both
models on a single figure:

```python
e_ft, f_ft = read_predictions("pred-ft.xyz")
parity(
    e_ref,
    f_ref,
    [("fine-tuned", e_ft, f_ft)],
    "fine-tuned",
    "parity-ft.png",
    indices=test_idx,
)

# overlay both models on one figure
parity(
    e_ref,
    f_ref,
    [("from scratch", e_scratch, f_scratch), ("fine-tuned", e_ft, f_ft)],
    "from scratch vs fine-tuned",
    "parity-compare.png",
    indices=test_idx,
)
```

Compare the two. Which model reproduces the reference energies and forces
better, and by how much? Remember that the fine-tuned one got there in a quarter
of the epochs, by reusing what PET-MAD already knew.

### (5) Run MD of a single ethanol molecule

Finally, use the fine-tuned model to run a short constant-temperature trajectory
of one ethanol molecule, *selecting your fine-tuned head*.

You can use whichever engine you like: a metatomic model plugs into ASE, i-PI
and LAMMPS through the same interface. The water-md recipe drives all of its MD
with **LAMMPS** (the `metatomic` pair style reading an exported `.pt` model), so
if you prefer LAMMPS you can follow the input files shown there, just pointing
the pair style at your own `model-ft.pt`. Below we stay in Python and use
**ASE** instead, which is why that code is provided for you.

Whichever engine you choose, the piece to get right is selecting the head you
just trained instead of the pretrained one (in the ASE code below that is the
`variants` argument). In LAMMPS, it would be something like `pair_style metatomic
... variant finetune`.

```python
import ase.build
import numpy as np
from ase import units
from ase.md.bussi import Bussi
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.io.trajectory import Trajectory
from metatomic.torch.ase_calculator import MetatomicCalculator

mol = ase.build.molecule("CH3CH2OH")      # a single ethanol molecule, gas phase

calc = MetatomicCalculator(
    "model-ft.pt",
    variants={"energy": "finetune"},      # <-- use the head you just trained
    device="cpu",
)
mol.calc = calc

rng = np.random.default_rng(0)
MaxwellBoltzmannDistribution(mol, temperature_K=300.0, rng=rng)
Stationary(mol)
dyn = Bussi(mol, 0.5 * units.fs, temperature_K=300.0, taut=100 * units.fs, rng=rng)
traj = Trajectory("ethanol.traj", "w", mol)
dyn.attach(traj.write, interval=1)

def report():
    # ase.build.molecule("CH3CH2OH") orders atoms as C C O H H H H H H,
    # so atom 2 is the oxygen and atom 3 the hydroxyl hydrogen.
    print(
        f"E = {mol.get_potential_energy():9.3f} eV"
        f"   T = {mol.get_temperature():6.1f} K"
        f"   O-H = {mol.get_distance(2, 3):.3f} A"
    )


dyn.attach(report, interval=200)
dyn.run(2000)

# Visualize with chemiscope
import chemiscope

traj = Trajectory("ethanol.traj")
e = [frame.get_potential_energy() for frame in traj]
T = [frame.get_temperature() for frame in traj]

chemiscope.show(
    traj,
    properties={
        "time (fs)": 0.5 * np.arange(len(e)),
        "energy (eV)": e,
        "temperature (K)": T,
    },
)
```

### What to check

- The temperature fluctuates around 300 K and does not run away.
- The molecule stays bound: the O–H distance stays near 1 Å, and the C–C, C–O
  and C–H bonds keep their equilibrium lengths. Nothing flies apart.

If the dynamics blow up, the usual culprits are: the wrong head (forgetting
`variants`, so the calculator falls back to the pretrained `energy` head instead
of yours: this can be catastrophically bad, as PET's internal features change
during fine-tuning but the original head does not), a timestep that is too large,
or a model that has not trained enough.

### Optional extensions

- **Run MD with the from-scratch model too.** Point the calculator at
  `model-scratch.pt` (no `variants` needed, it has a single `energy` head) and
  run the same dynamics. Is it as stable as the fine-tuned one?
