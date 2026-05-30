# Import the sscha code
import sscha, sscha.Ensemble, sscha.SchaMinimizer
import sscha.Relax, sscha.Utilities

# Import the cellconstructor library to manage phonons
import cellconstructor as CC, cellconstructor.Phonons
import cellconstructor.Structure, cellconstructor.calculators

# Import the force field of Gold
import ase, ase.calculators
from ase.calculators.emt import EMT

# Import numerical and general pourpouse libraries
import numpy as np, matplotlib.pyplot as plt
import sys, os

TEMPERATURE = 300
N_CONFIGS = 50
MAX_ITERATIONS = 20

# Get the force field for gold
calculator = EMT()

# Initialize the random ionic ensemble
nqirr = 13 # How many irreducible q points?
gold_harmonic_dyn = CC.Phonons.Phonons("gold_harmonic_dyn", nqirr)
ensemble = sscha.Ensemble.Ensemble(gold_harmonic_dyn, TEMPERATURE)

# Initialize the free energy minimizer
minim = sscha.SchaMinimizer.SSCHA_Minimizer(ensemble)
minim.set_minimization_step(0.01)

# Initialize the NVT simulation
relax = sscha.Relax.SSCHA(minim, calculator, N_configs = N_CONFIGS,
max_pop = MAX_ITERATIONS)

# Define the I/O operations
# To save info about the free energy minimization after each step
ioinfo = sscha.Utilities.IOInfo()
ioinfo.SetupSaving("minim_info")
relax.setup_custom_functions(custom_function_post = ioinfo.CFP_SaveAll)

# Run the NVT simulation
relax.relax(get_stress = True)

# Save the final dynamical matrix
relax.minim.dyn.save_qe("sscha_T300_dyn")
