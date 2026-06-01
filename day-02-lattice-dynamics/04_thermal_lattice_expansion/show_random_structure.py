import cellconstructor as CC, cellconstructor.Phonons
import sscha, sscha.Ensemble
import ase, ase.visualize

# Load the quantum espresso harmonic dynamical matrix
nqirr = 13 # How many irreducible q points
harmonic_dyn = CC.Phonons.Phonons("gold_harmonic_dyn", nqirr)

# Generate a random ensemble at 600 K
temperature = 600 # K
n_configs = 100
ensemble = sscha.Ensemble.Ensemble(harmonic_dyn, temperature)
ensemble.generate(n_configs)

# Use ASE to display the first structure
ase.visualize.view(ensemble.structures[0].get_ase_atoms())
