# Import the sscha code
import cellconstructor as CC, cellconstructor.Phonons
import cellconstructor.Structure, cellconstructor.calculators

# Import the force field of Gold
import ase, ase.calculators
from ase.calculators.emt import EMT

# Import numerical and general pourpouse libraries
import sys, os

gold_structure = CC.Structure.Structure()
gold_structure.read_generic_file("Au.cif")

# Get the force field for gold
calculator = EMT()

gold_harmonic_dyn = CC.Phonons.compute_phonons_finite_displacements(gold_structure, calculator, supercell = (4,4,4))

# Impose the symmetries and
# save the dynamical matrix in the quantum espresso format
gold_harmonic_dyn.Symmetrize()
gold_harmonic_dyn.save_qe("gold_harmonic_dyn")
