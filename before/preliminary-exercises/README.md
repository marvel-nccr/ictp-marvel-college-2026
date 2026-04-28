# Preliminary exercises

This exercise forms the very first part of the introduction to DFT hands-on session. **Please complete it before you arrive in Trieste.** It introduces running `Quantum ESPRESSO` calculations on bulk sodium chloride (NaCl). The remaining problems will follow during the in-person session (see [day-01](../../day-01/)).

> **Goal**
>
> The exercise is a task to complete (along with some related questions and notes). **Make sure you read the notes before attempting the exercise!**

## Files provided

- `NaCl_primitive.scf.in` — input file for the NaCl primitive cell
- `pseudopotentials/` — pseudopotentials for Na and Cl
- `script.sh` — example bash script for running a sweep over a parameter

To run a single calculation:

```bash
pw.x < NaCl_primitive.scf.in > NaCl_primitive.scf.out
```

---

## Problem 1: Running your first calculation

### Part A

Calculate the total energy of NaCl using the provided input file. Plot how the total energy changes during the self-consistent cycle.

> **Note**
>
> Command-line tools such as `grep` can help you quickly extract values from output files.

### Part B

What criterion does the code use to decide when the self-consistent cycle has converged?

> **Note**
>
> The `estimated scf accuracy` printed in `Quantum ESPRESSO` output files is defined in Appendix A.1 of Giannozzi *et al.*, *J. Phys. Condens. Matter* **21**, 395502 (2009).

### Part C

Does the energy decrease monotonically during the self-consistent cycle? Is this what we expect?
