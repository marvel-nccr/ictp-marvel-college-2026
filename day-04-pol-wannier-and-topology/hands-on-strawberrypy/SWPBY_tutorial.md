<div align="right">

**Joint ICTP MARVEL College: Materials simulations in the age of AI**  
Thursday, 4 June 2026

&nbsp;


</div>

<div align="center">

# Disorder in topological insulators with StraWBerryPy:  Hands-on session

</div>



**Tutors**: N. Baù, A. Marrazzo, R. Dal Molin, M. Dowlatabadi, R. D'Souza

**Time budget**: 30 minutes

This session continues from the [preliminary exercises](../hands-on-pythtb), where you defined the Haldane and Kane-Mele model and studied their topological phase diagrams. Topological properties are known to be robust against weak perturbation $-$ but how weak is "weak"? How does the topology look like in real space? In this tutorial, you will investigate the Haldane model with on-site (Anderson) disorder and vacancies to understand how these forms of disorder affect both the global and real-space topological properties of a material.

> **Documentation**
>
> If you need or are curious to have a look at the documentation of the code, you can find it at [this link](https://strawberrypy.readthedocs.io/en/latest/)! Consulting the documentation may also help you better understand the functions and flags used throughout the tutorial.

> **Note**
>
> Since you will run the same calculations multiple times with different parameters, using a Jupyter notebook is recommended.

---

## Problem 1: Real-space topology of clean systems

We begin with the same pristine Haldane model you studied in the previous tutorial. You may either reuse the model defined earlier or create it again from scratch. The next step is pass the model to StraWBerryPy, specifying that it should be interpreted as spinless:

```python
import numpy as np
from strawberrypy import *

haldane_unit_cell = example_models.haldane_pythtb(delta, t, t2, phi)
haldane = Model(model = haldane_unit_cell, spinful = False)
```

### Part A

Let us first examine how topology behaves in real space. Before running any code, suppose you had a probe capable of determining whether a given point in space belongs to a topological or a trivial phase. Would you expect any difference between performing this measurement in a finite sample and in a periodic system? If so, what is this due to?


### Part B

Let's now verify this numerically. Choose the set of parameters (`delta`, `t`, `t2`, `phi`) you prefer to define the Haldane model and create a supercell out of it. For example:
* for a topological phase you can use (`delta=0.1`, `t=1`, `t2=0.3`, `phi=np.pi/2`)
* for a trivial phase you can use (`delta=1.1`, `t=1`, `t2=0.3`, `phi=np.pi/8`)

Then, a supercell can be created with
```python
haldane_supercell = haldane.make_supercell(Lx = L, Ly = L)
```
where `L` is the number of unit cells repeated along the specified direction. Once the supercell is constructed, we can evaluate the single-point Chern number and compare it with the global topological invariant studied in the previous tutorial. The single-point Chern number is defined by
$$C_{sp}=-\frac{1}{\pi}\sum_n\langle \tilde u_{n\mathbf b_1}|\tilde u_{n\mathbf b_2}\rangle$$
where $|\tilde u_{n\mathbf b_j}\rangle$ are the eigenstates of the Hamiltonian at $\Gamma$ parallel transported at the Brillouin zone edge. In StraWBerryPy, this can be computed by calling
```python
haldane_supercell.single_point_chern()
```
Because the single-point invariant relies on a discretized derivative, it is affected by finite-size effects. How large must the supercell be for the result to converge? Does the single-point approach the expected value for your chosen parameters?

> **Hint:** The single-point actually converges really fast with the size of the system `L` away from a topological phase transition. Even when close to such boundary though, convergence is faster if `L` is a multiple of 3. Why?

> **Note:** The function to compute the single-point Chern number accept also the argument ``return_ham_gap`` to check the gap of the Hamiltonian at the $\Gamma$-point. Near the topological phase transition, where the gap closes, numerical instabilities arise and larger system sizes are required to obtain reliable results.

### Part C
Let's now look at the topology in real space. The local Chern marker in periodic boundary conditions is defined by
$$\mathcal{C}(\mathbf r)=-\frac{1}{2\pi}\mathrm{Im}\braket{\mathbf{r}|\big[ \mathcal P_{\mathbf b_1}, \mathcal P_{\mathbf b_2} \big] \mathcal P_{\Gamma} |\mathbf{r}}$$
where $\mathcal P_{\mathbf b_j}$ are the parallel transported ground-state projectors from the $\Gamma$-point to the Brillouin zone edges. This can be computed and plotted in StraWBerryPy by using
```python
pbc_marker = haldane_supercell.pbc_local_chern_marker()
fig, ax = postprocessing.plot_marker(
        model = haldane_supercell,
        marker = pbc_marker,
        savefig = True,
        filename = 'pbc_marker.pdf',

        # For better visualization
        vmin = -1, vmax = 1,
    )
```
The local Chern marker for a bounded system is defined by
$$ \mathfrak C(\mathbf r)=4\pi\text{Im} \braket{\mathbf{r}|\mathcal{P}\left[x,\mathcal{P}\right]\left[y,\mathcal{P}\right]|\mathbf{r}} $$
where $\mathcal P$ is the ground-state projector. To compute this quantity, we first construct a finite sample by "carving out" a bounded region from the bulk system. The finite system and its local Chern marker can be evaluated with
```python
haldane_finite = haldane.make_finite(Lx = L, Ly = L)
obc_marker = haldane_finite.local_chern_marker()
fig, ax = postprocessing.plot_marker(
        model = haldane_finite,
        marker = obc_marker,
        savefig = True,
        filename = 'obc_marker.pdf',

        # For better visualization
        vmin = -1, vmax = 1,
    )
```
Have a look at the two maps for a converged value of `L`, the only difference in the two are the boundary conditions requiring two different formulas for the local topological marker. Are your observations consistent with the expectations discussed in **part A**?

> **Hint:** You can increase the size of the points in the scatter plot by passing `s=[int]` to `plot_marker(...)`.

---

## Problem 2: Adding disorder
The single-point invariant and local topological markers do not rely on translational invariance and can therefore also be used to study disordered systems. Here, you will study how the topology is affected by on-site random (Anderson) disorder and vacancies. In particular, you will use local markers and single-point invariants to understand how disorder can modify the global and local topology of the system.

> **Note:** Playing with disorder accounts to use random numbers, so you may have to run the same code a few times to see the trend. Especially for the topological phase diagrams as a function of disorder, consider performing averages over disorder realizations to converge.

### Part A
We start from a topological phase, using for instance $\Delta=0$, $t=1$, $t_2=0.2$, $\phi=\pi/8$. Anderson disorder is a uniformly distributed random on-site term to be added to the Hamiltonian whose amplitude is $W$
$$\mathcal H = \mathcal H_0+\sum_iw_ic^{\dagger}_ic_i\qquad w_i\in\left[-\frac{W}{2},\frac{W}{2}\right]$$
and can be introduced in the Haldane model (both in the supercell and in the finite system) by using
```python
haldane_supercell.add_disorder_uniform(w = w)
```
Now, add disorder to the Haldane model in the supercell, try a few values ranging, for instance, from `w=0.1` to `w=10`. What happens to the topological phase as the disorder increases?

### Part B
For the same values of disorder `w`, plot also the PBC local Chern marker. What can you evince from those plots, based on the previous point? Does the picture change by going to open boundary conditions and bounded samples?

### Part C
Adding disorder to a system is not always detrimental to the topological phase. Construct an Haldane model in a trivial phase but close to the topological phase transition, for instance with parameters $\Delta=1.63$, $t=1$, $t_2=0.3$ and $\phi=-\pi/2$. Again, add Anderson disorder ranging from `w=0`to `w=10` and compute local markers and single-point invariants. What do you observe?

---

## Problem 3: Final take-aways (optional)

This final section requires more computational resources than the previous ones, thus is left as optional.

### Part A

Repeat the previous analysis introducing vacancies. To add vacancies into the model, we first create a list of `n_vac` unique vacancies in a $L\times L$ supercell/finite model, and then pass it to the model to remove those sites
```python
vacancies_list = utils.unique_vacancies(num = n_vac, Lx = L, Ly = L, basis = 2)
haldane_supercell.add_vacancies(vacancies_list)
haldane_finite.add_vacancies(vacancies_list)
```
For example, you may use the parameters $\Delta=0$, $t=1$, $t_2=0.2$, $\phi=-\pi/8$, $L=12$, and vary `n_vac` from $0$ to $100$ (density of vacancies $35$%). Do your conclusions still hold?

> **Note:** In StraWBerryPy, adding vacancies permanently modifies the model instance. If you wish to test a different number of vacancies, you must recreate the model from scratch.

> **Note:** Here, large system sizes are needed to properly converge the single-point Chern number as the gap of the Hamiltonian can be small. Averaging the topological invariant over multiple disorder realization is also recommended.

### Part B

Would you expect the topological phase transition as a function of the amplitude of Anderson disorder $W$ to be sharp or smooth? Why? You can test your intuition by running size scaling of the topological phase diagram by increasing the system size `L`. What do you see?

<div align="right">
&nbsp;

&nbsp;
---
**Hands-on by Nicolas Baù**

Thursday, 4 June 2026 - Trieste, Italy

</div>
