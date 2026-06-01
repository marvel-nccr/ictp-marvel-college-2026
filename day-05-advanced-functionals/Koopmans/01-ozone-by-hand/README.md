# Exercise 1: A screening parameter of ozone by hand

**Tutors**: Nicola Colonna and Edward Linscott

In this exercise you will perform a Koopmans-compliant KI calculation for the ozone molecule *by hand*, running each of the underlying `kcp.x` calculations yourself instead of letting the `koopmans` package orchestrate them.

The aim is to demystify what a Koopmans calculation actually does: you will compute the optimal screening parameter $\alpha$ for the HOMO of ozone, and then use it to run a final KI calculation. Exercise 2 repeats the same physics with the `koopmans` package driving everything automatically — so keep this exercise in mind as a reference for what is happening "under the hood" there.

## Files provided

- [`ozone_dft.in`](ozone_dft.in) — `kcp.x` input for the neutral ($N$-electron) DFT calculation
- [`ozone_dft_n-1.in`](ozone_dft_n-1.in) — `kcp.x` input for the constrained ($N{-}1$-electron) DFT calculation
- [`ozone_ki.in`](ozone_ki.in) — `kcp.x` input for the trial KI calculation
- [`get_alpha.sh`](get_alpha.sh) — a shell script that applies the screening formula from the lecture

> [!NOTE]
> All four calculations share the same `prefix` (`kc`) and `outdir` (`TMP/`), and communicate through restart files labelled by the `ndr`/`ndw` units. Run the steps below **in order** — each one reads the wavefunctions written by the previous one.

## Problem 1: Setting up the environment

Check that `kcp.x` is found:

```bash
which kcp.x
```

If `which` prints nothing, try running `workon koopmans`. If, after that, `which` still prints nothing, ask a tutor before continuing.

## Problem 2: The neutral ($N$-electron) DFT calculation

Run the DFT calculation for the neutral $N$-electron ozone molecule:

```bash
mpirun kcp.x -in ozone_dft.in | tee ozone_dft.out
```

Check that the calculation completed successfully before continuing.

### Part A

Open `ozone_dft.in` and inspect the `&SYSTEM` block. The molecule has 18 valence electrons (`nelec = 18`), split evenly between the two spin channels. What does `do_orbdep = .false.` mean? Why do we want to run such a calculation?

<details>
<summary><b>Solution</b></summary>

`do_orbdep = .false.` switches off the orbital-dependent part of the functional, leaving only the standard semi-local DFT (here, PBE) exchange–correlation.

We run a DFT calculation to give us the initial density and variational orbitals from which to start the Koopmans calculation. Doing this within a DFT framework before switching over to solving an ODDFT is more efficient. (This is especially true for the KI functional that we will use here, for which these will also be the _optimal_ density and variational orbitals!)

</details>

### Part B

After the neutral calculation finishes, copy the Kohn–Sham orbitals so that they can be reused as the initial *variational* orbitals of the KI calculation:

```bash
cp TMP/kc_90.save/K00001/evc1.dat        TMP/kc_90.save/K00001/evc01.dat
cp TMP/kc_90.save/K00001/evc2.dat        TMP/kc_90.save/K00001/evc02.dat
cp TMP/kc_90.save/K00001/evc_empty1.dat  TMP/kc_90.save/K00001/evc0_empty1.dat
cp TMP/kc_90.save/K00001/evc_empty2.dat  TMP/kc_90.save/K00001/evc0_empty2.dat
```

For molecules, the Kohn–Sham orbitals are used directly as the variational orbitals.

## Problem 3: The charged ($N{-}1$-electron) DFT calculation

Now run the DFT calculation for the charged system:

```bash
mpirun kcp.x -in ozone_dft_n-1.in | tee ozone_dft_n-1.out
```

Again, make sure the calculation finishes correctly.

### Part A

Open `ozone_dft_n-1.in` and look at the three keywords in the `&SYSTEM` block:

```fortran
fixed_state  = .true.
fixed_band   = 9
f_cutoff     = 1e-05
```

`f_cutoff` corresponds to the occupation of the variational orbital with index `fixed_band`.

What is this calculation doing?

<details>
<summary><b>Solution</b></summary>

The 9th variational orbital is the HOMO. We are setting its occupation to 0.

That means that this is the input for a $N-1$-electron calculation where we have constrained the HOMO of the $N$-electron solution to be empty.

</details>

### Part B

This calculation restarts from the neutral run (`restart_mode = 'restart'`, `ndr = 90`). From `ozone_dft.out` and `ozone_dft_n-1.out`, read off the total energies $E^\text{DFT}[N]$ and $E^\text{DFT}[N{-}1]$. Their difference is a ΔSCF estimate of the ionisation potential — what value do you get (in eV)?

<details>
<summary><b>Solution</b></summary>

You should get something close to 12.49 eV

</details>

## Problem 4: The trial KI calculation

The KI eigenvalue depends on the screening parameter $\alpha$. To pin down the optimal $\alpha$ we first need *one* trial KI calculation at a guessed value.

Open `ozone_ki.in` and replace the `<alpha>` placeholder in the `&NKSIC` block with the trial value `0.7`:

```fortran
&NKSIC
   nkscalfact         = 0.7
   which_orbdep       = 'nki'
   do_innerloop       = .false.
   esic_conv_thr      = 1.8000000000000002e-08
   do_innerloop_empty = .false.
/
```

Then run the trial KI calculation:

```bash
mpirun kcp.x -in ozone_ki.in | tee ozone_ki.out
```

### Part A

This input has `do_orbdep = .true.` and `which_orbdep = 'nki'`, whereas `ozone_dft.in` had `do_orbdep = .false.`. What do these keywords switch on?

<details>
<summary><b>Solution</b></summary>

Together they switch on the orbital-density-dependent KI correction: `do_orbdep = .true.` enables the orbital-dependent term in the functional, and `which_orbdep = 'nki'` picks the *Koopmans-integral* (KI) form for it.

</details>

### Part B

Compare the HOMO eigenvalue in `ozone_ki.out` with the one in `ozone_dft.out`. The KI eigenvalue at $\alpha = 0$ is, by construction, equal to the DFT eigenvalue. How has the trial correction (at $\alpha_0 = 0.7$) shifted the HOMO?

<details>
<summary><b>Solution</b></summary>

The KI HOMO is pushed *down* (more negative) compared to the PBE value, so the predicted ionisation potential $-\varepsilon_\text{HOMO}$ goes up. This is the expected qualitative effect of the Koopmans correction: it counteracts the self-interaction error that makes the PBE HOMO too shallow.

</details>

## Problem 5: Extracting the optimal screening parameter

The optimal $\alpha$ is the one that enforces the Koopmans condition on the HOMO,

$$\varepsilon^\text{KI}_\text{HOMO}(\alpha_\text{opt}) = E^\text{DFT}[N] - E^\text{DFT}[N{-}1].$$

(The right-hand side is minus the ΔSCF ionisation potential from Problem 3 — the HOMO eigenvalue should equal the energy *change* on removing an electron, which is negative.)

### Part A

Derive the equation for $\alpha_\text{opt}$, using the facts that...
- $\varepsilon^\text{KI}_\text{HOMO}(\alpha)$ is linear in $\alpha$ (why?)
- you have at your disposal $`\varepsilon^\text{KI}_\text{HOMO}(0) = \varepsilon^\text{DFT}_\text{HOMO}`$ and $\varepsilon^\text{KI}_\text{HOMO}(\alpha_0)$
- for $\alpha_\text{opt}$, the Koopmans condition above is satisfied

<details>
<summary><b>Solution</b></summary>

The KI correction does not change the ground-state density nor the occupied variational orbitals. Therefore the only  $\alpha$-dependence of the corrected eigenvalue is the explicit prefactor:

$`\varepsilon^\text{KI}_\text{HOMO}(\alpha)=\varepsilon^\text{DFT}_\text{HOMO} + \alpha\lambda_{\text{HOMO}}[\rho,\rho_\text{HOMO}],`$

with $\lambda_{\text{HOMO}}[\rho,\rho_\text{HOMO}]$ independent of $\alpha$. At $\alpha = 0$ the correction vanishes and we recover the DFT eigenvalue.

We know the line at $\alpha = 0$ and at $\alpha = \alpha_0$, so

$`\lambda_{\text{HOMO}} = \frac{\varepsilon^\text{KI}_\text{HOMO}(\alpha_0) - \varepsilon^\text{DFT}_\text{HOMO}}{\alpha_0}.`$

Setting $`\varepsilon^\text{KI}_\text{HOMO}(\alpha_\text{opt}) = E^\text{DFT}[N] - E^\text{DFT}[N{-}1]`$ and solving for $\alpha_\text{opt}$ we get

$`\alpha_\text{opt} = \alpha_0\,\frac{\big(E^\text{DFT}[N] - E^\text{DFT}[N{-}1]\big) - \varepsilon^\text{DFT}_\text{HOMO}}{\varepsilon^\text{KI}_\text{HOMO}(\alpha_0) - \varepsilon^\text{DFT}_\text{HOMO}}.`$

</details>

### Part B

Now plug numbers into the equation you derived in Part A. Read $`\varepsilon^\text{DFT}_\text{HOMO}`$ and $`\varepsilon^\text{KI}_\text{HOMO}(\alpha_0)`$ from `ozone_dft.out` and `ozone_ki.out`, use the ΔSCF target from Problem 3, and work out $\alpha_\text{opt}$ by hand.

<details>
<summary><b>Solution</b></summary>

First convert the ΔSCF target to eV:

$$E^\text{DFT}[N] - E^\text{DFT}[N{-}1] = (-47.5296) - (-47.0705) \,\text{Ha} = -.4591\,\text{Ha} \approx -12.49\,\text{eV}.$$

Then plug into the formula from Part A with $\alpha_0 = 0.7$:

$$\alpha_\text{opt} = 0.7 \times \frac{(-12.49) - (-7.92)}{(-12.00) - (-7.92)} = 0.7 \times \frac{-4.57}{-4.08} \approx 0.784.$$

</details>

### Part C

Verify your value by running the helper script, which applies exactly this formula to the output files:

```bash
sh get_alpha.sh
```

## Problem 6: The final KI calculation

Create a new input file `ozone_ki_opt.in` by copying `ozone_ki.in` and making two changes:

1. Update the restart units so that the final run reads the trial KI wavefunctions and writes to a fresh unit:

   ```fortran
   ndr = 90
   ndw = 92
   ```

2. Replace `nkscalfact` with the optimal value from Problem 5:

   ```fortran
   &NKSIC
      nkscalfact         = <optimal_alpha>
      which_orbdep       = 'nki'
      do_innerloop       = .false.
      esic_conv_thr      = 1.8000000000000002e-08
      do_innerloop_empty = .false.
   /
   ```

Then run the final KI calculation at the optimal screening parameter:

```bash
mpirun kcp.x -in ozone_ki_opt.in | tee ozone_ki_opt.out
```

### Part A

Read the HOMO eigenvalue from `ozone_ki_opt.out`. Does it now satisfy the Koopmans condition — that is, does $-\varepsilon^\text{KI}_\text{HOMO}$ agree with the ΔSCF ionisation potential from Problem 3?

<details>
<summary><b>Solution</b></summary>

It should agree. (Remember to convert to the correct units!)

</details>

### Part B

Compare the KI ionisation potential against the DFT (PBE) value and the experimental value of ozone, IP ≈ 12.5 eV.

## Problem 7: Take-aways

In this exercise you computed a *single* screening parameter, for the HOMO, by running three explicit `kcp.x` calculations (neutral DFT, constrained $N{-}1$ DFT, trial KI) and applying the same screening parameter to every variational orbital in the system.

A full Koopmans calculation needs one screening parameter for *every* orbital, and the trial-KI/constrained-DFT loop is iterated to self-consistency. How many calculations would that require for ozone's ten orbitals?

<details>
<summary><b>Solution</b></summary>

- 1 for the DFT initialzation
- 10 for the self-consistent loop (1 per orbital)
- 1 for the final KI calculation
... for a total of 12 calculations if we don't do any self-consistency, 22 for two iterations, _etc._


> [!NOTE]
> Self-consistency with KI is unnecessary, and we can often take advantage of symmetries to not have to compute from scratch screening parameters for every orbital.
</details>
<br>

This is exactly the bookkeeping that the `koopmans` package automates — which is what you will see in Exercise 2.
