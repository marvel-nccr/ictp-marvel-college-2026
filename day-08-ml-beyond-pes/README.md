# ML beyond the PES: Hands-on session

**Tutors**: Joseph W. Abbott, Matthias Kellner, Paolo Pegolo, Filippo Bigi

---

This folder contains three exercises on machine learning for targets beyond the
potential energy surface.

## Setup

Ensure you have updated your virtual machine to pull the latest changes. Then, enter the
directory for today's exercises:

```bash
update
cdd 8
```

Each exercise requires its own environment to avoid any versioning conflict. At the top
of each exercise markdown file there will be an indication about which environment to
use.

For exercises 1 and 2, use the `metatomic` environment:

```bash
cd exercise_1  # or exercise_2
workon metatomic
```

For exercise 3, the `nmr` environment is needed:

```bash
cd exercise_3
workon nmr
```

JupyterLab can be launched from the terminal as:

```bash
jupyter lab
```

---

## Exercises

Access each exercise folder as instructed above. Then, follow the
instructions in `.md` file contained in that folder, with roughly follows the structrue:

1. Open the atomistic cookbook recipe correpsonding to that exercise in the browser on
   your laptop

    * Exercise 1: https://atomistic-cookbook.org/examples/water-ir-spectrum/water-ir-spectrum.html
    * Exercise 2: https://atomistic-cookbook.org/examples/ml-density/ml-density.html
    * Exercise 3: https://atomistic-cookbook.org/examples/shiftml-structure-match/shiftml-structure-match.html

2. Follow and read the recipe
3. Complete the follow up hands-on exercise contained in the relevant subdirectory of
   this folder.
