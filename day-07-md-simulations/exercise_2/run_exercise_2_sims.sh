#!/bin/zsh
set -euo pipefail

cd "$(dirname "$0")"

mkdir -p simulations/spcfw_prod
cd simulations/spcfw_prod
lmp -in ../../in.spcfw.prod.lammps
cd ../..

mkdir -p simulations/tip3p_prod
cd simulations/tip3p_prod
lmp -in ../../in.tip3p.prod.lammps
cd ../..
