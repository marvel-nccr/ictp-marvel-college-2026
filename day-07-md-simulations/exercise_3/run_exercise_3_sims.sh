#!/bin/zsh
set -euo pipefail

cd "$(dirname "$0")"

mkdir -p simulations/spcfw_nvt_nh
cd simulations/spcfw_nvt_nh
lmp -in ../../in.spcfw.nvt_nh.lammps
cd ../..

mkdir -p simulations/spcfw_nvt_langevin
cd simulations/spcfw_nvt_langevin
lmp -in ../../in.spcfw.nvt_langevin.lammps
cd ../..

mkdir -p simulations/tip3p_nvt_nh
cd simulations/tip3p_nvt_nh
lmp -in ../../in.tip3p.nvt_nh.lammps
cd ../..

mkdir -p simulations/tip3p_nvt_langevin
cd simulations/tip3p_nvt_langevin
lmp -in ../../in.tip3p.nvt_langevin.lammps
cd ../..
