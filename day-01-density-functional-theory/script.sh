#!/usr/bin/env bash
# ^^^ This line tells the system to run this file using bash

# Stop the script immediately if any command fails
set -e

# --- Input parameters ---
ecuts="20 30"  # Space-separated list of plane-wave cutoffs (in Ry) to try

# --- Directory paths ---
tmp_dir="./tmp"                 # Temporary data created during calculations
pseudo_dir="./pseudopotentials" # Where pseudopotential files are stored
out_dir="./out"                 # Where input/output files will be saved

# The command used to run Quantum ESPRESSO's pw.x
pw_launch="pw.x"

# Create directories if they don't already exist
# (-p means "no error if it already exists, and create parent dirs if needed")
mkdir -p "$tmp_dir"
mkdir -p "$out_dir"

# Loop over each plane-wave cutoff value
# (bash splits the string on spaces, so this iterates: 20, then 30)
for ecut in $ecuts; do

   # Build the input/output file paths for this cutoff
   input="$out_dir/NaCl.ecut=${ecut}.in"
   output="$out_dir/NaCl.ecut=${ecut}.out"

   # Write the input file with the current settings
   # (the values of $ecut, $pseudo_dir, and $tmp_dir are substituted automatically)
   cat > "$input" << EOF
&CONTROL
   calculation      = 'scf'
   tstress          = .true.
   tprnfor          = .true.
   outdir           = '$tmp_dir'
   prefix           = 'NaCl.$ecut'
   pseudo_dir       = '$pseudo_dir'
/
&SYSTEM
   ibrav            = 2
   ecutwfc          = $ecut
   ntyp             = 2
   nat              = 2
   celldm(1)        = 10.6
/
&ELECTRONS
   conv_thr         = 1e-08
   mixing_mode      = 'plain'
   mixing_beta      = 0.7
   diagonalization  = 'david'
/
&IONS
/
&CELL
/

ATOMIC_SPECIES
Na 22.98976928 Na_pseudo_dojo_v0.5.upf
Cl 35.45 Cl_pseudo_dojo_v0.5.upf

K_POINTS automatic
4 4 4  0 0 0

ATOMIC_POSITIONS crystal
Na -0.0000000000 0.0000000000 -0.0000000000 
Cl 0.5000000000 0.5000000000 0.5000000000
EOF

   # Run the calculation
   # '<' feeds the input file to pw.x (stdin redirection)
   # '>' sends the output to a file (stdout redirection)
   echo "Running: $pw_launch < $input > $output"
   $pw_launch < "$input" > "$output"

done

# Clean up temporary files created during calculations
rm -r "$tmp_dir"/*
