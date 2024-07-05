#!/bin/bash
#SBATCH -A m2651
#SBATCH -J etch_mem
#SBATCH -o etch_mem.log
#SBATCH -C cpu
#SBATCH --qos=shared
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1          # Number of tasks (MPI processes) per node
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=32GB
##SBATCH --mail-user=tommy_lin@berkeley.edu
##SBATCH --mail-type=ALL

dir="CALCS_diam15.4_thick3.5/-4.0_-4.0/"

rm -f "${dir}mprofile_output.dat"
mprof run --output "${dir}mprofile_output.dat" main.py "${dir}" > "${dir}run.dat"
mprof plot -o "${dir}mprofile.png" "${dir}mprofile_output.dat"