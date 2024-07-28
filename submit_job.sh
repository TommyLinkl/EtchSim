#!/bin/bash
#SBATCH -A m2651
#SBATCH -J etch
#SBATCH -o etch.log
#SBATCH -C cpu
#SBATCH --qos=shared
#SBATCH -t 3:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1          # Number of tasks (MPI processes) per node
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=32GB
##SBATCH --mail-user=tommy_lin@berkeley.edu
##SBATCH --mail-type=ALL

source ~/.bashrc

dir="CALCS_diam15.4_thick3.5/test/"  # "CALCS_diam15.4_thick3.5/-4.0_-4.0/"

scratch_dir="$SCRATCH/EtchSim/$dir"
echo $scratch_dir
mkdir -p "$scratch_dir"

nohup python main.py "${scratch_dir}" > "${dir}run.dat" &
# nohup python postproc.py "${dir}" > "${dir}run.dat" &
wait
cp -r "$scratch_dir"/* "$dir"