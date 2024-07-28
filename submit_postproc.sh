#!/bin/bash
#SBATCH -A m2651
#SBATCH -J etch_postproc
#SBATCH -o etch_postproc.log
#SBATCH -C cpu
#SBATCH --qos=shared
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1          # Number of tasks (MPI processes) per node
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=32GB
##SBATCH --mail-user=tommy_lin@berkeley.edu
##SBATCH --mail-type=END,FAIL    # ALL, BEGIN

source ~/.bashrc

calcGroupDir="CALCS_diam15.4_thick3.5/"
calcName="PleaseREPLACE/"     # "-4.0_-4.0/"   # "test/"
homeDir="$calcGroupDir$calcName"

nohup python postproc.py "$homeDir" > "$homeDir/run_postproc.dat" &
wait