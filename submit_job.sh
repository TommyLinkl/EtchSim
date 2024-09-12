#!/bin/bash
#SBATCH -A m2651
#SBATCH -J etch
#SBATCH -o etch.log
#SBATCH -C cpu
#SBATCH --qos=shared
#SBATCH -t 8:00:00
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

scratchDir="$SCRATCH/EtchSim/$calcGroupDir$calcName"
mkdir -p "$scratchDir"
cp -r "$homeDir"* "$scratchDir"
cp "$calcGroupDir"init* "$SCRATCH/EtchSim/$calcGroupDir"
cp "$calcGroupDir"sites.pkl "$SCRATCH/EtchSim/$calcGroupDir"

nohup python main.py "$scratchDir" > "$scratchDir/run.dat" &
wait
cp -r "$scratchDir"/* "$homeDir"