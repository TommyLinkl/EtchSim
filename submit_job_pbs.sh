#!/bin/bash
#PBS -N etch
#PBS -l nodes=1:ppn=1:turtle
#PBS -q batch
##PBS -m abe
##PBS -M tommy_lin@berkeley.edu

# source /opt/intel/oneapi/setvars.sh
# source ~/.bashrc

calcGroupDir="CALCS_diam15.4_thick3.5/"
calcName="PleaseREPLACE/"     # "-4.0_-4.0/"   # "test/"
homeDir="$calcGroupDir$calcName"

scratchDir="/scratch/tommylin/EtchSim/$calcGroupDir$calcName"
mkdir -p "$scratchDir"
cp -r "$homeDir"* "$scratchDir"
cp "$calcGroupDir"init* "$SCRATCH/EtchSim/$calcGroupDir"
cp "$calcGroupDir"sites.pkl "$SCRATCH/EtchSim/$calcGroupDir"

nohup python main.py "$scratchDir" > "$scratchDir/run.dat" &
wait
cp -r "$scratchDir"/* "$homeDir"