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

nohup python postproc.py "$homeDir" > "$homeDir/run_postproc.dat" &
wait