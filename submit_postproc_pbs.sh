#!/bin/bash
#PBS -N etch_postproc_diam15.4_thick3.5_PleaseREPLACE
#PBS -l nodes=1:ppn=4:turtle
#PBS -q batch
#PBS -j oe
#PBS -o etch_postproc_diam15.4_thick3.5_PleaseREPLACE.log
##PBS -m abe
##PBS -M tommy_lin@berkeley.edu

start_time=$(date +%s)

source ~/.bashrc
export PATH="/home/tommylin/miniconda3/bin:$PATH"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate myenv

echo $HOSTNAME

calcGroupDir="CALCS_diam15.4_thick3.5/"
calcName="PleaseREPLACE/"     # "-4.0_-4.0/"   # "test/"
echo $calcGroupDir$calcName

homeDir="$PBS_O_WORKDIR/$calcGroupDir$calcName"
cd $PBS_O_WORKDIR

nohup python postproc.py "$homeDir" > "$homeDir"run_postproc.dat 2>&1 &
wait

end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Total runtime: $runtime seconds"