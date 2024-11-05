#!/bin/bash
#PBS -N etch_lowTemp_PleaseREPLACE
#PBS -l nodes=1:ppn=2:turtle
#PBS -q batch
#PBS -j oe
#PBS -o etch_lowTemp_PleaseREPLACE.log
##PBS -m abe
##PBS -M tommy_lin@berkeley.edu

start_time=$(date +%s)

source ~/.bashrc
export PATH="/home/tommylin/miniconda3/bin:$PATH"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate myenv

echo $HOSTNAME

calcGroupDir="CALCS_lowTemp/"
calcName="PleaseREPLACE/"     # "-4.0_-4.0/"   # "test/"
echo $calcGroupDir$calcName

homeDir="$PBS_O_WORKDIR/$calcGroupDir$calcName"
cd $PBS_O_WORKDIR

scratchDir="/scratch/tommylin/EtchSim/$calcGroupDir$calcName"
mkdir -p "$scratchDir"
cp -r "$homeDir"input.json "$scratchDir"
# cp "$PBS_O_WORKDIR/$calcGroupDir"init* "/scratch/tommylin/EtchSim/$calcGroupDir"
cp "$PBS_O_WORKDIR/$calcGroupDir"*.pkl "/scratch/tommylin/EtchSim/$calcGroupDir"

nohup python main.py "$scratchDir" > "$scratchDir"run.dat 2>&1 &
wait
cp -r "$scratchDir"* "$homeDir"

end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Total runtime: $runtime seconds"