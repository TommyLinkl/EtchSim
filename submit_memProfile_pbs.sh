#!/bin/bash
#PBS -N etch
#PBS -l nodes=1:ppn=1:turtle
#PBS -q batch
##PBS -m abe
##PBS -M tommy_lin@berkeley.edu.edu
##SBATCH --mail-type=ALL

source ~/.bashrc
export PATH="/home/tommylin/miniconda3/bin:$PATH"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate myenv

dir="CALCS_diam15.4_thick3.5/-4.0_-4.0/"

rm -f "${dir}mprofile_output.dat"
mprof run --output "${dir}mprofile_output.dat" main.py "${dir}" > "${dir}run.dat"
mprof plot -o "${dir}mprofile.png" "${dir}mprofile_output.dat"