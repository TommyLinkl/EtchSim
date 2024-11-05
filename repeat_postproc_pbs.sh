#!/bin/bash

# calcNames=(
#     '-2.0_-2.0_repeat_10'
#     '-2.0_-2.0_repeat_13'
#     '-2.0_-2.0_repeat_16'
# )

calcNames=()
for i in $(seq 1 16); do
    calcNames+=("-2.5_-2.0_repeat_${i}")
    calcNames+=("-3.5_-2.0_repeat_${i}")
    calcNames+=("-3.5_-2.5_repeat_${i}")
    calcNames+=("-3.5_-3.0_repeat_${i}")
done

calcGroupDir="CALCS_uneven"

for calcName in "${calcNames[@]}"; do
    echo $calcName
    
    # Replace placeholder and create job script, then submit it
    # sed "s/PleaseREPLACE/${calcName}/g" submit_postproc_pbs.sh > submit_postproc_pbs_temp.sh
    # qsub submit_postproc_pbs_temp.sh
    # rm submit_postproc_pbs_temp.sh

    # Create a custom SLURM script for each job name
    cat > submit_postproc_pbs_temp.sh << EOF
#!/bin/bash
#PBS -N etch_postproc_${calcGroupDir}_${calcName}
#PBS -l nodes=1:ppn=2:turtle
#PBS -q batch
#PBS -j oe
#PBS -o etch_postproc_${calcGroupDir}_${calcName}.log
##PBS -m abe
##PBS -M tommy_lin@berkeley.edu

start_time=$(date +%s)

source ~/.bashrc
export PATH="/home/tommylin/miniconda3/bin:$PATH"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate myenv

echo $HOSTNAME
echo $calcGroupDir/$calcName/

homeDir="$PBS_O_WORKDIR/$calcGroupDir/$calcName/"
cd $PBS_O_WORKDIR

nohup python postproc.py "$homeDir" > "$homeDir"run_postproc.dat 2>&1 &
wait

end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Total runtime: $runtime seconds"
EOF

    qsub submit_postproc_pbs_temp.sh
    rm submit_postproc_pbs_temp.sh

done
