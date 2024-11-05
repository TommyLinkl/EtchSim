#!/bin/bash

calcNames=(
    "-1.7_-1.7_repeat_1"
    "-1.8_-1.8_repeat_1"
    "-1.9_-1.9_repeat_1"
    "-2.0_-2.0_repeat_1"
    "-2.1_-2.1_repeat_1"
    "-2.2_-2.2_repeat_1"
    "-2.3_-2.3_repeat_1"
    "-2.4_-2.4_repeat_1"
    "-2.5_-2.5_repeat_1"
    "-2.6_-2.6_repeat_1"
    "-2.7_-2.7_repeat_1"
    "-2.8_-2.8_repeat_1"
    "-2.9_-2.9_repeat_1"
    "-3.0_-3.0_repeat_1"
    "-3.1_-3.1_repeat_1"
    "-3.5_-3.5_repeat_1"
    "-4.0_-4.0_repeat_1"
)

# calcNames=()
# for i in $(seq 1 8); do
#     calcNames+=("-2.0_-3.0_repeat_${i}")
#     calcNames+=("-2.0_-4.0_repeat_${i}")
#     calcNames+=("-3.0_-4.0_repeat_${i}")
# done

calcGroupDir="CALCS_lowTemp"

for new_calcName in "${calcNames[@]}"; do
    echo $new_calcName
    
    # Replace placeholder and create job script, then submit it
    # sed "s/PleaseREPLACE/${new_calcName}/g" submit_job_pbs.sh > submit_job_pbs_temp.sh
    # qsub submit_job_pbs_temp.sh
    # rm submit_job_pbs_temp.sh
    
    # Create a custom SLURM script for each job name
    cat > submit_job_pbs_temp.sh << EOF
#!/bin/bash
#PBS -N etch_${calcGroupDir}_${calcName}
#PBS -l nodes=1:ppn=2:turtle
#PBS -q batch
#PBS -j oe
#PBS -o etch_${calcGroupDir}_${calcName}.log
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

scratchDir="/scratch/tommylin/EtchSim/$calcGroupDir/$calcName/"
mkdir -p "$scratchDir"
cp -r "$homeDir"input.json "$scratchDir"
# cp "$PBS_O_WORKDIR/$calcGroupDir/"init* "/scratch/tommylin/EtchSim/$calcGroupDir/"
cp "$PBS_O_WORKDIR/$calcGroupDir/"*.pkl "/scratch/tommylin/EtchSim/$calcGroupDir/"

nohup python main.py "$scratchDir" > "$scratchDir"run.dat 2>&1 &
wait
cp -r "$scratchDir"* "$homeDir"

end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Total runtime: $runtime seconds"
EOF

    qsub submit_job_pbs_temp.sh
    rm submit_job_pbs_temp.sh

done
