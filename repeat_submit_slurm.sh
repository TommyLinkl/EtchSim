#!/bin/bash

# calcNames=(
#     "-2.0_-2.0_repeat_1"
#     "-2.0_-2.0_repeat_2"
#     "-2.0_-2.0_repeat_3"
# )

calcNames=()
for i in $(seq 1 16); do
    calcNames+=("-2.0_-2.0_repeat_${i}")
done

for new_calcName in "${calcNames[@]}"; do
    echo $new_calcName
    
    # Replace placeholder and create job script, then submit it
    sed "s/PleaseREPLACE/${new_calcName}/g" submit_job.sh > submit_job_temp.sh
    sbatch submit_job_temp.sh
    rm submit_job_temp.sh
done
