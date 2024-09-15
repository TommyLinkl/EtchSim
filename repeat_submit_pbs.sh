#!/bin/bash

calcNames=(
    "-4.0_-4.0_repeat_10"
    "-4.0_-4.0_repeat_11"
    "-4.0_-4.0_repeat_12"
    "-4.0_-4.0_repeat_13"
    "-4.0_-4.0_repeat_14"
    "-4.0_-4.0_repeat_15"
)

# calcNames=()
# for i in $(seq 1 16); do
#     calcNames+=("-4.0_-4.0_repeat_${i}")
# done

for new_calcName in "${calcNames[@]}"; do
    echo $new_calcName
    
    # Replace placeholder and create job script, then submit it
    sed "s/PleaseREPLACE/${new_calcName}/g" submit_job_pbs.sh > submit_job_pbs_temp.sh
    qsub submit_job_pbs_temp.sh
    rm submit_job_pbs_temp.sh
done
