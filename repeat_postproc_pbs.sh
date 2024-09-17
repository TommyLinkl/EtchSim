#!/bin/bash

calcNames=(
    '-2.0_-2.0_repeat_10'
    '-2.0_-2.0_repeat_13'
    '-2.0_-2.0_repeat_16'
    '-2.0_-2.0_repeat_2'
    '-2.0_-2.0_repeat_3'
    '-2.0_-2.0_repeat_5'
    '-2.0_-2.0_repeat_6'
    '-2.0_-2.0_repeat_7'
    '-2.5_-2.5_repeat_13'
    '-3.0_-3.0_repeat_3'
)

# calcNames=()
# for i in $(seq 1 16); do
#     calcNames+=("-4.0_-4.0_repeat_${i}")
# done

for new_calcName in "${calcNames[@]}"; do
    echo $new_calcName
    
    # Replace placeholder and create job script, then submit it
    sed "s/PleaseREPLACE/${new_calcName}/g" submit_postproc_pbs.sh > submit_postproc_pbs_temp.sh
    qsub submit_postproc_pbs_temp.sh
    rm submit_postproc_pbs_temp.sh
done
