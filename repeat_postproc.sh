#!/bin/bash

for i in $(seq 1 16); do
    new_calcName="-2.0_-2.0_repeat_${i}"
    echo $new_calcName
    
    sed "s/PleaseREPLACE/${new_calcName}/g" submit_postproc.sh > submit_postproc_repeat_${i}.sh
    sbatch submit_postproc_repeat_${i}.sh
    rm submit_postproc_repeat_${i}.sh
done
