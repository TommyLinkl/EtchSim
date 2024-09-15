#!/bin/bash

mode="check"
calc="-3.5_-3.5"
calcPrefix="CALCS_defect1_2layer/"

if [[ "$mode" == "check" ]]; then
    for i in {1..16}; do
        diff "${calcPrefix}${calc}_repeat_1/init_site_list.pkl" "${calcPrefix}${calc}_repeat_${i}/init_site_list.pkl"
        diff "${calcPrefix}${calc}_repeat_1/init_siteXY_list.pkl" "${calcPrefix}${calc}_repeat_${i}/init_siteXY_list.pkl"
        diff "${calcPrefix}${calc}_repeat_1/init_vacXY_list.pkl" "${calcPrefix}${calc}_repeat_${i}/init_vacXY_list.pkl"
    done
elif [[ "$mode" == "remove" ]]; then
    for i in {1..16}; do
        mv "${calcPrefix}${calc}_repeat_${i}/init_site_list.pkl" "${calcPrefix}${calc}/init_site_list.pkl"
        mv "${calcPrefix}${calc}_repeat_${i}/init_siteXY_list.pkl" "${calcPrefix}${calc}/init_siteXY_list.pkl"
        mv "${calcPrefix}${calc}_repeat_${i}/init_vacXY_list.pkl" "${calcPrefix}${calc}/init_vacXY_list.pkl"
    done
elif [[ "$mode" == "putBack" ]]; then
    for i in {1..16}; do
        cp "${calcPrefix}${calc}/init_site_list.pkl" "${calcPrefix}${calc}_repeat_${i}/init_site_list.pkl"
        cp "${calcPrefix}${calc}/init_siteXY_list.pkl" "${calcPrefix}${calc}_repeat_${i}/init_siteXY_list.pkl"
        cp "${calcPrefix}${calc}/init_vacXY_list.pkl" "${calcPrefix}${calc}_repeat_${i}/init_vacXY_list.pkl"
    done
fi