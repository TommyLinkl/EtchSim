#!/bin/bash

mode="remove"
calc="-4.0_-4.0"
calcPrefix="CALCS_diam15.4_thick3.5/"

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