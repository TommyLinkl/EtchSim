#!/bin/bash

mode="putBack"
calc="-2.0_-2.0"

if [[ "$mode" == "check" ]]; then
    for i in {1..32}; do
        diff "CALCS_diam15.4_thick3.5/${calc}_repeat_1/init_site_list.pkl" "CALCS_diam15.4_thick3.5/${calc}_repeat_${i}/init_site_list.pkl"
        diff "CALCS_diam15.4_thick3.5/${calc}_repeat_1/init_siteXY_list.pkl" "CALCS_diam15.4_thick3.5/${calc}_repeat_${i}/init_siteXY_list.pkl"
        diff "CALCS_diam15.4_thick3.5/${calc}_repeat_1/init_vacXY_list.pkl" "CALCS_diam15.4_thick3.5/${calc}_repeat_${i}/init_vacXY_list.pkl"
    done
elif [[ "$mode" == "remove" ]]; then
    for i in {1..32}; do
        mv "CALCS_diam15.4_thick3.5/${calc}_repeat_${i}/init_site_list.pkl" "CALCS_diam15.4_thick3.5/${calc}/init_site_list.pkl"
        mv "CALCS_diam15.4_thick3.5/${calc}_repeat_${i}/init_siteXY_list.pkl" "CALCS_diam15.4_thick3.5/${calc}/init_siteXY_list.pkl"
        mv "CALCS_diam15.4_thick3.5/${calc}_repeat_${i}/init_vacXY_list.pkl" "CALCS_diam15.4_thick3.5/${calc}/init_vacXY_list.pkl"
    done
elif [[ "$mode" == "putBack" ]]; then
    for i in {1..32}; do
        cp "CALCS_diam15.4_thick3.5/${calc}/init_site_list.pkl" "CALCS_diam15.4_thick3.5/${calc}_repeat_${i}/init_site_list.pkl"
        cp "CALCS_diam15.4_thick3.5/${calc}/init_siteXY_list.pkl" "CALCS_diam15.4_thick3.5/${calc}_repeat_${i}/init_siteXY_list.pkl"
        cp "CALCS_diam15.4_thick3.5/${calc}/init_vacXY_list.pkl" "CALCS_diam15.4_thick3.5/${calc}_repeat_${i}/init_vacXY_list.pkl"
    done
fi