#!/bin/bash

dir="CALCS/-4.0_-4.0/"

rm -f "${dir}mprofile_output.dat"
mprof run --output "${dir}mprofile_output.dat" run_calcs.py "${dir}" > "${dir}run.dat"
mprof plot -o "${dir}mprofile.png" "${dir}mprofile_output.dat"