import os
import json

large_calc_names = ['CALCS_defect1_2layer', 'CALCS_defect2_2layer', 'CALCS_defect3_2layer', 'CALCS_diam15.4_thick3.5']
chemPot_settings = ['-2.0', '-2.5', '-3.0', '-3.5', '-4.0']

calc_names = []
for i in large_calc_names:
    for j in chemPot_settings: 
        for k in range(1, 17): 
            calc_names.append(f"{i}/{j}_{j}_repeat_{k}/")
# print(calc_names)
invalid_calcs = []

for calc_name in calc_names:
    json_file = f"{calc_name}stats.json"

    if not os.path.exists(json_file) or os.path.getsize(json_file) == 0:
        invalid_calcs.append(calc_name)  # Add calc_name to invalid list
    else:
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
                if not data: 
                    invalid_calcs.append(calc_name)
        except json.JSONDecodeError:
            invalid_calcs.append(calc_name)

print(f"Calculations to re-do:\n{invalid_calcs}")
