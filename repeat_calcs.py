import os
import shutil
import json
import random

def update_json_file(file_path, random_seed):
    with open(file_path, 'r') as f:
        data = json.load(f)
    
    data['calc_setting']['random_seed'] = random_seed
    
    with open(file_path, 'w') as f:
        json.dump(data, f, indent=4)

def repeat_calcs(nRepeat, source_dir, destination_template, json_filename='input.json'):
    for i in range(1, nRepeat+1):
        destination_dir = destination_template.format(i)
        
        shutil.copytree(source_dir, destination_dir)
        
        random_seed = random.randint(0, 1000)
        json_file_path = os.path.join(destination_dir, json_filename)
        update_json_file(json_file_path, random_seed)

########################
if __name__ == "__main__":
    for mu in [-1.7, -1.8, -1.9, -2.0, -2.1, -2.2, -2.3, -2.4, -2.5, -2.6, -2.7, -2.8, -2.9, -3.0, -3.1, -3.5, -4.0]: 
        source_dir = f"./CALCS_lowTemp/{mu:.1f}_{mu:.1f}"
        destination_template = source_dir + "_repeat_{}"

        repeat_calcs(1, source_dir, destination_template)


