import os
import shutil
import json
import random

# Constants
source_dir = './CALCS_diam15.4_thick3.5/-4.0_-4.0'
destination_template = './CALCS_diam15.4_thick3.5/-4.0_-4.0_repeat_{}'
json_filename = 'input.json'

def update_json_file(file_path, random_seed):
    with open(file_path, 'r') as f:
        data = json.load(f)
    
    data['calc_setting']['random_seed'] = random_seed
    
    with open(file_path, 'w') as f:
        json.dump(data, f, indent=4)

def repeat_calcs(nRepeat):
    for i in range(1, nRepeat+1):
        destination_dir = destination_template.format(i)
        
        shutil.copytree(source_dir, destination_dir)
        
        random_seed = random.randint(0, 1000)
        json_file_path = os.path.join(destination_dir, json_filename)
        update_json_file(json_file_path, random_seed)

########################
repeat_calcs(32)
