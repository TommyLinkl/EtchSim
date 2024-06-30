import sys
import json
from src.constants import *
from src.sites import initialize_wz_lattice, initialize_hex_NPL
from src.kmc import kmc_step, kmc_run

def load_input_from_json(filename):
    with open(filename, 'r') as file:
        loaded_data = json.load(file)
    
    sim_params.update(loaded_data.get("sim_params", {}))
    npl_params.update(loaded_data.get("NPL_setting", {}))


def main():
    if len(sys.argv) != 2:
        print("Usage: python run_calc.py CALC_DIR/")
        return
    
    np.random.seed(32)
    sim_params['calc_dir'] = sys.argv[1]
    load_input_from_json(f"{sys.argv[1]}input.json")
    # print(sim_params)
    # print(npl_params)

    buffer = npl_params['buffer']      # all in AA
    NPL_hex_diameter = npl_params['NPL_hex_diameter']
    NPL_thickness = npl_params['NPL_thickness']
    wz_lattice = initialize_wz_lattice(NPL_hex_diameter+buffer, NPL_thickness+buffer, sim_params, verbosity=0)
    initialize_hex_NPL(wz_lattice, NPL_hex_diameter, NPL_thickness, sim_params, verbosity=0)

    kmc_run(wz_lattice, sim_params, f'{sim_params['calc_dir']}traj.xyz', runtime_flag=True)


######################################################
if __name__ == "__main__":
    main()