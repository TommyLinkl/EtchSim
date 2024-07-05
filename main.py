import sys, time, os
import json
import pickle
from src.constants import *
from src.sites import initialize_wz_lattice, initialize_hex_NPL, project_to_XY, write_projXY
from src.kmc import kmc_run

def load_input_from_json(filename):
    with open(filename, 'r') as file:
        loaded_data = json.load(file)
    
    sim_params.update(loaded_data.get("sim_params", {}))
    npl_params.update(loaded_data.get("NPL_setting", {}))
    calc_setting.update(loaded_data.get("calc_setting", {}))


def main():
    if len(sys.argv) != 2:
        print("Usage: python run_calc.py CALC_DIR/")
        return
    
    np.random.seed(32)
    sim_params['calc_dir'] = sys.argv[1]
    load_input_from_json(f"{sys.argv[1]}input.json")
    # print(sim_params)
    # print(npl_params)

    if npl_params['read_sites_from'] is None: 
        buffer = npl_params['buffer']      # all in AA
        NPL_hex_diameter = npl_params['NPL_hex_diameter']
        NPL_thickness = npl_params['NPL_thickness']
        wz_lattice = initialize_wz_lattice(NPL_hex_diameter+buffer, NPL_thickness+buffer, sim_params, verbosity=calc_setting['verbosity'])
        initialize_hex_NPL(wz_lattice, NPL_hex_diameter, NPL_thickness, sim_params, verbosity=calc_setting['verbosity'])

        sites_pklFileNormPath = os.path.normpath(os.path.join(sim_params['calc_dir'], sites_pklFileRelPath))
        print(f"\nWe are writing the wz_lattice sites list to file: {sites_pklFileNormPath}\n")
        with open(sites_pklFileNormPath, 'wb') as f:
            pickle.dump(wz_lattice, f)

    else: 
        sites_pklFileNormPath = os.path.normpath(os.path.join(sim_params['calc_dir'], npl_params['read_sites_from']))
        print(f"\nWARNING: We are reading the site_list from file {sites_pklFileNormPath}. Please make sure that this desired. NPL_setting -> 'buffer', 'NPL_hex_diameter', 'NPL_thickness' are all ignored. \n")
        with open(sites_pklFileNormPath, 'rb') as f:
            wz_lattice = pickle.load(f)

    # map to siteXY_list
    start_time = time.time()
    wz_lattice_XY = project_to_XY(wz_lattice)
    end_time = time.time()
    print(f"Project down to XY, elapsed time: {(end_time - start_time)*1000:.2f}ms. ")
    write_projXY(wz_lattice_XY, writeXY_filename=f"{sim_params['calc_dir']}init_projXY_atoms.xyz", mode='all_atoms')
    write_projXY(wz_lattice_XY, writeXY_filename=f"{sim_params['calc_dir']}init_projXY_atoms_a1.xyz", mode='a1_label')
    write_projXY(wz_lattice_XY, writeXY_filename=f"{sim_params['calc_dir']}init_projXY_atoms_a2.xyz", mode='a2_label')
    write_projXY(wz_lattice_XY, writeXY_filename=f"{sim_params['calc_dir']}init_projXY_atoms_a3.xyz", mode='a3_label')

    start_kmc_time = time.time()
    if not calc_setting['write_traj']: 
        trajFileName = None
    else: 
        trajFileName = f"{sim_params['calc_dir']}traj.xyz"
    kmc_run(wz_lattice, wz_lattice_XY, sim_params, trajFileName, calc_setting['write_every'], runtime_flag=calc_setting['runtime_flag'])
    end_kmc_time = time.time()
    print(f"\nDone with all KMC steps. Total elapsed time: {(end_kmc_time - start_kmc_time):.2f}s")


######################################################
if __name__ == "__main__":
    main()