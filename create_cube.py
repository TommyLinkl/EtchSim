import sys, time, os
import json
import pickle
from src.constants import *
from src.sites import initialize_wz_lattice, aggregate_to_xyz, update_whole_lattice_iteration, NPL_within_sphere, NPL_within_cube, check_neighbors
from main import load_input_from_json

def create_cube():
    if len(sys.argv) != 2:
        print("Usage: python create_cube.py CALC_DIR/")
        return
    
    calc_setting['calc_dir'] = sys.argv[1]
    load_input_from_json(f"{sys.argv[1]}input.json")

    # initialize_wz_lattice -> cube
    buffer = npl_params['buffer']      # all in AA
    cube_edge = npl_params["cube_edge"]
    wz_lattice = initialize_wz_lattice(cube_edge+buffer, cube_edge+buffer, sim_params, verbosity=calc_setting['verbosity'])
    

    # initialize_hex_NPL -> cube
    start_time = time.time()
    for site in wz_lattice:
        if NPL_within_cube(site.real_space_coord, cube_edge):
            site.has_atom = True
            site.iteration_changed = True

            # Change its neighbors' iteration_changed to True
            for neighbor_idx in site.neighbor_sites_idx:
                wz_lattice[neighbor_idx].iteration_changed = True

    update_whole_lattice_iteration(wz_lattice, sim_params)

    end_time = time.time()
    print(f"Initialize NPL atoms, elapsed time: {(end_time - start_time):.2f}s")

    check_neighbors(wz_lattice, just_atoms=True)

    aggregate_to_xyz(wz_lattice, write_site=True, write_atoms=True, write_filename=f"{calc_setting['calc_dir']}init_lattice_atoms.xyz")
    aggregate_to_xyz(wz_lattice, write_site=False, write_atoms=True, write_filename=f"{calc_setting['calc_dir']}init_atoms.xyz")


    sites_pklFileNormPath = os.path.normpath(os.path.join(calc_setting['calc_dir'], sites_pklFileRelPath))
    print(f"\nWe are writing the wz_lattice sites list to file: {sites_pklFileNormPath}\n")
    with open(sites_pklFileNormPath, 'wb') as f:
        pickle.dump(wz_lattice, f)


######################################################
if __name__ == "__main__":
    create_cube()