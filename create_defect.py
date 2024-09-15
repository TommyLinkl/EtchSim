import copy
import pickle
import gc
from src.constants import *
from src.sites import update_whole_lattice_iteration, aggregate_to_xyz

def create_defect(layer=1):
    # Reading in sites from the .pkl file
    with open('CALCS_defect_construction/sites.pkl', 'rb') as f:
        wz_lattice = pickle.load(f)

    if layer == 1: 
        z_cutoff = 16.0
    elif layer == 2: 
        z_cutoff = 13.0
    elif layer == 3: 
        z_cutoff = 9.5
    else: 
        raise ValueError("Layer input not allowed. ")

    ##########################################
    # defect1
    defect = copy.deepcopy(wz_lattice)
    for site in defect: 
        coord = site.real_space_coord
        if (coord[2] > z_cutoff) and (coord[1]>=np.sqrt(3)*coord[0]) and (coord[1]>=-np.sqrt(3)*coord[0]): 
            site.has_atom = False
            site.iteration_changed = True

    update_whole_lattice_iteration(defect, sim_params)
    aggregate_to_xyz(defect, write_site=False, write_atoms=True, write_filename=f"CALCS_defect_construction/defect1_{layer}layer.xyz")
    with open(f"CALCS_defect_construction/defect1_{layer}layer.pkl", 'wb') as f:
        pickle.dump(defect, f)

    del defect
    gc.collect()

    ##########################################
    # defect2
    defect = copy.deepcopy(wz_lattice)
    for site in defect: 
        coord = site.real_space_coord
        if (coord[2] > z_cutoff) and (coord[1]>=0) and (coord[1]>=-np.sqrt(3)*coord[0]): 
            site.has_atom = False
            site.iteration_changed = True

    update_whole_lattice_iteration(defect, sim_params)
    aggregate_to_xyz(defect, write_site=False, write_atoms=True, write_filename=f"CALCS_defect_construction/defect2_{layer}layer.xyz")
    with open(f"CALCS_defect_construction/defect2_{layer}layer.pkl", 'wb') as f:
        pickle.dump(defect, f)

    del defect
    gc.collect()


    ##########################################
    # defect3
    defect = copy.deepcopy(wz_lattice)
    for site in defect: 
        coord = site.real_space_coord
        if (coord[2] > z_cutoff) and (coord[1]>=-np.sqrt(3)*coord[0]): 
            site.has_atom = False
            site.iteration_changed = True

    update_whole_lattice_iteration(defect, sim_params)
    aggregate_to_xyz(defect, write_site=False, write_atoms=True, write_filename=f"CALCS_defect_construction/defect3_{layer}layer.xyz")
    with open(f"CALCS_defect_construction/defect3_{layer}layer.pkl", 'wb') as f:
        pickle.dump(defect, f)

    del defect
    gc.collect()



######################################################
if __name__ == "__main__":
    for l in range(1, 4): 
        create_defect(layer=l)