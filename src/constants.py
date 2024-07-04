import numpy as np

# About the wurtzite lattice
wz_num_neighbors = 4

a = 4.17  # AA
c = 6.84  # AA
bond_length_max = 2.65  # AA, should be either 2.55 or 2.56 AA
unitCellVectors = np.array([[a, 0, 0], 
                            [-0.5 * a, np.sqrt(3) / 2 * a, 0], 
                            [0, 0, c]])
cationPos = np.array([[2 / 3, 1 / 3, -7 / 16], 
                        [1 / 3, 2 / 3, 1 / 16]]) @ unitCellVectors
anionPos = np.array([[2 / 3, 1 / 3, -1 / 16], 
                        [1 / 3, 2 / 3, 7 / 16]]) @ unitCellVectors

sites_pklFileRelPath = "../sites.pkl"

# Simulation parameters, default placeholder
sim_params = {
    'calc_dir': './', 
    'epsilon': 1.0, 
    'mu_In': -2.0,
    'mu_P': -2.0,
    'T': 300.0, 
    'max_time': None, 
    'max_steps': None
}

npl_params = {
    'read_sites_from': None, 
    'buffer': 20,
    'NPL_hex_diameter': 40,
    'NPL_thickness': 20
}

calc_setting = {
    'verbosity': 0, 
    'runtime_flag': 0, 
    'write_traj': 0, 
    'write_every': 50
}