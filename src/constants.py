import numpy as np

# About the wurtzite lattice
wz_num_neighbors = 4

a = 4.17  # AA
c = 6.84  # AA
bond_length_max = 2.65  # AA, should be either 2.55 or 2.56 AA
bond_length_max_XY = 2.60  # AA, should be 2.40755 AA. Must be smaller than 4.17000 AA. 
XY_neighbor_dist = 2.40755  # AA
veryFar = 300.0 
VacXY_neighbor_dist = 4.17000   # AA
VacXY_hex_pixel_area = 2*np.sqrt(3)*(VacXY_neighbor_dist/2)**2

unitCellVectors = np.array([[a, 0, 0], 
                            [-0.5 * a, np.sqrt(3) / 2 * a, 0], 
                            [0, 0, c]])
cationPos = np.array([[2 / 3, 1 / 3, -7 / 16], 
                        [1 / 3, 2 / 3, 1 / 16]]) @ unitCellVectors
anionPos = np.array([[2 / 3, 1 / 3, -1 / 16], 
                        [1 / 3, 2 / 3, 7 / 16]]) @ unitCellVectors

sites_pklFileRelPath = "../sites.pkl"

# Simulation parameters, default placeholder
calc_setting = {
    'calc_dir': './', 
    'verbosity': 0, 
    'runtime_flag': 0, 
    'write_every': 50, 
    'process_stats_now': 0,       # Default, store information and do post processing
    'random_seed': 32
}

npl_params = {
    'read_sites_from': None, 
    'buffer': 30,
    'NPL_hex_diameter': 160,
    'NPL_thickness': 40
}

sim_params = {
    'epsilon': 1.0, 
    'mu_In': -2.0,
    'mu_P': -2.0,
    'T': 300.0, 
    'max_time': None, 
    'max_steps': None
}