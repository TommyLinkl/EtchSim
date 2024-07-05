import time, os
import numpy as np
from typing import List
from .constants import *

class Site:
    def __init__(self, wz_lattice_idx, real_space_coord, cation_bool, neighbor_sites_idx, sim_params, has_atom=False, ready_to_attach=False, ready_to_detach=False):
        """
        Parameters:
            wz_lattice_idx (int): The index of the site in the list 'wz_lattice'.
            real_space_coord (array-like): The real space coordinates of the site.
            cation_bool (bool): True if the site is a cation site, False if it's an anion site.
            neighbor_sites_idx (array-like): Indices of neighboring sites in the list 'wz_lattice'.
        """
        self.wz_lattice_idx = wz_lattice_idx
        self.real_space_coord = np.array(real_space_coord)
        self.cation_bool = cation_bool
        self.neighbor_sites_idx = np.array(neighbor_sites_idx)
        
        self.neighbor_atoms_bool = np.zeros(len(neighbor_sites_idx), dtype=bool)
        self.num_atom_neighbors = np.sum(self.neighbor_atoms_bool)
        
        self.has_atom = has_atom
        self.ready_to_attach = ready_to_attach
        self.ready_to_detach = ready_to_detach
        self.update_ready_to_attach()
        self.update_ready_to_detach()
        
        self.iteration_changed = False
        
        self.atom_k_attach = None
        self.atom_k_detach = None
        self.update_rates(sim_params)


    def update_ready_to_attach(self): 
        if (not self.has_atom) and (self.num_atom_neighbors>0): 
            self.ready_to_attach = True
        else: 
            self.ready_to_attach = False


    def update_ready_to_detach(self): 
        if (self.has_atom) and (self.num_atom_neighbors<wz_num_neighbors): 
            self.ready_to_detach = True
        else: 
            self.ready_to_detach = False


    def update_site_occupation_only(self, occ_bool, verbosity=0): 
        if occ_bool == self.has_atom:
            pass
        elif occ_bool:
            print(f'\tSite {self.wz_lattice_idx} is getting occupied by an atom.') if verbosity>0 else None
            self.has_atom = occ_bool
            self.iteration_changed = True
        else: 
            print(f'Atom at site {self.wz_lattice_idx} is detaching.') if verbosity>0 else None
            self.has_atom = occ_bool
            self.iteration_changed = True


    def update_neighbor_occupation_only(self, neighbor_wz_lattice_idx, neighbor_occ_bool):
        """
        Parameters:
        neighbor_wz_lattice_idx (int): The index of the neighbor site (that has been updated) in the list 'wz_lattice'.
        neighbor_occ_bool (bool): The new occupation status of this neighbor site.
        """
        if neighbor_wz_lattice_idx not in self.neighbor_sites_idx: 
            raise ValueError(f'Site {neighbor_wz_lattice_idx} is not a neighbor of site {self.wz_lattice_idx}.')
        else: 
            neighbor_idx = np.where(self.neighbor_sites_idx == neighbor_wz_lattice_idx)[0][0]
            self.neighbor_atoms_bool[neighbor_idx] = neighbor_occ_bool
            self.num_atom_neighbors = np.sum(self.neighbor_atoms_bool)
            self.iteration_changed = True


    def update_rates(self, sim_params):
        epsilon = sim_params['epsilon']
        if self.cation_bool:
            mu = sim_params['mu_In']
        else: 
            mu = sim_params['mu_P']
        T = sim_params['T']
        self.atom_k_attach = np.exp( (epsilon * self.num_atom_neighbors + mu) / (2 * T) ) 
        self.atom_k_detach = np.exp(-(epsilon * self.num_atom_neighbors + mu) / (2 * T) ) 


class SiteXY:
    def __init__(self, siteXY_id, site3D_ids, coordXY: np.ndarray, has_atom, neighborXY_sites_idx):
        self.siteXY_id = siteXY_id
        self.site3D_ids = site3D_ids
        self.coordXY = coordXY
        self.has_atom = has_atom

        self.neighborXY_sites_idx = np.array(neighborXY_sites_idx)
        self.neighborXY_atoms_bool = np.zeros(len(neighborXY_sites_idx), dtype=bool)
        self.num_atom_neighborsXY = np.sum(self.neighborXY_atoms_bool)

        self.measure_a1_a2_a3()


    def measure_a1_a2_a3(self): 
        x = self.coordXY[0]
        y = self.coordXY[1]
        self.a1_measure = y
        self.a2_measure = y - np.sqrt(3)*x
        self.a3_measure = y + np.sqrt(3)*x


def initialize_wz_lattice(max_xy, max_z, sim_params, verbosity=0):
    """
    Initialize the large wurtzite lattice and all sites are initialized as empty. This function is specific to this usage. 

    Parameters:
    max_xy (float): Maximum xy-coordinate extent.
    max_z (float): Maximum z-coordinate extent.

    Returns:
    list of Site: The initialized lattice containing all Site instances.
    """
    wz_lattice = []

    max_lattice_numbers_x = int(np.ceil(max_xy/2 / a) * 1.7)
    max_lattice_numbers_y = int(np.ceil(max_xy/2 / (a * np.sqrt(3) / 2))* 1.7) 
    max_lattice_numbers_z = int(np.ceil(max_z/2 / c)* 1.7)

    # Generate lattice positions within the specified boundaries of a "square prism" 
    start_time = time.time()
    for i in np.arange(-max_lattice_numbers_x, max_lattice_numbers_x, 1):
        for j in np.arange(-max_lattice_numbers_y, max_lattice_numbers_y, 1):
            for k in np.arange(-max_lattice_numbers_z, max_lattice_numbers_z, 1):
                lattice_vector = i * unitCellVectors[0] + j * unitCellVectors[1] + k * unitCellVectors[2]

                if (np.abs(lattice_vector[0]) <= max_xy/2) and (np.abs(lattice_vector[1]) <= max_xy/2) and (np.abs(lattice_vector[2]) <= max_z/2): 
                    for pos in cationPos:
                        site_coord = lattice_vector + pos
                        wz_lattice.append(Site(len(wz_lattice), site_coord, True, [], sim_params, has_atom=False))
                    for pos in anionPos:
                        site_coord = lattice_vector + pos
                        wz_lattice.append(Site(len(wz_lattice), site_coord, False, [], sim_params, has_atom=False))
    end_time = time.time()
    print(f"Generate large lattice, elapsed time: {(end_time - start_time):.2f}s")

    # Partition the lattice into a 3D grid
    start_time = time.time()
    grid_size = 2.5 * bond_length_max   # Algorithmically, 2 * bond_length_max is enough. But I just want to be safe. 
    grid = {}
    for site in wz_lattice:
        grid_coord = tuple((site.real_space_coord // grid_size).astype(int))
        if grid_coord not in grid:
            grid[grid_coord] = []
        grid[grid_coord].append(site)

    # Find neighbors for each site
    for site in wz_lattice:
        grid_coord = tuple((site.real_space_coord // grid_size).astype(int))
        neighbors = []
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                for dz in [-1, 0, 1]:
                    neighbor_grid_coord = (grid_coord[0] + dx, grid_coord[1] + dy, grid_coord[2] + dz)
                    if neighbor_grid_coord in grid:
                        neighbors.extend(grid[neighbor_grid_coord])
        site.neighbor_sites_idx = np.array([neighbor.wz_lattice_idx for neighbor in neighbors if neighbor != site and np.linalg.norm(site.real_space_coord - neighbor.real_space_coord) <= bond_length_max])
    end_time = time.time()
    print(f"Neighbor information, elapsed time: {(end_time - start_time):.2f}s")

    if verbosity>0: 
        check_neighbors(wz_lattice) 
    aggregate_to_xyz(wz_lattice, write_site=True, write_atoms=False, write_filename=f"{sim_params['calc_dir']}init_lattice.xyz")

    return wz_lattice
    

def initialize_hex_NPL(site_list, NPL_hex_diameter, NPL_thickness, sim_params, verbosity=0):
    """
    Initialize the atoms within the NPL boundary. This function is specific to this usage.

    Parameters:
    site_list (list of Site): The initialized lattice containing all Site instances.
    NPL_hex_diamter (float): Hexagonal diameter for the NPL.
    NPL_thickness (float): Thickness of the NPL.
    """
    start_time = time.time()
    for site in site_list:
        if NPL_geom_boundary(site.real_space_coord, NPL_hex_diameter, NPL_thickness):
            site.has_atom = True
            site.iteration_changed = True

            # Change its neighbors' iteration_changed to True
            for neighbor_idx in site.neighbor_sites_idx:
                site_list[neighbor_idx].iteration_changed = True

    update_whole_lattice_iteration(site_list, sim_params)

    end_time = time.time()
    print(f"Initialize NPL atoms, elapsed time: {(end_time - start_time):.2f}s")

    if verbosity>0:
        check_neighbors(site_list, just_atoms=True)

    aggregate_to_xyz(site_list, write_site=True, write_atoms=True, write_filename=f"{sim_params['calc_dir']}init_lattice_atoms.xyz")
    aggregate_to_xyz(site_list, write_site=False, write_atoms=True, write_filename=f"{sim_params['calc_dir']}init_atoms.xyz")

    return


def update_whole_lattice_iteration(site_list, sim_params):
    """
    Update the entire lattice for the current iteration.
    
    Parameters:
    site_list (list of Site): The list of Site instances.
    """
    for site in site_list:
        if site.iteration_changed:
            site.neighbor_atoms_bool = np.array([site_list[neighbor_idx].has_atom for neighbor_idx in site.neighbor_sites_idx])
            site.num_atom_neighbors = np.sum(site.neighbor_atoms_bool)
            
            site.update_ready_to_attach()
            site.update_ready_to_detach()
            site.update_rates(sim_params)
            
            # Reset the iteration_changed status
            site.iteration_changed = False


def aggregate_to_xyz(site_list, write_site=True, write_atoms=True, write_filename='init_lattice_atoms.xyz'):
    """
    Write the site and atom information to file.
    
    Parameters:
    site_list (list of Site): List of Site instances.
    """

    if write_site and write_atoms: 
        totalNum = len(site_list) + sum(1 for site in site_list if site.has_atom)
        comment = "# Both lattice sites and atoms are written. Sites are represented as In and P, while atoms are represented using Tl and As for larger atoms."
    elif write_site and not write_atoms: 
        totalNum = len(site_list)
        comment = "# Only lattice sites are written. Sites are represented as In and P."
    elif not write_site and write_atoms: 
        totalNum = sum(1 for site in site_list if site.has_atom)
        comment = "# Only atoms are written. Atoms are represented as In and P."
    else: 
        print("You are not writing anything to file. ")
        return

    with open(write_filename, "w") as file:
        file.write(f"{totalNum}\n")
        file.write(f"{comment}\n")

        if write_site and write_atoms: 
            for site in site_list:
                element = "In" if site.cation_bool else "P"
                coord = site.real_space_coord
                file.write(f"{element}   {coord[0]:.6f}   {coord[1]:.6f}   {coord[2]:.6f}\n")
            for site in site_list:
                if site.has_atom: 
                    element = "Tl" if site.cation_bool else "As"
                    coord = site.real_space_coord
                    file.write(f"{element}   {coord[0]:.6f}   {coord[1]:.6f}   {coord[2]:.6f}\n")
            
        elif write_site and not write_atoms: 
            for site in site_list:
                element = "In" if site.cation_bool else "P"
                coord = site.real_space_coord
                file.write(f"{element}   {coord[0]:.6f}   {coord[1]:.6f}   {coord[2]:.6f}\n")

        elif not write_site and write_atoms: 
            for site in site_list:
                if site.has_atom: 
                    element = "In" if site.cation_bool else "P"
                    coord = site.real_space_coord
                    file.write(f"{element}   {coord[0]:.6f}   {coord[1]:.6f}   {coord[2]:.6f}\n")

    return


def check_neighbors(site_list, just_atoms=False):
    """
    Check the neighbor relationships in the wurtzite lattice or only among atoms.

    Parameters:
    site_list (list of Site): The list of Site instances representing the lattice.
    just_atoms (bool, optional): If True, only check atom neighbors; otherwise, check lattice site neighbors as well. Default is False.
    """
    if not just_atoms:
        print("Checking site neighbor relationships: ")

        neighbor_counts = [0, 0, 0, 0, 0, 0]  # To store the count of sites with 0, 1, 2, 3, 4, and >4 neighbors
        
        for site in site_list:
            num_neighbors = len(site.neighbor_sites_idx)
            if num_neighbors > 4:
                print(f"\tWARNING: Site {site.wz_lattice_idx} has more than 4 neighbors: {num_neighbors}")
                neighbor_counts[5] += 1
            else:
                neighbor_counts[num_neighbors] += 1
        
        for i in range(5):
            print(f"\tNumber of sites with {i} neighbors: {neighbor_counts[i]}")
        print(f"\tNumber of sites with 4+ neighbors: {neighbor_counts[5]}")
        
        # Check if neighbor list is reciprocal
        reciprocal_errors = 0
        for site in site_list:
            for neighbor_idx in site.neighbor_sites_idx:
                neighbor = site_list[neighbor_idx]
                if site.wz_lattice_idx not in neighbor.neighbor_sites_idx:
                    print(f"\tReciprocal error: Site {site.wz_lattice_idx} has neighbor {neighbor_idx}, but the reverse is not true.")
                    reciprocal_errors += 1
        print(f"\tTotal reciprocal errors: {reciprocal_errors}")

        # Check that cation neighbors are all anions and vice versa
        neighbor_type_errors = 0
        for site in site_list:
            for neighbor_idx in site.neighbor_sites_idx:
                neighbor = site_list[neighbor_idx]
                if site.cation_bool == neighbor.cation_bool:
                    print(f"\tNeighbor type error: Site {site.wz_lattice_idx} has the same type (cation-vs-anion) neighbor at site {neighbor_idx}.")
                    neighbor_type_errors += 1
        print(f"\tTotal neighbor type errors: {neighbor_type_errors}")

    elif just_atoms: 
        print(f"Checking only atom neighbors. Ignoring lattice structure.")

        neighbor_counts = [0, 0, 0, 0, 0, 0]
        total_atoms = 0
    
        for site in site_list:
            if site.has_atom:
                total_atoms += 1
                num_neighbors = np.sum(site.neighbor_atoms_bool)
                if num_neighbors > 4:
                    print(f"\tWARNING: Atom at site {site.wz_lattice_idx} has more than 4 neighbors: {num_neighbors}")
                    neighbor_counts[5] += 1
                else:
                    neighbor_counts[num_neighbors] += 1
        
        print(f"\tThere are {total_atoms} atoms in total. ")
        for i in range(5):
            print(f"\tNumber of atoms with {i} atom neighbors: {neighbor_counts[i]}")
        print(f"\tNumber of atoms with 4+ atom neighbors: {neighbor_counts[5]}")
    
    # Check ready_to_attach and ready_to_detach: 
    event_errors = 0
    for site in site_list: 
        if site.ready_to_attach and site.ready_to_detach:
            print(f"Site {site.wz_lattice_idx} is marked as both ready to attach and detach, which is unphysical.")
            event_errors += 1
    print(f"\tTotal event (attach-vs-detach) errors: {event_errors}")


def NPL_geom_boundary(site_coord, hex_diameter, thickness, NPL_center=np.array([0.0, 0.0, 0.0])):
    """
    Determine if a site is within the boundary of a hexagonal nanoplatelet (NPL).

    Parameters:
    site_coord (np.array): The coordinates of the site to check.
    hex_diameter (float): The diameter of the circle that circumscribes the hexagon. 
    thickness (float): The thickness of the NPL.
    NPL_center (np.array): The center of the NPL (default is the origin).
    """
    part_of_NPL = False
    
    # Translate site coordinates to NPL-centered coordinates
    site_coord = site_coord - NPL_center
    x, y, z = site_coord

    # Please the following link for details: https://www.desmos.com/calculator/orepbmxo3q
    part_of_hex = (np.abs(x) <= hex_diameter/2) and (np.abs(y) <= hex_diameter/2 * np.sqrt(3) / 2) and (np.abs(x) + np.abs(y) / np.sqrt(3) <= hex_diameter/2)

    part_of_thickness = np.abs(z) <= thickness / 2
    
    if part_of_hex and part_of_thickness:
        part_of_NPL = True

    return part_of_NPL


def project_to_XY(site_list):  # Designed to be called only once
    coord_dict = {}

    for site in site_list:
        x, y = site.real_space_coord[:2]
        coord = (x, y)
        if coord not in coord_dict:
            coord_dict[coord] = {"coordXY": np.array([x, y]), "has_atom": False, "site3D_ids": []}
        if site.has_atom:
            coord_dict[coord]["has_atom"] = True
        coord_dict[coord]["site3D_ids"].append(site.wz_lattice_idx)

    # Convert the dictionary to a list of SiteXY objects
    siteXY_list = []
    for i, (coord, info) in enumerate(coord_dict.items()):
        thisSiteXY = SiteXY(siteXY_id=i, site3D_ids=info["site3D_ids"], coordXY=info["coordXY"], has_atom=info["has_atom"], neighborXY_sites_idx=[])
        siteXY_list.append(thisSiteXY)

    # Assign layers in a1, a2, a3 directions
    group_XY_layers(siteXY_list)

    # Partition the lattice into a 2D grid
    grid_size = 2.5 * bond_length_max_XY   # Algorithmically, 2 * bond_length_max is enough. But I just want to be safe. 
    grid = {}
    for siteXY in siteXY_list:
        grid_coord = tuple((siteXY.coordXY // grid_size).astype(int))
        if grid_coord not in grid:
            grid[grid_coord] = []
        grid[grid_coord].append(siteXY)

    for siteXY in siteXY_list:
        grid_coord = tuple((siteXY.coordXY // grid_size).astype(int))
        neighbors = []
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                neighbor_grid_coord = (grid_coord[0] + dx, grid_coord[1] + dy)
                if neighbor_grid_coord in grid:
                    neighbors.extend(grid[neighbor_grid_coord])
        siteXY.neighborXY_sites_idx = np.array([neighbor.siteXY_id for neighbor in neighbors if neighbor != siteXY and np.linalg.norm(siteXY.coordXY - neighbor.coordXY) <= bond_length_max_XY])
    
    update_XY_projection(site_list, siteXY_list)

    return siteXY_list


def update_XY_projection(site_list, siteXY_list): 
    for siteXY in siteXY_list: 
        siteXY.has_atom = False
        for site3DIdx in siteXY.site3D_ids:
            if site_list[site3DIdx].has_atom: 
                siteXY.has_atom = True
                break
    
    for siteXY in siteXY_list: 
        siteXY.neighborXY_atoms_bool = np.array([siteXY_list[neighbor_idx].has_atom for neighbor_idx in siteXY.neighborXY_sites_idx])
        siteXY.num_atom_neighborsXY = np.sum(siteXY.neighborXY_atoms_bool)
    
    return siteXY_list
    

def write_projXY(siteXY_list, writeXY_filename='projXY_atoms.xyz', mode='all_atoms'):
    fake_atom_species = {
        0: 'Be', 1: 'He', 2: 'Al', 3: 'Ne', 4: 'Na'
    }

    nAtomsXY = sum(1 for siteXY in siteXY_list if siteXY.has_atom)
    comment = "# Only atoms are written. Atoms are represented as H."

    with open(writeXY_filename, "w") as file:
        file.write(f"{nAtomsXY}\n")
        file.write(f"{comment}\n")
        for siteXY in siteXY_list:
            if siteXY.has_atom: 
                # file.write(f"H   {siteXY.coordXY[0]:.6f}   {siteXY.coordXY[1]:.6f}   0.0   {siteXY.a1_label}\n")
                if mode=='a1_label':
                    file.write(f"{fake_atom_species[(siteXY.a1_label)%5]}   {siteXY.coordXY[0]:.6f}   {siteXY.coordXY[1]:.6f}   0.0\n")
                elif mode=='a2_label':
                    file.write(f"{fake_atom_species[(siteXY.a2_label)%5]}   {siteXY.coordXY[0]:.6f}   {siteXY.coordXY[1]:.6f}   0.0\n")
                elif mode=='a3_label':
                    file.write(f"{fake_atom_species[(siteXY.a3_label)%5]}   {siteXY.coordXY[0]:.6f}   {siteXY.coordXY[1]:.6f}   0.0\n")
                elif mode=='all_atoms':
                    file.write(f"Na   {siteXY.coordXY[0]:.6f}   {siteXY.coordXY[1]:.6f}   0.0\n")
                else: 
                    raise NotImplementedError("We currently only allow mode to be 'a1_label', 'a2_label', 'a3_label', or 'all_atoms'. ")
    return


def group_XY_layers(siteXY_list, a=XY_neighbor_dist):  # Designed to be called only once, at init
    a1_measure_list = [(siteXY.a1_measure, siteXY.siteXY_id) for siteXY in siteXY_list if siteXY.has_atom]
    sorted_a1_measure_list = sorted(a1_measure_list, key=lambda x: x[0])

    current_label = 0
    sorted_a1_measure_list[0] = (sorted_a1_measure_list[0][0], sorted_a1_measure_list[0][1], current_label)
    for i in range(1, len(sorted_a1_measure_list)):
        current_measure = sorted_a1_measure_list[i][0]
        previous_measure = sorted_a1_measure_list[i-1][0]

        if np.isclose(current_measure - previous_measure, a, atol=1e-4):
            current_label += 1

        sorted_a1_measure_list[i] = (current_measure, sorted_a1_measure_list[i][1], current_label)

    for siteXY in siteXY_list:
        for measure, siteXY_id, label in sorted_a1_measure_list:
            if siteXY.siteXY_id == siteXY_id:
                siteXY.a1_label = label

    # Repeat for a2
    a2_measure_list = [(siteXY.a2_measure, siteXY.siteXY_id) for siteXY in siteXY_list if siteXY.has_atom]
    sorted_a2_measure_list = sorted(a2_measure_list, key=lambda x: x[0])

    current_label = 0
    sorted_a2_measure_list[0] = (sorted_a2_measure_list[0][0], sorted_a2_measure_list[0][1], current_label)
    for i in range(1, len(sorted_a2_measure_list)):
        current_measure = sorted_a2_measure_list[i][0]
        previous_measure = sorted_a2_measure_list[i-1][0]

        if np.isclose(current_measure - previous_measure, 2*a, atol=1e-4):
            current_label += 1

        sorted_a2_measure_list[i] = (current_measure, sorted_a2_measure_list[i][1], current_label)

    for siteXY in siteXY_list:
        for measure, siteXY_id, label in sorted_a2_measure_list:
            if siteXY.siteXY_id == siteXY_id:
                siteXY.a2_label = label

    # Repeat for a3
    a3_measure_list = [(siteXY.a3_measure, siteXY.siteXY_id) for siteXY in siteXY_list if siteXY.has_atom]
    sorted_a3_measure_list = sorted(a3_measure_list, key=lambda x: x[0])

    current_label = 0
    sorted_a3_measure_list[0] = (sorted_a3_measure_list[0][0], sorted_a3_measure_list[0][1], current_label)
    for i in range(1, len(sorted_a3_measure_list)):
        current_measure = sorted_a3_measure_list[i][0]
        previous_measure = sorted_a3_measure_list[i-1][0]

        if np.isclose(current_measure - previous_measure, 2*a, atol=1e-4):
            current_label += 1

        sorted_a3_measure_list[i] = (current_measure, sorted_a3_measure_list[i][1], current_label)

    for siteXY in siteXY_list:
        for measure, siteXY_id, label in sorted_a3_measure_list:
            if siteXY.siteXY_id == siteXY_id:
                siteXY.a3_label = label

    return siteXY_list


