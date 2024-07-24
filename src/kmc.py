import time
import numpy as np
from multiprocessing import Pool, cpu_count
from .constants import *
from .sites import update_whole_lattice_iteration, check_neighbors, aggregate_to_xyz, update_XY_projection, write_projXY, update_XYvac

def kmc_step(site_list, sim_params, verbosity=0):
    """
    Perform a single kinetic Monte Carlo step. Currently doing a vanilla Gillespie algorithm implementation, since number of active sites (i.e., available events) is few. 

    Parameters:
    site_list (list of Site): The lattice containing all Site instances.
    sim_params (dict): Dictionary containing simulation parameters ('epsilon', 'mu', 'T').

    Returns:
    float or bool: The time increment if the step was successful, or False if no further events can occur.
    """

    # Initialize arrays for available events
    event_types = []
    event_siteID = []
    event_rates = []

    # Collect available events
    for site in site_list:
        if site.ready_to_attach and not site.ready_to_detach:
            event_types.append("a")                   # "a" for attachment
            event_siteID.append(site.wz_lattice_idx)  # Site index
            event_rates.append(site.atom_k_attach)    # Attachment rate
        elif site.ready_to_detach and not site.ready_to_attach:
            event_types.append("d")                   # "d" for detachment
            event_siteID.append(site.wz_lattice_idx)  # Site index
            event_rates.append(site.atom_k_detach)    # Detachment rate
        elif site.ready_to_attach and site.ready_to_detach:
            raise ValueError(f"Site {site.wz_lattice_idx} is marked as both ready to attach and detach, which is unphysical. ")

    event_rates = np.array(event_rates)
    if verbosity>0:
        print(f"\tNumber of available events: {len(event_rates)}")
    
    if len(event_rates) == 0:
        print("There are no available events. Thus, Gillespie algorithm can't proceed. ")
        return False

    # Calculate total rate and cumulative sum of rates
    total_rate = np.sum(event_rates)
    cumulative_sum = np.cumsum(event_rates)

    # Generate random numbers for event selection and time increment
    r1 = np.random.random()
    delta_t = -np.log(r1) / total_rate  # Time increment based on total rate
    
    r2 = np.random.random() * total_rate
    selected_event_idx = np.searchsorted(cumulative_sum, r2)  # Select event based on rates
    selected_event_type = event_types[selected_event_idx]
    selected_site_idx = event_siteID[selected_event_idx]
    selected_site = site_list[selected_site_idx]
    if verbosity>0:
        print(f"\tAt this iteration, we have selected to perform '{selected_event_type}' event on site {selected_site_idx}. ")

    # Perform the selected event (attachment or detachment)
    if selected_event_type == "a":
        selected_site.update_site_occupation_only(True, verbosity)
    elif selected_event_type == "d":
        selected_site.update_site_occupation_only(False, verbosity)
    selected_site.iteration_changed = True
    for neighbor_idx in selected_site.neighbor_sites_idx:
        site_list[neighbor_idx].iteration_changed = True

    # Update the entire lattice for the current iteration
    update_whole_lattice_iteration(site_list, sim_params)

    if verbosity>0:
        check_neighbors(site_list, just_atoms=True)

    return delta_t   # , total_rate (a float), cumulative_sum (a numpy array). 


def kmc_run(site_list, siteXY_list, vacXY_list, sim_params, trajectory_filename=None, write_every=20, runtime_flag=False):
    """
    Perform a KMC simulation and write the trajectory to a file compatible with VMD.

    Parameters:
    site_list (list of Site): The lattice containing all Site instances.
    sim_params (dict): Dictionary containing simulation parameters ('epsilon', 'mu', 'T', 'max_steps', 'max_time').
    trajectory_filename (str): The name of the output file for the trajectory.
    """

    # Initialize variables
    max_steps = sim_params.get('max_steps')
    max_time = sim_params.get('max_time')
    total_time = 0.0
    step_count = 0
    traj_In_coord = []
    traj_P_coord = []

    while True:
        start_time = time.time()
        delta_t = kmc_step(site_list, sim_params, calc_setting['verbosity'])
        end_time = time.time()
        if step_count % write_every == 0:
            print(f"\nKMC step {step_count}, sim time: {total_time:.5f}, this one KMC step elapsed time: {(end_time - start_time)*1000:.2f}ms") if runtime_flag else None
        if delta_t is False:
            print("No further events can occur. Simulation stopped.")
            break

        # Record the positions for this step
        if (trajectory_filename is not None) and (step_count % write_every == 0):
            start_time = time.time()
            frame_In_coord = []
            frame_P_coord = []
            for site in site_list:
                if (site.has_atom) and ((site.ready_to_detach) or any(site_list[neighborIdx].ready_to_detach for neighborIdx in site.neighbor_sites_idx)):
                    coord = site.real_space_coord
                    if site.cation_bool:
                        frame_In_coord.append(f"In {coord[0]} {coord[1]} {coord[2]}")
                    else:
                        frame_P_coord.append(f"P {coord[0]} {coord[1]} {coord[2]}")
            traj_In_coord.append(frame_In_coord)
            traj_P_coord.append(frame_P_coord)
            end_time = time.time()
            print(f"\tRecording frame for step {step_count}, elapsed time: {(end_time - start_time)*1000:.2f}ms") if runtime_flag else None
        
        # Accumulate stats
        if (step_count % write_every == 0):
            start_time = time.time()
            siteXY_list = update_XY_projection(site_list, siteXY_list)
            vacXY_list = update_XYvac(siteXY_list, vacXY_list)
            collect_stats(site_list, siteXY_list, vacXY_list, sim_params, writeProjXY_filePrefix=f"step_{step_count}")
            end_time = time.time()
            print(f"\tAccumulate stats, elapsed time: {(end_time - start_time)*1000:.2f}ms. ")

        if max_steps is not None and step_count >= max_steps:
            print(f"\nReached the maximum number of steps: {max_steps}. Simulation stopped.\n")
            break
        if max_time is not None and total_time >= max_time:
            print(f"\nReached the maximum simulation time: {max_time}. Simulation stopped.\n")
            break

        total_time += delta_t
        step_count += 1

    # Write the trajectory to the file
    max_In_num = max(len(inner_list) for inner_list in traj_In_coord)
    max_P_num = max(len(inner_list) for inner_list in traj_P_coord)
    if trajectory_filename is None:
        print("WARNING: The code is not writing out trajectory data. ")
    else:
        start_time = time.time()
        with open(trajectory_filename, 'w') as file:
            for iFrame in range(len(traj_In_coord)): 
                file.write(f"{(max_In_num + max_P_num)}\n")
                file.write("Frame\n")
                for iIn in range(max_In_num): 
                    if iIn < len(traj_In_coord[iFrame]):
                        file.write(f"{traj_In_coord[iFrame][iIn]}\n")
                    else: 
                        file.write(f"In 200.0 200.0 200.0\n")
                for iP in range(max_P_num): 
                    if iP < len(traj_P_coord[iFrame]):
                        file.write(f"{traj_P_coord[iFrame][iP]}\n")
                    else: 
                        file.write(f"P 200.0 200.0 200.0\n") 
        end_time = time.time()
        print(f"Trajectory written to '{trajectory_filename}', elapsed time: {(end_time - start_time):.2f}s")

    aggregate_to_xyz(site_list, write_site=False, write_atoms=True, write_filename=f"{sim_params['calc_dir']}final_atoms.xyz")
    siteXY_list = update_XY_projection(site_list, siteXY_list)
    vacXY_list = update_XYvac(siteXY_list, vacXY_list)
    collect_stats(site_list, siteXY_list, vacXY_list, sim_params, writeProjXY_filePrefix="final")
    results = {
        'total_steps': step_count,
        'total_time': total_time
    }
    return results


def collect_stats(site_list, siteXY_list, vacXY_list, sim_params, writeProjXY_filePrefix="init"): 
    num_cation = 0
    num_anion = 0
    numNeighborDouble = 0
    neighbor_3D_counts = [0, 0, 0, 0, 0]  # Store the num of sites with 0, 1, 2, 3, 4 neighbors
    neighbor_XY_counts = [0, 0, 0, 0]  # Store the num of sites with 0, 1, 2, 3 neighbors in the XY projection
    min_z = float('inf') 
    max_z = float('-inf')
    XY_a1_layer_counts = [0] * (max([siteXY.a1_label for siteXY in siteXY_list if siteXY.a1_label is not None])+1)   # should be 38 layers
    XY_a2_layer_counts = [0] * (max([siteXY.a2_label for siteXY in siteXY_list if siteXY.a2_label is not None])+1)
    XY_a3_layer_counts = [0] * (max([siteXY.a3_label for siteXY in siteXY_list if siteXY.a3_label is not None])+1)

    for site in site_list: 
        if site.has_atom and site.cation_bool: 
            num_cation += 1
        
        if site.has_atom and not site.cation_bool: 
            num_anion += 1

        if site.has_atom: 
            numNeighborDouble += site.num_atom_neighbors
            neighbor_3D_counts[site.num_atom_neighbors] += 1

            z = site.real_space_coord[2]
            if z < min_z:
                min_z = z
            if z > max_z:
                max_z = z

    for siteXY in siteXY_list: 
        if siteXY.has_atom: 
            neighbor_XY_counts[siteXY.num_atom_neighborsXY] += 1

            XY_a1_layer_counts[siteXY.a1_label] += 1
            XY_a2_layer_counts[siteXY.a2_label] += 1
            XY_a3_layer_counts[siteXY.a3_label] += 1

    energy = - sim_params['mu_In'] * num_cation - sim_params['mu_P'] * num_anion - sim_params['epsilon']*numNeighborDouble/2

    write_projXY(siteXY_list, writeXY_filename=f"{sim_params['calc_dir']}{writeProjXY_filePrefix}_projXY_atoms.xyz", mode='all_atoms')
    if (calc_setting['verbosity']>0) or (writeProjXY_filePrefix=="init") or (writeProjXY_filePrefix=="final"):
        write_projXY(siteXY_list, writeXY_filename=f"{sim_params['calc_dir']}{writeProjXY_filePrefix}_projXY_atoms_a1.xyz", mode='a1_label')
        write_projXY(siteXY_list, writeXY_filename=f"{sim_params['calc_dir']}{writeProjXY_filePrefix}_projXY_atoms_a2.xyz", mode='a2_label')
        write_projXY(siteXY_list, writeXY_filename=f"{sim_params['calc_dir']}{writeProjXY_filePrefix}_projXY_atoms_a3.xyz", mode='a3_label')

    # Vacancy stats to match experiments
    vacNeighbors_counts = [0, 0, 0, 0, 0, 0, 0]  # Store the num of sites with 0, 1, 2, 3, 4, 5, 6 neighbors
    for vacXY in vacXY_list: 
        if vacXY.light_up: 
            print(vacXY.num_litUp_neighborsVacXY)
            vacNeighbors_counts[vacXY.num_litUp_neighborsVacXY] += 1

    
    stats_dict = {
        "num_cation": num_cation, 
        "num_anion": num_anion, 
        "energy": energy, 
        "neighbor_3D_counts": neighbor_3D_counts, 
        "thickness": max_z - min_z, 
        "neighbor_XY_counts": neighbor_XY_counts, 
        "XY_a1_layer_counts": XY_a1_layer_counts, 
        "XY_a2_layer_counts": XY_a2_layer_counts, 
        "XY_a3_layer_counts": XY_a3_layer_counts, 
        "vacNeighbors_counts": vacNeighbors_counts, 
    }
    print(f"\t {stats_dict}")
    return stats_dict