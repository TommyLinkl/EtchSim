import time
import numpy as np
import json, csv, gzip
from multiprocessing import Pool, cpu_count
from scipy.spatial import ConvexHull, QhullError
from .constants import *
from .sites import update_whole_lattice_iteration, check_neighbors, aggregate_to_xyz, update_XY_projection, write_projXY, update_XYvac, write_XY_sites_vac

def kmc_step(site_list, sim_params, verbosity=0):
    """
    Perform a single kinetic Monte Carlo step. Currently doing a vanilla Gillespie algorithm implementation, since number of active sites (i.e., available events) is few. 

    Parameters:
    site_list (list of Site): The lattice containing all Site instances.
    sim_params (dict): Dictionary containing simulation parameters ('epsilon', 'mu', 'T').

    Returns:
    float or bool: The time increment if the step was successful, or False if no further events can occur.
    """

    # Collect available events
    start_time = time.time()

    event_types = []
    event_siteID = []
    event_rates = []

    for site in site_list:
        if not site.ready_to_attach and not site.ready_to_detach:
            continue

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
    end_time = time.time()
    print(f"\tCollect available events: {(end_time - start_time)*1000:.2f}ms") if verbosity > 1 else None

    # Calculate total rate and cumulative sum of rates
    start_time = time.time()
    total_rate = np.sum(event_rates)
    cumulative_sum = np.cumsum(event_rates)
    end_time = time.time()
    print(f"\tCalculate total rate and cumulative sum of rates: {(end_time - start_time)*1000:.2f}ms") if verbosity > 1 else None

    # Generate random numbers for event selection and time increment
    start_time = time.time()
    r1 = np.random.random()
    delta_t = -np.log(r1) / total_rate  # Time increment based on total rate
    
    r2 = np.random.random() * total_rate
    selected_event_idx = np.searchsorted(cumulative_sum, r2)  # Select event based on rates
    selected_event_type = event_types[selected_event_idx]
    selected_site_idx = event_siteID[selected_event_idx]
    selected_site = site_list[selected_site_idx]
    if verbosity>0:
        print(f"\tAt this iteration, we have selected to perform '{selected_event_type}' event on site {selected_site_idx}. ")
    end_time = time.time()
    print(f"\tChoose event: {(end_time - start_time)*1000:.2f}ms") if verbosity > 1 else None

    # Perform the selected event (attachment or detachment)
    start_time = time.time()
    if selected_event_type == "a":
        selected_site.update_site_occupation_only(True, verbosity)
    elif selected_event_type == "d":
        selected_site.update_site_occupation_only(False, verbosity)
    selected_site.iteration_changed = True
    for neighbor_idx in selected_site.neighbor_sites_idx:
        site_list[neighbor_idx].iteration_changed = True
    end_time = time.time()
    print(f"\tPerform event: {(end_time - start_time)*1000:.2f}ms") if verbosity > 1 else None

    # Update the entire lattice for the current iteration
    start_time = time.time()
    update_whole_lattice_iteration(site_list, sim_params)
    end_time = time.time()
    print(f"\tUpdate the entire lattice for the current iteration: {(end_time - start_time)*1000:.2f}ms\n") if verbosity > 1 else None

    if verbosity>0:
        check_neighbors(site_list, just_atoms=True)

    return delta_t   # , total_rate (a float), cumulative_sum (a numpy array). 


def kmc_run(site_list, siteXY_list, vacXY_list, sim_params, write_every=20, runtime_flag=False):
    """
    Perform a KMC simulation and write the trajectory to a file compatible with VMD.

    Parameters:
    site_list (list of Site): The lattice containing all Site instances.
    sim_params (dict): Dictionary containing simulation parameters ('epsilon', 'mu', 'T', 'max_steps', 'max_time').
    """
    ################################# initialize #################################
    # Initialize variables
    max_steps = sim_params.get('max_steps')
    max_time = sim_params.get('max_time')
    total_time = 0.0
    step_count = 0

    trajFileName_zip = f"{calc_setting['calc_dir']}traj.xyz.gz"
    # trajXYFileName_zip = f"{calc_setting['calc_dir']}trajXY.xyz.gz"
    trajVacFileName_zip = f"{calc_setting['calc_dir']}trajVac.xyz.gz"
    with open(trajFileName_zip, 'w') as f:
        pass
    # with open(trajXYFileName_zip, 'w') as f:
    #     pass
    with open(trajVacFileName_zip, 'w') as f:
        pass
    traj_In_coord = []
    traj_P_coord = []
    stats_list = []

    # Store information for post-processing
    forLaterStatsFileName = f"{calc_setting['calc_dir']}forStatsLater.csv.gz"
    with gzip.open(forLaterStatsFileName, 'wt', newline='') as file:
        writer = csv.writer(file)
        header = ['stepNum', 'time'] + [f'site_{i}' for i in range(len(site_list))]
        writer.writerow(header)

    ################################# KMC step 0 #################################
    update_whole_lattice_iteration(site_list, sim_params)
    if calc_setting['process_stats_now']!=0:
        # Accumulate stats
        start_time = time.time()

        siteXY_list = update_XY_projection(site_list, siteXY_list)
        vacXY_list = update_XYvac(siteXY_list, vacXY_list)

        stats = collect_stats(site_list, siteXY_list, vacXY_list, sim_params, writeProjXY_filePrefix=f"step_{step_count}")
        stats['stepNum'] = step_count
        stats['simTime'] = total_time
        stats_list.append(stats)

        end_time = time.time()
        print(f"\tAccumulate stats, elapsed time: {(end_time - start_time)*1000:.2f}ms. ")

    # Store data (3D positions) for post-processing
    row = [step_count, total_time] + [int(site.has_atom) for site in site_list]
    with gzip.open(forLaterStatsFileName, 'at', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(row)

    # Record the positions of this step for plottable trajectory
    start_time = time.time()
    # Record 3D positions 
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

    if calc_setting['process_stats_now']!=0:
        # Record vacXY_trajectory
        with gzip.open(trajVacFileName_zip, 'at') as file:
            file.write(f"{len(siteXY_list) + len(vacXY_list)}\n")
            file.write("Frame\n")
            for siteXY in siteXY_list: 
                if siteXY.has_atom: 
                    file.write(f"H {siteXY.coordXY[0]} {siteXY.coordXY[1]} 0.0\n")
                else: 
                    file.write(f"H {veryFar} {veryFar} 0.0\n")
            for vacXY in vacXY_list: 
                if vacXY.light_up: 
                    file.write(f"Be {vacXY.coordXY[0]} {vacXY.coordXY[1]} 0.0\n")
                else: 
                    file.write(f"Be {veryFar} {veryFar} 0.0\n")

    end_time = time.time()
    print(f"\tRecording frame for step {step_count}, elapsed time: {(end_time - start_time)*1000:.2f}ms") if runtime_flag else None

    ################################# KMC step 1+ #################################
    while True: 
        # Make a kmc_step move. 
        start_time = time.time()
        delta_t = kmc_step(site_list, sim_params, calc_setting['verbosity'])
        end_time = time.time()
        if step_count % write_every == 0:
            print(f"\nKMC step {step_count}, sim time: {total_time:.5f}, this one KMC step elapsed time: {(end_time - start_time)*1000:.2f}ms") if runtime_flag else None
        if delta_t is False:
            print("No further events can occur. Simulation stopped.")
            break
        total_time += delta_t
        step_count += 1

        if (step_count % write_every == 0):
            if calc_setting['process_stats_now']!=0:
                # Accumulate stats
                start_time = time.time()

                siteXY_list = update_XY_projection(site_list, siteXY_list)
                vacXY_list = update_XYvac(siteXY_list, vacXY_list)

                stats = collect_stats(site_list, siteXY_list, vacXY_list, sim_params, writeProjXY_filePrefix=f"step_{step_count}")
                stats['stepNum'] = step_count
                stats['simTime'] = total_time
                stats_list.append(stats)

                end_time = time.time()
                print(f"\tAccumulate stats, elapsed time: {(end_time - start_time)*1000:.2f}ms. ")

            # Store data (3D positions) for post-processing
            row = [step_count, total_time] + [int(site.has_atom) for site in site_list]
            with gzip.open(forLaterStatsFileName, 'at', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(row)

            # Record the positions of this step for plottable trajectory
            start_time = time.time()

            # Record 3D positions 
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

            if calc_setting['process_stats_now']!=0:
                ''' # Now redundant
                # Record projected XY_trajectory
                with gzip.open(trajXYFileName_zip, 'at') as file:
                    file.write(f"{len(siteXY_list)}\n")
                    file.write("Frame\n")
                    for siteXY in siteXY_list: 
                        if siteXY.has_atom: 
                            file.write(f"Na {siteXY.coordXY[0]} {siteXY.coordXY[1]} 0.0\n")
                        else: 
                            file.write(f"Na {veryFar} {veryFar} 0.0\n")
                '''

                # Record vacXY_trajectory
                with gzip.open(trajVacFileName_zip, 'at') as file:
                    file.write(f"{len(siteXY_list) + len(vacXY_list)}\n")
                    file.write("Frame\n")
                    for siteXY in siteXY_list: 
                        if siteXY.has_atom: 
                            file.write(f"H {siteXY.coordXY[0]} {siteXY.coordXY[1]} 0.0\n")
                        else: 
                            file.write(f"H {veryFar} {veryFar} 0.0\n")
                    for vacXY in vacXY_list: 
                        if vacXY.light_up: 
                            file.write(f"Be {vacXY.coordXY[0]} {vacXY.coordXY[1]} 0.0\n")
                        else: 
                            file.write(f"Be {veryFar} {veryFar} 0.0\n")

            end_time = time.time()
            print(f"\tRecording frame for step {step_count}, elapsed time: {(end_time - start_time)*1000:.2f}ms") if runtime_flag else None

        if max_steps is not None and step_count >= max_steps:
            print(f"\nReached the maximum number of steps: {max_steps}. Simulation stopped.\n")
            break
        if max_time is not None and total_time >= max_time:
            print(f"\nReached the maximum simulation time: {max_time}. Simulation stopped.\n")
            break

    # Write the trajectory to the file
    max_In_num = max(len(inner_list) for inner_list in traj_In_coord)
    max_P_num = max(len(inner_list) for inner_list in traj_P_coord)
    start_time = time.time()
    with gzip.open(trajFileName_zip, 'wt') as file:
        for iFrame in range(len(traj_In_coord)): 
            file.write(f"{(max_In_num + max_P_num)}\n")
            file.write("Frame\n")
            for iIn in range(max_In_num): 
                if iIn < len(traj_In_coord[iFrame]):
                    file.write(f"{traj_In_coord[iFrame][iIn]}\n")
                else: 
                    file.write(f"In {veryFar} {veryFar} {veryFar}\n")
            for iP in range(max_P_num): 
                if iP < len(traj_P_coord[iFrame]):
                    file.write(f"{traj_P_coord[iFrame][iP]}\n")
                else: 
                    file.write(f"P {veryFar} {veryFar} {veryFar}\n") 
    end_time = time.time()
    print(f"Trajectory (3D) written to '{trajFileName_zip}', elapsed time: {(end_time - start_time):.2f}s")

    aggregate_to_xyz(site_list, write_site=False, write_atoms=True, write_filename=f"{calc_setting['calc_dir']}final_atoms.xyz") if calc_setting['verbosity']>0 else None
    if calc_setting['process_stats_now']!=0:
        siteXY_list = update_XY_projection(site_list, siteXY_list)
        vacXY_list = update_XYvac(siteXY_list, vacXY_list)
        stats = collect_stats(site_list, siteXY_list, vacXY_list, sim_params, writeProjXY_filePrefix="final")
        stats['stepNum'] = step_count
        stats['simTime'] = total_time
        stats_list.append(stats)

    if calc_setting['process_stats_now']!=0:
        # Dump the stats in a json file
        with open(f"{calc_setting['calc_dir']}stats.json", 'w') as file:
            # json.dump(stats_list, file, indent=4)
            json.dump(stats_list, file, separators=(',', ':'))


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
    min_z = float('inf') 
    max_z = float('-inf')
    neighbor_XY_counts = [0, 0, 0, 0]  # Store the num of sites with 0, 1, 2, 3 neighbors in the XY projection
    XY_a1_layer_counts = [0] * (max([siteXY.a1_label for siteXY in siteXY_list if siteXY.a1_label is not None])+1)   # should be 38 layers
    XY_a2_layer_counts = [0] * (max([siteXY.a2_label for siteXY in siteXY_list if siteXY.a2_label is not None])+1)
    XY_a3_layer_counts = [0] * (max([siteXY.a3_label for siteXY in siteXY_list if siteXY.a3_label is not None])+1)
    vacNeighbors_counts = [0, 0, 0, 0, 0, 0, 0]  # Store the num of sites with 0, 1, 2, 3, 4, 5, 6 neighbors
    vac_a1_layer_counts = [0] * (max([vacXY.a1_label for vacXY in vacXY_list if vacXY.a1_label is not None])+1)
    vac_a2_layer_counts = [0] * (max([vacXY.a2_label for vacXY in vacXY_list if vacXY.a2_label is not None])+1)
    vac_a3_layer_counts = [0] * (max([vacXY.a3_label for vacXY in vacXY_list if vacXY.a3_label is not None])+1)

    ###############################################
    update_whole_lattice_iteration(site_list, sim_params)
    update_XY_projection(site_list, siteXY_list)
    update_XYvac(siteXY_list, vacXY_list)

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

    # Vacancy stats to match experiments
    for vacXY in vacXY_list: 
        if vacXY.light_up: 
            vacNeighbors_counts[vacXY.num_litUp_neighborsVacXY] += 1

            vac_a1_layer_counts[vacXY.a1_label] += 1
            vac_a2_layer_counts[vacXY.a2_label] += 1
            vac_a3_layer_counts[vacXY.a3_label] += 1

    energy = - sim_params['mu_In'] * num_cation - sim_params['mu_P'] * num_anion - sim_params['epsilon']*numNeighborDouble/2
    binding_energy = - sim_params['epsilon']*numNeighborDouble/2
    surface_energy = 2.0 * num_cation + 2.0 * num_anion - sim_params['epsilon']*numNeighborDouble/2

    if (calc_setting['verbosity']>0) and ((writeProjXY_filePrefix=="init") or (writeProjXY_filePrefix=="final")):
        write_projXY(siteXY_list, writeFilename=f"{calc_setting['calc_dir']}{writeProjXY_filePrefix}_projXY_atoms.xyz", mode='all_atoms')
        write_XY_sites_vac(siteXY_list, vacXY_list, writeFilename=f"{calc_setting['calc_dir']}{writeProjXY_filePrefix}_vacXY_lit.xyz", mode='lit_vac')
    if (calc_setting['verbosity']>1) and (writeProjXY_filePrefix=="init") :
        write_projXY(siteXY_list, writeFilename=f"{calc_setting['calc_dir']}{writeProjXY_filePrefix}_projXY_atoms_a1.xyz", mode='a1_label')
        write_projXY(siteXY_list, writeFilename=f"{calc_setting['calc_dir']}{writeProjXY_filePrefix}_projXY_atoms_a2.xyz", mode='a2_label')
        write_projXY(siteXY_list, writeFilename=f"{calc_setting['calc_dir']}{writeProjXY_filePrefix}_projXY_atoms_a3.xyz", mode='a3_label')

        write_XY_sites_vac(siteXY_list, vacXY_list, writeFilename=f"{calc_setting['calc_dir']}{writeProjXY_filePrefix}_vacXY_lit_a1.xyz", mode='lit_vac_a1')
        write_XY_sites_vac(siteXY_list, vacXY_list, writeFilename=f"{calc_setting['calc_dir']}{writeProjXY_filePrefix}_vacXY_lit_a2.xyz", mode='lit_vac_a2')
        write_XY_sites_vac(siteXY_list, vacXY_list, writeFilename=f"{calc_setting['calc_dir']}{writeProjXY_filePrefix}_vacXY_lit_a3.xyz", mode='lit_vac_a3')    

    ###############################################
    # XY projection geometry analysis
    occupied_coordXY, boundary_coordXY = extract_XY_occPoints_boundPoints(vacXY_list)
    print(f"Number of vacXY occupied sites {len(occupied_coordXY)}. Its boundary: {len(boundary_coordXY)}. ")

    if boundary_coordXY.size <=4:
        VacXY_roughness_perimeter = -1.0
        VacXY_rel_std = -1.0
    elif boundary_coordXY.ndim > 1:
        centroid = np.mean(boundary_coordXY, axis=0)
        distances = np.linalg.norm(boundary_coordXY - centroid, axis=1)
        VacXY_roughness_perimeter = np.std(distances)
        VacXY_rel_std = np.std(distances) / np.mean(distances)
    else:
        VacXY_roughness_perimeter = -1.0
        VacXY_rel_std = -1.0
    # print(f"{VacXY_roughness_perimeter:.5f}")

    if len(boundary_coordXY) <= 5:
        print(f"Insufficient points to compute ConvexHull: {len(boundary_coordXY)} points")
        VacXY_roughness_hull_perimeter = 0.0
        VacXY_roughness_hull_area = 0.0
        VacXY_cir_ratio = 0.0
    else: 
        try:
            hull = ConvexHull(boundary_coordXY)
            # print(hull.points[hull.vertices])

            # hull_perimeter = np.sum(np.linalg.norm(np.diff(hull.points[hull.vertices], axis=0), axis=1)) + np.linalg.norm(hull.points[hull.vertices][-1] - hull.points[hull.vertices][0])   # This is an equivalent way to calculate the hull perimeter
            hull_perimeter = hull.area
            hull_area = hull.volume
            # print(f"hull_perimeter = {hull_perimeter:.5f}")
            # print(f"hull_area = {hull_area:.5f}")

            # actual_perimeter = len(boundary_coordXY) * VacXY_neighbor_dist   # This is an equivalent way to calculate the actual perimeter
            actual_perimeter = VacXY_neighbor_dist * (vacNeighbors_counts[1] + vacNeighbors_counts[2] + vacNeighbors_counts[3] + vacNeighbors_counts[4] + vacNeighbors_counts[5])
            # print(f"actual_perimeter = {actual_perimeter:.5f}")

            # actual_area = len(occupied_coordXY) * VacXY_hex_pixel_area
            actual_area = VacXY_hex_pixel_area * (1/6*vacNeighbors_counts[2] + 1/3*vacNeighbors_counts[3] + 1/2*vacNeighbors_counts[4] + 2/3*vacNeighbors_counts[5] + vacNeighbors_counts[6])
            # print(f"actual_area = {actual_area:.5f}")

            VacXY_roughness_hull_perimeter = actual_perimeter / hull_perimeter 
            VacXY_roughness_hull_area = actual_area / hull_area
            VacXY_cir_ratio = 4*np.pi*actual_area / (actual_perimeter)**2
            print(f"VacXY_roughness_hull_perimeter (should always be >=1) = {VacXY_roughness_hull_perimeter:.5f}")
            print(f"VacXY_roughness_hull_area (should always be <=1) = {VacXY_roughness_hull_area:.5f}")
            print(f"VacXY_cir_ratio = {VacXY_cir_ratio:.5f}")

        except QhullError as e:
            print(f"ConvexHull failed. (Expected behavior for end of simulation) We are setting VacXY_xxx to placeholder values. ")
            VacXY_roughness_hull_perimeter = 0.0
            VacXY_roughness_hull_area = 0.0
            VacXY_cir_ratio = 0.0

    ###############################################
    # XZ projection geometry analysis
    occupied_coordXZ = []
    for site in site_list: 
        if site.has_atom: 
            occupied_coordXZ.append(site.real_space_coord[[0, 2]])
    occupied_coordXZ = np.array(occupied_coordXZ)

    if len(occupied_coordXZ) <= 3:
        XZ_bb_AR = 0.0
    else:
        min_x = np.min(occupied_coordXZ[:, 0])
        max_x = np.max(occupied_coordXZ[:, 0])
        min_z = np.min(occupied_coordXZ[:, 1])
        max_z = np.max(occupied_coordXZ[:, 1])
        XZ_bb_AR = (max_x - min_x) / (max_z - min_z)

    if len(occupied_coordXZ) <= 5:
        print(f"Insufficient points to compute ConvexHull: {len(occupied_coordXZ)} points")
        XZ_area_ratio = 0.0
    else: 
        try: 
            hull = ConvexHull(occupied_coordXZ)
            hull_area = hull.volume
            XZ_bb_area = (max_x - min_x) * (max_z - min_z)
            XZ_area_ratio = hull_area / XZ_bb_area
            print(f"XZ_area_ratio (should always be <=1) = {XZ_area_ratio:.5f}")
        except QhullError as e:
            print(f"ConvexHull failed. (Expected behavior for end of simulation) We are setting XZ_area_ratio to a placeholder value.")
            XZ_area_ratio = 0.0


    ###############################################
    # YZ projection geometry analysis
    occupied_coordYZ = []
    for site in site_list: 
        if site.has_atom: 
            occupied_coordYZ.append(site.real_space_coord[[1, 2]])
    occupied_coordYZ = np.array(occupied_coordYZ)

    if len(occupied_coordYZ) <= 2:
        YZ_bb_AR = 0.0
    else:
        min_y = np.min(occupied_coordYZ[:, 0])
        max_y = np.max(occupied_coordYZ[:, 0])
        min_z = np.min(occupied_coordYZ[:, 1])
        max_z = np.max(occupied_coordYZ[:, 1])
        YZ_bb_AR = (max_y - min_y) / (max_z - min_z)

    if len(occupied_coordYZ) <= 3:
        print(f"Insufficient points to compute ConvexHull: {len(occupied_coordYZ)} points")
        YZ_area_ratio = 0.0
    else: 
        try: 
            hull = ConvexHull(occupied_coordYZ)
            hull_area = hull.volume
            YZ_bb_area = (max_y - min_y) * (max_z - min_z)
            YZ_area_ratio = hull_area / YZ_bb_area
            print(f"YZ_area_ratio (should always be <=1) = {YZ_area_ratio:.5f}")
        except QhullError as e:
            print(f"ConvexHull failed. (Expected behavior for end of simulation) We are setting YZ_area_ratio to a placeholder value.")
            YZ_area_ratio = 0.0


    stats_dict = {
        "num_cation": num_cation, 
        "num_anion": num_anion, 
        "energy": energy, 
        "binding_energy": binding_energy, 
        "surface_energy": surface_energy, 
        "neighbor_3D_counts": neighbor_3D_counts, 
        "thickness": max_z - min_z, 
        "neighbor_XY_counts": neighbor_XY_counts, 
        "XY_a1_layer_counts": XY_a1_layer_counts, 
        "XY_a2_layer_counts": XY_a2_layer_counts, 
        "XY_a3_layer_counts": XY_a3_layer_counts, 
        "vacNeighbors_counts": vacNeighbors_counts, 
        "vac_a1_layer_counts": vac_a1_layer_counts, 
        "vac_a2_layer_counts": vac_a2_layer_counts, 
        "vac_a3_layer_counts": vac_a3_layer_counts, 
        "VacXY_roughness_perimeter": VacXY_roughness_perimeter, 
        "VacXY_rel_std": VacXY_rel_std, 
        "VacXY_roughness_hull_perimeter": VacXY_roughness_hull_perimeter, 
        "VacXY_roughness_hull_area": VacXY_roughness_hull_area, 
        "VacXY_cir_ratio": VacXY_cir_ratio,
        "XZ_bb_AR": XZ_bb_AR, 
        "XZ_area_ratio": XZ_area_ratio, 
        "YZ_bb_AR": YZ_bb_AR, 
        "YZ_area_ratio": YZ_area_ratio, 
    }
    # print(f"\t {stats_dict}")
    return stats_dict


def extract_XY_occPoints_boundPoints(vacXY_list): 
    occupied_coord = []
    boundary_coord = []
    for vacXY in vacXY_list: 
        if vacXY.light_up: 
            occupied_coord.append(vacXY.coordXY)

            if np.sum(np.array(vacXY.neighborXY_vac_litUp_bool))!=6:
                boundary_coord.append(vacXY.coordXY)

    return np.array(occupied_coord), np.array(boundary_coord)

