import sys, time, os
import pickle, csv, gzip, json
import pandas as pd
from src.constants import *
from main import load_input_from_json
from src.sites import update_whole_lattice_iteration, update_XY_projection, update_XYvac
from src.kmc import collect_stats

def postprocessing():
    if len(sys.argv) != 2:
        print("Usage: python postproc.py CALC_DIR/")
        return
    calc_setting['calc_dir'] = sys.argv[1]
    load_input_from_json(f"{sys.argv[1]}input.json")

    # Read in site_list, siteXY_list, vacXY_list from pickle files and construct these lists
    start_time = time.time()
    with open(f'{calc_setting['calc_dir']}init_site_list.pkl', 'rb') as f:
        wz_lattice = pickle.load(f)
    with open(f'{calc_setting['calc_dir']}init_siteXY_list.pkl', 'rb') as f:
        wz_lattice_XY = pickle.load(f)
    with open(f'{calc_setting['calc_dir']}init_vacXY_list.pkl', 'rb') as f:
        wz_lattice_vacXY = pickle.load(f)

    # Read in forStatsLater.csv.gz using pandas
    data_df = pd.read_csv(f"{calc_setting['calc_dir']}forStatsLater.csv.gz", compression='gzip')
    # print(data_df.columns.tolist())
    
    # Check if each data row is of length len(wz_lattice) + 2
    expected_length = len(wz_lattice) + 2
    if not all(len(row) == expected_length for row in data_df.itertuples(index=False)):
        raise ValueError("Error: Data row length does not match expected length.")
    end_time = time.time()
    print(f"\nDone with reading struct_lists and data. Elapsed time: {(end_time - start_time):.5f}s")

    # For each frame, update the lattice and collect stats
    start_time = time.time()
    stats_list = []
    for row in data_df.itertuples(index=False):
        step_num = int(row[0])
        simTime = float(row[1])

        for site in wz_lattice:
            site.iteration_changed = True
            site.has_atom = bool(int(row[2+site.wz_lattice_idx]))
        update_whole_lattice_iteration(wz_lattice, sim_params)
        update_XY_projection(wz_lattice, wz_lattice_XY)
        update_XYvac(wz_lattice_XY, wz_lattice_vacXY)
        
        # Collect stats
        stats = collect_stats(wz_lattice, wz_lattice_XY, wz_lattice_vacXY, sim_params, writeProjXY_filePrefix=f"step_{step_num}")
        stats['stepNum'] = step_num
        stats['simTime'] = simTime
        stats_list.append(stats)

    # Dump the stats in a json file
    with open(f"{calc_setting['calc_dir']}postproc_stats.json", 'w') as file:
        json.dump(stats_list, file, separators=(',', ':'))
    end_time = time.time()
    print(f"\nDone with post-processing stats. Elapsed time: {(end_time - start_time):.5f}s")


######################################################
postprocessing()