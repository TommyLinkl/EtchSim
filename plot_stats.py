import sys, time, os
import pickle, csv, gzip, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_singleValue(key, allData, num_bins, bins, bin_centers, ax, lineStyle, lineColor, fillColor): 
    # Initialize arrays to store sum of values, sum of squares, and count of values in each bin
    bin_sums = np.zeros(num_bins)
    bin_squares = np.zeros(num_bins)
    bin_counts = np.zeros(num_bins)
    
    # Assign values to bins and accumulate sums and counts
    for traj in allData:
        times = np.array([item['simTime'] for item in traj])
        values = np.array([item[key] for item in traj])
        
        bin_indices = np.digitize(times, bins) - 1  # Find the bin index for each time point

        # Filter valid bin indices
        valid_indices = (0 <= bin_indices) & (bin_indices < num_bins)

        # Accumulate sums, squares, and counts for valid indices
        np.add.at(bin_sums, bin_indices[valid_indices], values[valid_indices])
        np.add.at(bin_squares, bin_indices[valid_indices], values[valid_indices]**2)
        np.add.at(bin_counts, bin_indices[valid_indices], 1)
    
    bin_averages = np.divide(bin_sums, bin_counts, out=np.zeros_like(bin_sums), where=bin_counts != 0)

    bin_counts_non_zero = np.where(bin_counts > 1, bin_counts, np.nan)  # Avoid division by zero
    bin_squares_divided = np.divide(bin_squares, bin_counts_non_zero, out=np.zeros_like(bin_squares), where=~np.isnan(bin_counts_non_zero))
    bin_variance = np.where(bin_counts > 1, bin_squares_divided - bin_averages**2, 0)
    bin_variance = np.maximum(bin_variance, 0)
    bin_std_dev = np.sqrt(bin_variance)
    
    ax.plot(bin_centers, bin_averages, linestyle=lineStyle, color=lineColor, label=key)
    ax.fill_between(bin_centers, bin_averages - bin_std_dev, bin_averages + bin_std_dev, color=fillColor, alpha=0.4)
    ax.grid(alpha=0.6)


def plot_num_neighbors(key, allData, num_bins, bins, bin_centers, ax): 
    num_neighbors = len(allData[0][0][key])
    cmap = plt.get_cmap('tab20')

    for i in range(num_neighbors): 
        # Initialize arrays to store sum of values, sum of squares, and count of values in each bin
        bin_sums = np.zeros(num_bins)
        bin_squares = np.zeros(num_bins)
        bin_counts = np.zeros(num_bins)
        
        # Assign values to bins and accumulate sums and counts
        for traj in allData:
            times = np.array([item['simTime'] for item in traj])
            values = np.array([item[key][i] for item in traj])
            
            bin_indices = np.digitize(times, bins) - 1  # Find the bin index for each time point

            # Filter valid bin indices
            valid_indices = (0 <= bin_indices) & (bin_indices < num_bins)

            # Accumulate sums, squares, and counts for valid indices
            np.add.at(bin_sums, bin_indices[valid_indices], values[valid_indices])
            np.add.at(bin_squares, bin_indices[valid_indices], values[valid_indices]**2)
            np.add.at(bin_counts, bin_indices[valid_indices], 1)
            
        bin_averages = np.divide(bin_sums, bin_counts, out=np.zeros_like(bin_sums), where=bin_counts != 0)

        bin_counts_non_zero = np.where(bin_counts > 1, bin_counts, np.nan)  # Avoid division by zero
        bin_squares_divided = np.divide(bin_squares, bin_counts_non_zero, out=np.zeros_like(bin_squares), where=~np.isnan(bin_counts_non_zero))
        bin_variance = np.where(bin_counts > 1, bin_squares_divided - bin_averages**2, 0)
        bin_variance = np.maximum(bin_variance, 0)
        bin_std_dev = np.sqrt(bin_variance)
        
        ax.plot(bin_centers, bin_averages, linestyle='-', color=cmap.colors[(2*i)%20], label=f"{i}_neighbors")
        ax.fill_between(bin_centers, bin_averages - bin_std_dev, bin_averages + bin_std_dev, color=cmap.colors[(2*i)%20], alpha=0.4)
    ax.grid(alpha=0.6)


def plot_theta(key, allData, num_bins, bins, bin_centers, deltaTBin, ax): 
    num_values = len(allData[0][0][key])
    cmap = plt.get_cmap('Reds')

    for i in range(num_values):  
        # Initialize arrays to store sum of values, sum of squares, and count of values in each bin
        bin_sums = np.zeros(num_bins)
        bin_counts = np.zeros(num_bins)
        
        # Assign values to bins and accumulate sums and counts
        for traj in allData:
            times = np.array([item['simTime'] for item in traj])
            values = np.array([item[key][i] for item in traj])
            
            bin_indices = np.digitize(times, bins) - 1  # Find the bin index for each time point

            # Filter valid bin indices
            valid_indices = (0 <= bin_indices) & (bin_indices < num_bins)

            # Accumulate sums and counts for valid indices
            np.add.at(bin_sums, bin_indices[valid_indices], values[valid_indices])
            np.add.at(bin_counts, bin_indices[valid_indices], 1)
            
        bin_averages = np.divide(bin_sums, bin_counts, out=np.zeros_like(bin_sums), where=bin_counts != 0)
        # print(len(bin_averages))
        # print(num_bins)
        

        for binIdx in range(num_bins): 
            if bin_averages[binIdx]==1.0:
                color = 'b'
            else:
                color = cmap(bin_averages[binIdx])
            rect = plt.Rectangle((bins[binIdx], i*deltaTBin), deltaTBin, deltaTBin, facecolor=color, edgecolor='grey', lw=0.05)
            ax.add_patch(rect)

    ax.set(xlim=(min(bins)-10*deltaTBin, max(bins)+10*deltaTBin), ylim=(0-5*deltaTBin, num_values*deltaTBin+5*deltaTBin))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical')
    cbar.set_label(r"$\theta$")
    ax.axis('equal')


def plot_all_stats(dirPrefix, numRepeats, writeToDir, deltaTBin):
    max_simTime = 0.0
    allData = []
    for repeat in range(numRepeats):
        json_file_path = f"{dirPrefix}_repeat_{repeat+1}/postproc_stats.json"

        with open(json_file_path, 'r') as file:
            stats_list = json.load(file)
        allData.append(stats_list)

        if stats_list[-1]['simTime'] > max_simTime: 
            max_simTime = stats_list[-1]['simTime']

    num_bins = int(np.ceil((max_simTime) / deltaTBin)) + 10
    bins = np.linspace(0.0, max_simTime, num_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    print(f"\nThe maximium simulation time is {max_simTime:.3e}. And we have split it into {num_bins} bins for statistics plotting. \n")
    # print(allData[0][0])
    # print(allData[0][-4])


    # plot number of cations, number of anions
    fig, axs = plt.subplots(2,1, figsize=(8, 10))
    cmap = plt.get_cmap('tab20')
    plot_singleValue('num_cation', allData, num_bins, bins, bin_centers, axs[0], '-', cmap.colors[0], cmap.colors[1])
    plot_singleValue('num_anion', allData, num_bins, bins, bin_centers, axs[1], '-', cmap.colors[2], cmap.colors[3])
    fig.suptitle(f"Num atoms vs. Time averaged over {numRepeats} trajectories")
    axs[0].set(xlabel='Time', ylabel='Number of Atoms')
    axs[0].legend()
    axs[1].set(xlabel='Time', ylabel='Number of Atoms')
    axs[1].legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_nAtoms.pdf")


    # plot number of cations, number of anions
    fig, ax = plt.subplots(1,1, figsize=(8,5))
    plot_singleValue('energy', allData, num_bins, bins, bin_centers, ax, ':', cmap.colors[4], cmap.colors[5])
    ax.set(xlabel='Time', ylabel='Energy', title=f"Energy vs. Time averaged over {numRepeats} trajectories")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_energy.pdf")
    

    # plot thickness over time
    fig, ax = plt.subplots(1,1, figsize=(8,5))
    plot_singleValue('thickness', allData, num_bins, bins, bin_centers, ax, '-', cmap.colors[6], cmap.colors[7])
    ax.set(xlabel='Time', ylabel='Thickness (nm)', title=f"Thickness vs. Time averaged over {numRepeats} trajectories")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_thickness.pdf")


    # plot 3D number of neighbors
    fig, ax = plt.subplots(1,1, figsize=(8,5))
    plot_num_neighbors('neighbor_3D_counts', allData, num_bins, bins, bin_centers, ax)
    ax.set(xlabel='Time', ylabel='Number Atoms', title=f"Averaged over {numRepeats} trajectories")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_3D_neigh_count.pdf")


    # plot 2D number of neighbors
    fig, ax = plt.subplots(1,1, figsize=(8,5))
    plot_num_neighbors('neighbor_XY_counts', allData, num_bins, bins, bin_centers, ax)
    ax.set(xlabel='Time', ylabel='Number Atoms', title=f"Averaged over {numRepeats} trajectories")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_2D_neigh_count.pdf")


    # plot 2D vacancy number of neighbors
    fig, ax = plt.subplots(1,1, figsize=(8,5))
    plot_num_neighbors('vacNeighbors_counts', allData, num_bins, bins, bin_centers, ax)
    ax.set(xlabel='Time', ylabel='Number Atoms', title=f"Averaged over {numRepeats} trajectories")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_2D_Vac_neigh_count.pdf")


    # plot 2D vacancy monolayer coverage, averaged over a1, a2, a3 directions
    for traj in allData:
        vac_a1_init = np.array(traj[0]['vac_a1_layer_counts'])
        vac_a2_init = np.array(traj[0]['vac_a2_layer_counts'])
        vac_a3_init = np.array(traj[0]['vac_a3_layer_counts'])
        for frame in traj: 
            frame['vac_theta'] = (( np.array(frame['vac_a1_layer_counts'])/vac_a1_init + np.array(frame['vac_a2_layer_counts'])/vac_a2_init + np.array(frame['vac_a3_layer_counts'])/vac_a3_init ) / 3.0).tolist()
    # print(allData[0][-4])
    fig, ax = plt.subplots(1,1, figsize=(15,5))
    plot_theta('vac_theta', allData, num_bins, bins, bin_centers, deltaTBin, ax)
    ax.set(xlabel='Time') # , ylabel='Number Atoms', title=f"Averaged over {numRepeats} trajectories")
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_2D_Vac_theta.pdf")

    
    return


######################################################
# plot_all_stats('CALCS_diam15.4_thick3.5/-4.0_-4.0', 32, 'CALCS_diam15.4_thick3.5/-4.0_-4.0/', 0.00001)
plot_all_stats('CALCS_diam15.4_thick3.5/-3.0_-3.0', 32, 'CALCS_diam15.4_thick3.5/-3.0_-3.0/', 0.0005)