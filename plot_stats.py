import sys, time, os
import pickle, csv, gzip, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_singleValue(key, allData, num_bins, bins, bin_centers, ax, lineStyle, lineColor, fillColor, label=""): 
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
    
    ax.plot(bin_centers, bin_averages, linestyle=lineStyle, color=lineColor, label=label+key)
    ax.fill_between(bin_centers, bin_averages - bin_std_dev, bin_averages + bin_std_dev, color=fillColor, alpha=0.4)
    ax.grid(alpha=0.6)


def plot_num_neighbors(key, allData, num_bins, bins, bin_centers, ax): 
    num_neighbors = len(allData[0][0][key]) - 1
    print(f"{num_neighbors+1} number of neighbors means fully-bonded/occupied. This is thus omitted from the plot.")
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


def plot_theta(key, oneTrajData, num_bins, bins, bin_centers, deltaTBin, ax): 
    num_values = len(oneTrajData[0][key])
    cmap = plt.get_cmap('Reds')

    for i in range(num_values):  
        # Initialize arrays to store sum of values, sum of squares, and count of values in each bin
        bin_sums = np.zeros(num_bins)
        bin_counts = np.zeros(num_bins)
        
        # Assign values to bins and accumulate sums and counts
        times = np.array([item['simTime'] for item in oneTrajData])
        values = np.array([item[key][i] for item in oneTrajData])
        
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


def avg_one_dataset(dirPrefix, numRepeats, writeToDir, deltaTBin):
    max_simTime = 0.0
    allData = []
    for repeat in range(numRepeats):
        json_file_path = f"{dirPrefix}_repeat_{repeat+1}/stats.json"

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
    fig, axs = plt.subplots(2,1, figsize=(6,8))
    cmap = plt.get_cmap('tab20')
    plot_singleValue('num_cation', allData, num_bins, bins, bin_centers, axs[0], '-', cmap.colors[0], cmap.colors[1])
    plot_singleValue('num_anion', allData, num_bins, bins, bin_centers, axs[1], '-', cmap.colors[2], cmap.colors[3])
    fig.suptitle(f"Num atoms vs. Time, avg over {numRepeats} traj")
    axs[0].set(xlabel='Time', ylabel='Number of Atoms')
    axs[0].legend()
    axs[1].set(xlabel='Time', ylabel='Number of Atoms')
    axs[1].legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_nAtoms.pdf")
    plt.close()


    # plot energy
    fig, ax = plt.subplots(1,1, figsize=(6,4))
    plot_singleValue('energy', allData, num_bins, bins, bin_centers, ax, ':', cmap.colors[4], cmap.colors[5])
    ax.set(xlabel='Time', ylabel='Energy', title=f"Energy vs. Time, avg over {numRepeats} traj")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_energy.pdf")
    plt.close()
    

    # plot thickness over time
    fig, ax = plt.subplots(1,1, figsize=(8,5))
    plot_singleValue('thickness', allData, num_bins, bins, bin_centers, ax, '-', cmap.colors[6], cmap.colors[7])
    ax.set(xlabel='Time', ylabel='Thickness (nm)', title=f"Thickness vs. Time, avg over {numRepeats} traj")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_thickness.pdf")
    plt.close()


    # plot 3D number of neighbors
    fig, ax = plt.subplots(1,1, figsize=(6,5))
    plot_num_neighbors('neighbor_3D_counts', allData, num_bins, bins, bin_centers, ax)
    ax.set(xlabel='Time', ylabel='Number Atoms', title=f"3D num_neigh, avg over {numRepeats} traj")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_site_neigh_count.pdf")
    plt.close()

    '''
    # plot 2D number of neighbors
    fig, ax = plt.subplots(1,1, figsize=(6,5))
    plot_num_neighbors('neighbor_XY_counts', allData, num_bins, bins, bin_centers, ax)
    ax.set(xlabel='Time', ylabel='Number Atoms', title=f"siteXY num_neigh, avg over {numRepeats} traj")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_siteXY_neigh_count.pdf")
    '''

    # plot 2D vacancy number of neighbors
    fig, ax = plt.subplots(1,1, figsize=(6,5))
    plot_num_neighbors('vacNeighbors_counts', allData, num_bins, bins, bin_centers, ax)
    ax.set(xlabel='Time', ylabel='Number Atoms', title=f"VacXY num_neigh, avg over {numRepeats} traj")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_VacXY_neigh_count.pdf")
    plt.close()


    # plot 2D vacancy monolayer coverage, avg over a1, a2, a3 directions
    traj = allData[0]
    vac_a1_init = np.array(traj[0]['vac_a1_layer_counts'])
    # vac_a2_init = np.array(traj[0]['vac_a2_layer_counts'])
    # vac_a3_init = np.array(traj[0]['vac_a3_layer_counts'])
    for frame in traj: 
        # frame['vac_theta'] = (( np.array(frame['vac_a1_layer_counts'])/vac_a1_init + np.array(frame['vac_a2_layer_counts'])/vac_a2_init + np.array(frame['vac_a3_layer_counts'])/vac_a3_init ) / 3.0).tolist()

        counts = np.array(frame['vac_a1_layer_counts'])
        init = np.array(vac_a1_init)
        vac_theta = np.divide(counts, init, out=np.zeros_like(counts, dtype=float), where=init!=0)
        
        # Apply the special conditions
        vac_theta[(init == 0) & (counts == 0)] = 0
        vac_theta[(init == 0) & (counts != 0)] = 1
        
        frame['vac_theta'] = vac_theta.tolist()

    fig, ax = plt.subplots(1,1, figsize=(15,4))
    plot_theta('vac_theta', traj, num_bins, bins, bin_centers, deltaTBin, ax)
    ax.set(xlabel='Time', title=r"Traj repeat_1, $/theta$ in a1 direction")
    ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_VacXY_theta.pdf")
    plt.close()

    return


def plot_all_datasets(all_chemPot_sets, deltaT_list, numTraj, dirPrefix, writeToDir):
    # all_chemPot_sets = ['-4.0_-4.0', '-3.5_-3.5', '-3.0_-3.0', '-2.5_-2.5', '-2.0_-2.0']
    # deltaT_list = [0.00001, 0.0001, 0.0005, 0.01, 0.01]
    # numTraj = 16
    # dirPrefix = 'CALCS_diam15.4_thick3.5/'
    # writeToDir = './'

    fig1, axs1 = plt.subplots(2,1, figsize=(6,8))
    fig2, ax2 = plt.subplots(1,1, figsize=(6,4))
    fig3, ax3 = plt.subplots(1,1, figsize=(8,5))
    fig4, ax4 = plt.subplots(1,1, figsize=(6,4))

    cmap = plt.get_cmap('tab20')
    
    for idx, chemPot in enumerate(all_chemPot_sets): 
        max_simTime = 0.0
        allData = []
        for repeat in range(numTraj):
            json_file_path = f"{dirPrefix}{chemPot}_repeat_{repeat+1}/stats.json"

            with open(json_file_path, 'r') as file:
                stats_list = json.load(file)
            allData.append(stats_list)

            if stats_list[-1]['simTime'] > max_simTime: 
                max_simTime = stats_list[-1]['simTime']

        num_bins = int(np.ceil((max_simTime) / deltaT_list[idx])) + 10
        bins = np.linspace(0.0, max_simTime, num_bins + 1)
        bin_centers = (bins[:-1] + bins[1:]) / 2
        print(f"\nThe maximium simulation time is {max_simTime:.3e}. And we have split it into {num_bins} bins for statistics plotting. \n")

        # plot number of cations, number of anions, fig1
        plot_singleValue('num_cation', allData, num_bins, bins, bin_centers, axs1[0], '-', cmap.colors[idx*4], cmap.colors[idx*4+1], label=f"{chemPot} ")
        plot_singleValue('num_anion', allData, num_bins, bins, bin_centers, axs1[1], '-', cmap.colors[idx*4+2], cmap.colors[idx*4+3], label=f"{chemPot} ")
        axs1[0].set(xlabel='Time', ylabel='Number of Atoms')
        axs1[0].set_xscale('log')
        axs1[0].legend()
        axs1[1].set(xlabel='Time', ylabel='Number of Atoms')
        axs1[1].set_xscale('log')
        axs1[1].legend()

        # plot surface area, fig4
        for traj in allData: 
            for frame in traj: 
                frame['VacXY_tot'] = np.sum(np.array(frame['vacNeighbors_counts']))
        plot_singleValue('VacXY_tot', allData, num_bins, bins, bin_centers, ax4, '-', cmap.colors[idx*2], cmap.colors[idx*2+1], label=f"{chemPot} ")
        ax4.set(xlabel='Time', ylabel='Surface Area [unit]', title=f"Surface Area vs. Time, avg over {numTraj} traj")
        ax4.set_xscale('log')
        ax4.legend()


        # plot energy, fig2
        plot_singleValue('energy', allData, num_bins, bins, bin_centers, ax2, '-', cmap.colors[idx*2], cmap.colors[idx*2+1], label=f"{chemPot} ")
        ax2.set(xlabel='Time', ylabel='Energy', title=f"Energy vs. Time, avg over {numTraj} traj")
        ax2.set_xscale('log')
        ax2.legend()
        
        # plot thickness over time, fig3
        plot_singleValue('thickness', allData, num_bins, bins, bin_centers, ax3, '-', cmap.colors[idx*2], cmap.colors[idx*2+1], label=f"{chemPot} ")
        ax3.set(xlabel='Time', ylabel='Thickness (nm)', title=f"Thickness vs. Time, avg over {numTraj} traj")
        ax3.set_xscale('log')
        ax3.legend()

    fig1.suptitle(f"Num atoms vs. Time, avg over {numTraj} traj")
    fig1.tight_layout()
    fig1.savefig(f"{writeToDir}plot_nAtoms.pdf")
    fig4.tight_layout()
    fig4.savefig(f"{writeToDir}plot_surfaceArea.pdf")
    fig2.tight_layout()
    fig2.savefig(f"{writeToDir}plot_energy.pdf")
    fig3.tight_layout()
    fig3.savefig(f"{writeToDir}plot_thickness.pdf")
    return


######################################################
all_chemPot_sets = ['-4.0_-4.0', '-3.0_-3.0', '-2.5_-2.5', '-2.0_-2.0'] 
# ['-4.0_-4.0', '-3.5_-3.5', '-3.0_-3.0', '-2.5_-2.5', '-2.0_-2.0']
deltaT_list = [0.00002, 0.002, 0.02, 1.5]
# [0.00001, 0.0001, 0.0005, 0.01, 1.0]
numTraj = 16

for idx, chemPot in enumerate(all_chemPot_sets): 
    avg_one_dataset(f'CALCS_diam15.4_thick3.5/{chemPot}', numTraj, f'CALCS_diam15.4_thick3.5/{chemPot}/', deltaT_list[idx])

plot_all_datasets(all_chemPot_sets, deltaT_list, numTraj, 'CALCS_diam15.4_thick3.5/', './')
