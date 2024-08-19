import sys, time, os, shutil
import pickle, csv, gzip, json
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import rc
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['font.sans-serif'] = "Liberation Sans"
matplotlib.rcParams['font.size'] = 10
# matplotlib.rcParams['mathtext.fontset'] = "custom"
# matplotlib.rcParams['mathtext.rm'] = "Liberation Sans"
# matplotlib.rcParams['mathtext.it'] = "Liberation Sans:italic"
# matplotlib.rcParams['mathtext.bf'] = "Liberation Sans:bold"
rc('axes', linewidth=1.0)
matplotlib.rcParams['xtick.major.width'] = 0.5
matplotlib.rcParams['ytick.major.width'] = 0.5
matplotlib.rcParams['xtick.major.size'] = 7
matplotlib.rcParams['ytick.major.size'] = 7
matplotlib.rcParams['xtick.minor.size'] = 4
matplotlib.rcParams['ytick.minor.size'] = 4
matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

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
    
    # Remove bins where counts are zero
    non_zero_mask = bin_counts > 0
    bin_averages = np.divide(bin_sums[non_zero_mask], bin_counts[non_zero_mask])
    bin_squares_divided = np.divide(bin_squares[non_zero_mask], bin_counts[non_zero_mask]) 
    bin_variance = bin_squares_divided - bin_averages**2
    bin_variance = np.where(np.abs(bin_variance) < 1e-10, 0.0, bin_variance)
    bin_std_dev = np.sqrt(bin_variance)
    
    ax.plot(bin_centers[non_zero_mask], bin_averages, linestyle=lineStyle, color=lineColor, label=label)
    ax.fill_between(bin_centers[non_zero_mask], bin_averages - bin_std_dev, bin_averages + bin_std_dev, color=fillColor, alpha=0.4)
    ax.grid(alpha=0.6)


def plot_num_neighbors(key, allData, num_bins, bins, bin_centers, ax): 
    num_neighbors = len(allData[0][0][key]) - 1
    # print(f"{num_neighbors+1} number of neighbors means fully-bonded/occupied. This is thus omitted from the plot.")
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
            
        # Remove bins where counts are zero
        non_zero_mask = bin_counts > 0
        bin_averages = np.divide(bin_sums[non_zero_mask], bin_counts[non_zero_mask])
        bin_squares_divided = np.divide(bin_squares[non_zero_mask], bin_counts[non_zero_mask]) 
        bin_variance = bin_squares_divided - bin_averages**2
        bin_variance = np.where(np.abs(bin_variance) < 1e-10, 0.0, bin_variance)
        bin_std_dev = np.sqrt(bin_variance)
        
        ax.plot(bin_centers[non_zero_mask], bin_averages, linestyle='-', color=cmap.colors[(2*i)%20], label=f"{i}_neighbors")
        ax.fill_between(bin_centers[non_zero_mask], bin_averages - bin_std_dev, bin_averages + bin_std_dev, color=cmap.colors[(2*i)%20], alpha=0.4)
    ax.grid(alpha=0.6)


def plot_theta(key, oneTrajData, num_bins, bins, bin_centers, deltaTBin, ax): 
    num_values = len(oneTrajData[0][key])
    cmap = plt.get_cmap('Reds')

    for i in range(4, num_values-4):  
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


def reverse_plot_theta(key, oneTrajData, num_bins, bins, bin_centers, deltaTBin, ax): 
    num_values = len(oneTrajData[0][key])
    # cmap = plt.get_cmap('bwr')
    jet = plt.get_cmap('jet')
    jet_colors = jet(np.linspace(0, 1, 256))
    custom_colors = np.concatenate([
        jet_colors[::-1], 
        jet_colors
    ])
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_jet", custom_colors, N=512)

    norm = mcolors.Normalize(vmin=4, vmax=num_values-5)

    for i in range(4, num_values-4):  
        bin_sums = np.zeros(num_bins)
        bin_counts = np.zeros(num_bins)
        
        # Assign values to bins and accumulate sums and counts
        times = np.array([item['simTime'] for item in oneTrajData])
        values = np.array([item[key][i] for item in oneTrajData])
        # print(values)
        
        bin_indices = np.digitize(times, bins) - 1  # Find the bin index for each time point

        # Filter valid bin indices
        valid_indices = (0 <= bin_indices) & (bin_indices < num_bins)

        # Accumulate sums and counts for valid indices
        np.add.at(bin_sums, bin_indices[valid_indices], values[valid_indices])
        np.add.at(bin_counts, bin_indices[valid_indices], 1)
        
        bin_averages = np.divide(bin_sums, bin_counts, out=np.zeros_like(bin_sums), where=bin_counts != 0)

        color = cmap(norm(i))
        # print(bin_centers)
        # print(bin_averages)
        if not np.all(bin_averages == 0.0): 
            if i>=num_values/2:
                ax.plot(bin_centers, bin_averages, "-", color=color, linewidth=1.0)
            else: 
                ax.plot(bin_centers, bin_averages, ":", color=color, linewidth=1.0)

    ax.set(ylabel=r"$\theta$", xlim=(min(bin_centers)-10*deltaTBin, max(bin_centers)+10*deltaTBin), ylim=(-0.1, 1.1))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, orientation='vertical')
    cbar.set_label(r"$a_1$")


def avg_one_dataset(dirPrefix, numRepeats, writeToDir, deltaTBin, color_line, color_shade): 
    matplotlib.rcParams['font.size'] = 8
    min_simTime = np.inf
    allData = []
    for repeat in range(numRepeats):
        json_file_path = f"{dirPrefix}_repeat_{repeat+1}/stats.json"

        with open(json_file_path, 'r') as file:
            stats_list = json.load(file)
        allData.append(stats_list)

        if stats_list[-1]['simTime'] < min_simTime: 
            min_simTime = stats_list[-1]['simTime']

    num_bins = int(np.floor((min_simTime) / deltaTBin))
    bins = np.linspace(0.0, min_simTime, num_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    print(f"{dirPrefix}\nThe shortest simulation time is {min_simTime:.3e}. And we have split it into {num_bins} bins for statistics plotting. \n")
    # print(allData[0][0])
    # print(allData[0][-4])

    # plot number of cations, number of anions
    fig, ax = plt.subplots(1,1, figsize=(3,3))
    cmap = plt.get_cmap('tab20')
    plot_singleValue('num_cation', allData, num_bins, bins, bin_centers, ax, '-', color_line, color_shade)
    # plot_singleValue('num_anion', allData, num_bins, bins, bin_centers, axs[1], '-', cmap.colors[2], cmap.colors[3])
    fig.suptitle(f"Num cations (or anions), {dirPrefix.split('/')[-1]}")
    ax.set(xlabel='Time', ylabel='Number of Atoms')
    # ax.legend()
    # axs[1].set(xlabel='Time', ylabel='Number of Atoms')
    # axs[1].legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_nAtoms.pdf")
    plt.close()

    # plot energy
    fig, ax = plt.subplots(1,1, figsize=(3,3))
    plot_singleValue('energy', allData, num_bins, bins, bin_centers, ax, '-', color_line, color_shade)
    ax.set(xlabel='Time', ylabel='Energy', title=f"Energy, {dirPrefix.split('/')[-1]}")
    # ax.legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_energy.pdf")
    plt.close()
    
    # plot surface area
    fig, ax = plt.subplots(1,1, figsize=(3,3))
    for traj in allData: 
        for frame in traj: 
            frame['VacXY_tot'] = np.sum(np.array(frame['vacNeighbors_counts']))
    plot_singleValue('VacXY_tot', allData, num_bins, bins, bin_centers, ax, '-', color_line, color_shade)
    ax.set(xlabel='Time', ylabel='Surface Area [unit]', title=f"Surface Area, {dirPrefix.split('/')[-1]}")
    # ax.legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_surfaceArea.pdf")
    plt.close()

    # plot thickness over time
    fig, ax = plt.subplots(1,1, figsize=(3,3))
    plot_singleValue('thickness', allData, num_bins, bins, bin_centers, ax, '-', color_line, color_shade)
    ax.set(xlabel='Time', ylabel='Thickness (nm)', title=f"Thickness, {dirPrefix.split('/')[-1]}")
    # ax.legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_thickness.pdf")
    plt.close()

    # plot 3D number of neighbors
    fig, ax = plt.subplots(1,1, figsize=(3,3))
    plot_num_neighbors('neighbor_3D_counts', allData, num_bins, bins, bin_centers, ax)
    ax.set(xlabel='Time', ylabel='Number Atoms', title=f"3D num_neigh, {dirPrefix.split('/')[-1]}")
    # ax.legend()
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
    fig, ax = plt.subplots(1,1, figsize=(4,4))
    plot_num_neighbors('vacNeighbors_counts', allData, num_bins, bins, bin_centers, ax)
    ax.set(xlabel='Time', ylabel='Number Atoms', title=f"VacXY num_neigh, {dirPrefix.split('/')[-1]}")
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

    fig, ax = plt.subplots(1,1, figsize=(11,3.5))
    plot_theta('vac_theta', traj, num_bins, bins, bin_centers, deltaTBin, ax)
    ax.set(xlabel='Time', title=rf"Traj repeat_1, $\theta$ in $a_1$ direction. {dirPrefix.split('/')[-1]}")
    ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_VacXY_theta.pdf")
    plt.close()

    fig, ax = plt.subplots(1,1, figsize=(8,4))
    reverse_plot_theta('vac_theta', traj, num_bins, bins, bin_centers, deltaTBin, ax)
    ax.set(xlabel='Time', title=rf"Traj repeat_1, $\theta$ in $a_1$ direction. {dirPrefix.split('/')[-1]}")
    ax.set_yticks([])
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_VacXY_theta_reverse.pdf")
    plt.close()
    matplotlib.rcParams['font.size'] = 10

    return


def plot_all_datasets(all_chemPot_sets, deltaT_list, numTraj, dirPrefix, writeToDir):
    fig1, ax1 = plt.subplots(1,1, figsize=(5,3))
    fig2, ax2 = plt.subplots(1,1, figsize=(5,3))
    fig3, ax3 = plt.subplots(1,1, figsize=(5,3))
    fig4, ax4 = plt.subplots(1,1, figsize=(5,3))
    print("Summarizing all plots... ")

    cmap = plt.get_cmap('tab20')
    
    for idx, chemPot in enumerate(all_chemPot_sets): 
        min_simTime = np.inf
        allData = []
        for repeat in range(numTraj):
            json_file_path = f"{dirPrefix}{chemPot}_repeat_{repeat+1}/stats.json"

            with open(json_file_path, 'r') as file:
                stats_list = json.load(file)
            allData.append(stats_list)

            if stats_list[-1]['simTime'] < min_simTime: 
                min_simTime = stats_list[-1]['simTime']

        num_bins = int(np.floor((min_simTime) / deltaT_list[idx]))
        bins = np.linspace(0.0, min_simTime, num_bins + 1)
        bin_centers = (bins[:-1] + bins[1:]) / 2

        # plot number of cations, number of anions, fig1
        plot_singleValue('num_cation', allData, num_bins, bins, bin_centers, ax1, '-', cmap.colors[idx*2], cmap.colors[idx*2+1], label=f"{chemPot}")
        # plot_singleValue('num_anion', allData, num_bins, bins, bin_centers, ax1, ':', cmap.colors[idx*4], cmap.colors[idx*4+1], label=f"{chemPot}")
        ax1.set(xlabel='Time', ylabel=f"Num cations (or anions), avg over {numTraj} traj", ylim=(0, 12000))
        ax1.set_xscale('log')
        ax1.legend()

        # plot surface area, fig4
        for traj in allData: 
            for frame in traj: 
                frame['VacXY_tot'] = np.sum(np.array(frame['vacNeighbors_counts']))
        plot_singleValue('VacXY_tot', allData, num_bins, bins, bin_centers, ax4, '-', cmap.colors[idx*2], cmap.colors[idx*2+1], label=f"{chemPot} ")
        ax4.set(xlabel='Time', ylabel='Surface Area [unit]', title=f"Surface Area, avg over {numTraj} traj", ylim=(0, 1100))
        ax4.set_xscale('log')
        ax4.legend()

        # plot energy, fig2
        plot_singleValue('energy', allData, num_bins, bins, bin_centers, ax2, '-', cmap.colors[idx*2], cmap.colors[idx*2+1], label=f"{chemPot} ")
        ax2.set(xlabel='Time', ylabel='Energy', title=f"Energy, avg over {numTraj} traj", ylim=(0, 50000))
        ax2.set_xscale('log')
        ax2.legend()
        
        # plot thickness over time, fig3
        plot_singleValue('thickness', allData, num_bins, bins, bin_centers, ax3, '-', cmap.colors[idx*2], cmap.colors[idx*2+1], label=f"{chemPot} ")
        ax3.set(xlabel='Time', ylabel='Thickness (nm)', title=f"Thickness, avg over {numTraj} traj", ylim=(0,40))
        ax3.set_xscale('log')
        ax3.legend()

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
all_chemPot_sets = ['-4.0_-4.0', '-3.5_-3.5', '-3.0_-3.0', '-2.5_-2.5', '-2.0_-2.0'] 
deltaT_list = [0.00003, 0.0003, 0.003, 0.05, 2.8]
numTraj = 16
cmap = plt.get_cmap('tab20')

for idx, chemPot in enumerate(all_chemPot_sets): 
    avg_one_dataset(f'CALCS_diam15.4_thick3.5/{chemPot}', numTraj, f'CALCS_diam15.4_thick3.5/{chemPot}/', deltaT_list[idx], cmap.colors[idx*2], cmap.colors[idx*2+1])

    shutil.copy(f'CALCS_diam15.4_thick3.5/{chemPot}/plot_energy.pdf', f"PLOTS/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_energy.pdf")
    shutil.copy(f'CALCS_diam15.4_thick3.5/{chemPot}/plot_nAtoms.pdf', f"PLOTS/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_nAtoms.pdf")
    shutil.copy(f'CALCS_diam15.4_thick3.5/{chemPot}/plot_site_neigh_count.pdf', f"PLOTS/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_site_neigh_count.pdf")
    shutil.copy(f'CALCS_diam15.4_thick3.5/{chemPot}/plot_surfaceArea.pdf', f"PLOTS/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_surfaceArea.pdf")
    shutil.copy(f'CALCS_diam15.4_thick3.5/{chemPot}/plot_thickness.pdf', f"PLOTS/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_thickness.pdf")
    shutil.copy(f'CALCS_diam15.4_thick3.5/{chemPot}/plot_VacXY_neigh_count.pdf', f"PLOTS/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_VacXY_neigh_count.pdf")
    shutil.copy(f'CALCS_diam15.4_thick3.5/{chemPot}/plot_VacXY_theta.pdf', f"PLOTS/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_VacXY_theta.pdf")
    shutil.copy(f'CALCS_diam15.4_thick3.5/{chemPot}/plot_VacXY_theta_reverse.pdf', f"PLOTS/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_VacXY_theta_reverse.pdf")

plot_all_datasets(all_chemPot_sets, deltaT_list, numTraj, 'CALCS_diam15.4_thick3.5/', 'PLOTS/')