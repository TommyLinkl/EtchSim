import shutil
import json
import numpy as np
from scipy.ndimage import uniform_filter1d
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import rc
# matplotlib.rcParams['font.family'] = "sans-serif"
# matplotlib.rcParams['font.sans-serif'] = "Liberation Sans"
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
from src.constants import *
from plot_stats import plot_singleValue, plot_num_neighbors, plot_theta, reverse_plot_theta


def uneven_avg_one_dataset(dirPrefix, numRepeats, writeToDir, deltaTBin, color_line, color_shade): 
    # This function is for creating plots for each parameter set, where all 16 trajectories are time-averaged using bins. 
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
    print(f"{dirPrefix}\nThe shortest simulation time is {min_simTime:.3e}. And we have split it into {num_bins} bins for statistics plotting.")
    # print(allData[0][0])
    # print(allData[0][-4])

    # plot number of cations, number of anions
    fig, ax = plt.subplots(1,1, figsize=(3,3))
    cmap = plt.get_cmap('tab20')
    plot_singleValue('num_cation', allData, num_bins, bins, bin_centers, ax, '-', color_line, color_shade)
    plot_singleValue('num_anion', allData, num_bins, bins, bin_centers, ax, ':', color_line, color_shade)
    fig.suptitle(f"Num cations (or anions), {dirPrefix.split('/')[-1]}")
    ax.set(xlabel='Time', ylabel='Number of Atoms', ylim=(0, 12000))
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_nAtoms.pdf")
    plt.close()

    # plot energy
    fig, ax = plt.subplots(1,1, figsize=(3,3))
    plot_singleValue('energy', allData, num_bins, bins, bin_centers, ax, '-', color_line, color_shade)
    ax.set(xlabel='Time', ylabel='Total Energy', title=f"Total Energy, {dirPrefix.split('/')[-1]}")
    # ax.legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_energy.pdf")
    plt.close()
    fig, ax = plt.subplots(1,1, figsize=(3,3))
    plot_singleValue('binding_energy', allData, num_bins, bins, bin_centers, ax, '-', color_line, color_shade)
    ax.set(xlabel='Time', ylabel='Binding Energy', title=f"Binding Energy, {dirPrefix.split('/')[-1]}")
    # ax.legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_binding_energy.pdf")
    plt.close()
    fig, ax = plt.subplots(1,1, figsize=(3,3))
    plot_singleValue('surface_energy', allData, num_bins, bins, bin_centers, ax, '-', color_line, color_shade)
    ax.set(xlabel='Time', ylabel='Surface Energy', title=f"Surface Energy, {dirPrefix.split('/')[-1]}")
    # ax.legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_surface_energy.pdf")
    plt.close()
    
    # plot surface area
    fig, ax = plt.subplots(1,1, figsize=(3,3))
    for traj in allData: 
        for frame in traj: 
            frame['VacXY_tot'] = np.sum(np.array(frame['vacNeighbors_counts']))
            frame['surfaceArea'] = frame['VacXY_tot'] * VacXY_hex_pixel_area / 100
    plot_singleValue('surfaceArea', allData, num_bins, bins, bin_centers, ax, '-', color_line, color_shade)
    ax.set(xlabel='Time', ylabel=r"Surface Area (${nm}^2$)", title=f"Surface Area, {dirPrefix.split('/')[-1]}", ylim=(0, 160))
    # ax.legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_surfaceArea.pdf")
    plt.close()
    
    '''
    # plot thickness over time
    fig, ax = plt.subplots(1,1, figsize=(3,3))
    plot_singleValue('thickness', allData, num_bins, bins, bin_centers, ax, '-', color_line, color_shade)
    ax.set(xlabel='Time', ylabel='Thickness (nm)', title=f"Thickness, {dirPrefix.split('/')[-1]}", ylim=(0, 42))
    # ax.legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_thickness.pdf")
    plt.close()

    # plot 3D number of neighbors
    fig, ax = plt.subplots(1,1, figsize=(3,3))
    plot_num_neighbors('neighbor_3D_counts', allData, num_bins, bins, bin_centers, ax)
    ax.set(xlabel='Time', ylabel='Number of atoms', title=f"3D num_neigh, {dirPrefix.split('/')[-1]}", ylim=(0, 3500))
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_site_neigh_count.pdf")
    plt.close()

    # plot 2D vacancy number of neighbors
    fig, ax = plt.subplots(1,1, figsize=(4,4))
    plot_num_neighbors('vacNeighbors_counts', allData, num_bins, bins, bin_centers, ax)
    # Do a running average and then extract the maximum. 
    ax.set(xlabel='Time', ylabel='Number of atoms', title=f"VacXY num_neigh, {dirPrefix.split('/')[-1]}", ylim=(0, 110))
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{writeToDir}plot_VacXY_neigh_count.pdf")
    plt.close()
    bin_sums = np.zeros(num_bins)
    bin_counts = np.zeros(num_bins)
    for traj in allData:
        times = np.array([item['simTime'] for item in traj])
        values = np.array([item['vacNeighbors_counts'][3] for item in traj])
        bin_indices = np.digitize(times, bins) - 1  # Find the bin index for each time point

        valid_indices = (0 <= bin_indices) & (bin_indices < num_bins)

        np.add.at(bin_sums, bin_indices[valid_indices], values[valid_indices])
        np.add.at(bin_counts, bin_indices[valid_indices], 1)
        
    non_zero_mask = bin_counts > 0
    bin_averages = np.divide(bin_sums[non_zero_mask], bin_counts[non_zero_mask])
    x = bin_centers[non_zero_mask]
    y = bin_averages
    window_size = 16
    y_running_avg = uniform_filter1d(y, size=window_size)
    max_index = np.argmax(y_running_avg)
    x_max = x[max_index]
    y_max = y_running_avg[max_index]
    print(f"3-neighbors maximum is at time = {(x_max):.5f}, the overall time = {(x[-1]):.5f}. The ratio is: {(x_max/x[-1]):.5f}")

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
    '''
    return


def uneven_plot_individual_trajectories(all_chemPot_sets, numTraj, cmap, dirPrefix):
    # This function create '_indivTraj.pdf' plots for each parameter set, where each trajectory is plotted individually. 
    for idx, chemPot in enumerate(all_chemPot_sets): 
        allData = []
        for repeat in range(numTraj):
            json_file_path = f"{dirPrefix}{chemPot}_repeat_{repeat+1}/stats.json"
            with open(json_file_path, 'r') as file:
                stats_list = json.load(file)
            allData.append(stats_list)

        #################################################################
        # plot number of cations, number of anions
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        cmap = plt.get_cmap('tab20')

        for trajIdx in range(numTraj): 
            times = np.array([item['simTime'] for item in allData[trajIdx]])
            values = np.array([item['num_cation'] for item in allData[trajIdx]])
            values2 = np.array([item['num_anion'] for item in allData[trajIdx]])
            if trajIdx == 0: 
                setAlpha = 1.0
            else: 
                setAlpha = 0.15
            ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[idx*2])
            ax.plot(times, values2, ':', alpha=setAlpha, color=cmap.colors[idx*2])
            # print(f"Initial num of atoms: {values[0]}, {np.max(values)}")
        ax.set(xlabel='Time', ylabel='Number of Atoms', ylim=(0, 12000))
        ax.grid(alpha=0.6)
        fig.suptitle(f"Num cations (or anions), {chemPot}, 16 traj")
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/plot_nAtoms_indivTraj.pdf")
        plt.close()

        #################################################################
        # plot surface area
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        for traj in allData: 
            for frame in traj: 
                frame['VacXY_tot'] = np.sum(np.array(frame['vacNeighbors_counts']))
                frame['surfaceArea'] = frame['VacXY_tot'] * VacXY_hex_pixel_area / 100
        
        for trajIdx in range(numTraj): 
            times = np.array([item['simTime'] for item in allData[trajIdx]])
            values = np.array([item['surfaceArea'] for item in allData[trajIdx]])
            if trajIdx == 0: 
                setAlpha = 1.0
            else: 
                setAlpha = 0.15
            ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[idx*2])
            # print(f"Initial surface area: {values[0]}, {np.max(values)}")
        ax.set(xlabel='Time', ylabel=r"Surface Area (${nm}^2$)", title=f"Surface Area, {chemPot}, 16 traj", ylim=(0, 160))
        ax.grid(alpha=0.6)
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/plot_surfaceArea_indivTraj.pdf")
        plt.close()

        #################################################################
        # plot energy
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        for trajIdx in range(numTraj): 
            times = np.array([item['simTime'] for item in allData[trajIdx]])
            values = np.array([item['energy'] for item in allData[trajIdx]])
            if trajIdx == 0: 
                setAlpha = 1.0
            else: 
                setAlpha = 0.15
            ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[idx*2])
        ax.set(xlabel='Time', ylabel='Total Energy', title=f"Total Energy, {chemPot}, 16 traj")
        ax.grid(alpha=0.6)
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/plot_energy_indivTraj.pdf")
        plt.close()
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        for trajIdx in range(numTraj): 
            times = np.array([item['simTime'] for item in allData[trajIdx]])
            values = np.array([item['binding_energy'] for item in allData[trajIdx]])
            if trajIdx == 0: 
                setAlpha = 1.0
            else: 
                setAlpha = 0.15
            ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[idx*2])
        ax.set(xlabel='Time', ylabel='Binding Energy', title=f"Binding Energy, {chemPot}, 16 traj")
        ax.grid(alpha=0.6)
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/plot_binding_energy_indivTraj.pdf")
        plt.close()
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        for trajIdx in range(numTraj): 
            times = np.array([item['simTime'] for item in allData[trajIdx]])
            values = np.array([item['surface_energy'] for item in allData[trajIdx]])
            if trajIdx == 0: 
                setAlpha = 1.0
            else: 
                setAlpha = 0.15
            ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[idx*2])
        ax.set(xlabel='Time', ylabel='Surface Energy', title=f"Surface Energy, {chemPot}, 16 traj")
        ax.grid(alpha=0.6)
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/plot_surface_energy_indivTraj.pdf")
        plt.close()

        '''
        #################################################################
        # plot thickness over time
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        for trajIdx in range(numTraj): 
            times = np.array([item['simTime'] for item in allData[trajIdx]])
            values = np.array([item['thickness'] for item in allData[trajIdx]])
            if trajIdx == 0: 
                setAlpha = 1.0
            else: 
                setAlpha = 0.15
            ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[idx*2])
        ax.set(xlabel='Time', ylabel='Thickness (nm)', title=f"Thickness, {chemPot}, 16 traj", ylim=(0, 42))
        ax.grid(alpha=0.6)
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/plot_thickness_indivTraj.pdf")
        plt.close()

        #################################################################
        # plot 3D number of neighbors
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        num_neighbors = len(allData[0][0]['neighbor_3D_counts']) - 1
        for trajIdx in range(numTraj): 
            for i in range(num_neighbors): 
                times = np.array([item['simTime'] for item in allData[trajIdx]])
                values = np.array([item['neighbor_3D_counts'][i] for item in allData[trajIdx]])
                if trajIdx == 0: 
                    setAlpha = 1.0
                    ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[(2*i)%20], label=f"{i}_neighbors")
                else: 
                    setAlpha = 0.15
                    ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[(2*i)%20])
        ax.set(xlabel='Time', ylabel='Number of atoms', title=f"3D num_neigh, {chemPot}, 16 traj", ylim=(0, 3500))
        ax.grid(alpha=0.6)
        ax.legend()
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/plot_site_neigh_count_indivTraj.pdf")
        plt.close()
        '''
        #################################################################
        # plot 2D vacancy number of neighbors
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        num_neighbors = len(allData[0][0]['vacNeighbors_counts']) - 1
        for trajIdx in range(numTraj): 
            for i in range(num_neighbors): 
                times = np.array([item['simTime'] for item in allData[trajIdx]])
                values = np.array([item['vacNeighbors_counts'][i] for item in allData[trajIdx]])
                if trajIdx == 0: 
                    setAlpha = 1.0
                    ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[(2*i)%20], label=f"{i}_neighbors")
                else: 
                    setAlpha = 0.15
                    ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[(2*i)%20])
        ax.set(xlabel='Time', ylabel='Number of atoms', title=f"VacXY num_neigh, {chemPot}, 16 traj", ylim=(0, 110))
        ax.grid(alpha=0.6)
        ax.legend()
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/plot_VacXY_neigh_count_indivTraj.pdf")
        plt.close()

        shutil.copy(f'{dirPrefix}{chemPot}/plot_energy_indivTraj.pdf', f"PLOTS_uneven/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_energy_indivTraj.pdf")
        shutil.copy(f'{dirPrefix}{chemPot}/plot_binding_energy_indivTraj.pdf', f"PLOTS_uneven/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_binding_energy_indivTraj.pdf")
        shutil.copy(f'{dirPrefix}{chemPot}/plot_surface_energy_indivTraj.pdf', f"PLOTS_uneven/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_surface_energy_indivTraj.pdf")
        shutil.copy(f'{dirPrefix}{chemPot}/plot_nAtoms_indivTraj.pdf', f"PLOTS_uneven/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_nAtoms_indivTraj.pdf")
        # shutil.copy(f'{dirPrefix}{chemPot}/plot_site_neigh_count_indivTraj.pdf', f"PLOTS_uneven/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_site_neigh_count_indivTraj.pdf")
        shutil.copy(f'{dirPrefix}{chemPot}/plot_surfaceArea_indivTraj.pdf', f"PLOTS_uneven/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_surfaceArea_indivTraj.pdf")
        # shutil.copy(f'{dirPrefix}{chemPot}/plot_thickness_indivTraj.pdf', f"PLOTS_uneven/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_thickness_indivTraj.pdf")
        shutil.copy(f'{dirPrefix}{chemPot}/plot_VacXY_neigh_count_indivTraj.pdf', f"PLOTS_uneven/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_VacXY_neigh_count_indivTraj.pdf")


def uneven_plot_all_datasets(all_chemPot_sets, deltaT_list, numTraj, dirPrefix, writeToDir, rescaleTime=False):
    # This function is for creating 1 plot for all parameter sets. 
    fig1, ax1 = plt.subplots(1,1, figsize=(5,3))
    fig2, ax2 = plt.subplots(1,1, figsize=(5,3))
    fig3, ax3 = plt.subplots(1,1, figsize=(5,3))
    fig4, ax4 = plt.subplots(1,1, figsize=(5,3))
    fig5, ax5 = plt.subplots(1,1, figsize=(5,3))
    fig6, ax6 = plt.subplots(1,1, figsize=(5,3))
    print("Summarizing all plots... ")

    cmap = plt.get_cmap('tab20')
    
    for idx, chemPot in enumerate(all_chemPot_sets): 
        min_simTime = np.inf
        allData = []
        simTimes = []
        for repeat in range(numTraj):
            json_file_path = f"{dirPrefix}{chemPot}_repeat_{repeat+1}/stats.json"

            with open(json_file_path, 'r') as file:
                stats_list = json.load(file)
            allData.append(stats_list)

            simTimes.append(stats_list[-1]['simTime'])

            if stats_list[-1]['simTime'] < min_simTime: 
                min_simTime = stats_list[-1]['simTime']

        num_bins = int(np.floor((min_simTime) / deltaT_list[idx]))
        bins = np.linspace(0.0, min_simTime, num_bins + 1)
        bin_centers = (bins[:-1] + bins[1:]) / 2

        simTime_std = np.std(simTimes)
        simTime_mean = np.mean(simTimes)
        print(f"{idx}, {chemPot} simTime mean and std: {simTime_mean:.6f}, {simTime_std:.6f}. Rel std = {(simTime_std/simTime_mean):.6f}")

        # plot number of cations, number of anions, fig1
        plot_singleValue('num_cation', allData, num_bins, bins, bin_centers, ax1, '-', cmap.colors[idx*2], cmap.colors[idx*2+1], label=f"{chemPot}", rescale=rescaleTime)
        plot_singleValue('num_anion', allData, num_bins, bins, bin_centers, ax1, ':', cmap.colors[idx*2], cmap.colors[idx*2+1], label=f"{chemPot}", rescale=rescaleTime)
        ax1.set(xlabel='Time', ylabel=f"Num cations (or anions), avg over {numTraj} traj", ylim=(0, 12000))
        if not rescaleTime: 
            ax1.set_xscale('log')
        ax1.legend()

        # plot surface area, fig4
        for traj in allData: 
            for frame in traj: 
                frame['VacXY_tot'] = np.sum(np.array(frame['vacNeighbors_counts']))
                frame['surfaceArea'] = frame['VacXY_tot'] * VacXY_hex_pixel_area / 100
        plot_singleValue('surfaceArea', allData, num_bins, bins, bin_centers, ax4, '-', cmap.colors[idx*2], cmap.colors[idx*2+1], label=f"{chemPot} ", rescale=rescaleTime)
        ax4.set(xlabel='Time', ylabel=r"Surface Area (${nm}^2$)", title=f"Surface Area, avg over {numTraj} traj", ylim=(0, 160))
        if not rescaleTime: 
            ax4.set_xscale('log')
        ax4.legend()

        # plot energy, fig2
        plot_singleValue('energy', allData, num_bins, bins, bin_centers, ax2, '-', cmap.colors[idx*2], cmap.colors[idx*2+1], label=f"{chemPot} ", rescale=rescaleTime)
        ax2.set(xlabel='Time', ylabel='Energy', title=f"Energy, avg over {numTraj} traj", ylim=(0, 50000))
        if not rescaleTime: 
            ax2.set_xscale('log')
        ax2.legend()
        plot_singleValue('binding_energy', allData, num_bins, bins, bin_centers, ax5, '-', cmap.colors[idx*2], cmap.colors[idx*2+1], label=f"{chemPot} ", rescale=rescaleTime)
        ax5.set(xlabel='Time', ylabel='Binding Energy', title=f"Binding Energy, avg over {numTraj} traj")
        if not rescaleTime: 
            ax5.set_xscale('log')
        ax5.legend()
        plot_singleValue('surface_energy', allData, num_bins, bins, bin_centers, ax6, '-', cmap.colors[idx*2], cmap.colors[idx*2+1], label=f"{chemPot} ", rescale=rescaleTime)
        ax6.set(xlabel='Time', ylabel='Surface Energy', title=f"Surface Energy, avg over {numTraj} traj")
        if not rescaleTime: 
            ax6.set_xscale('log')
        ax6.legend()
        
        # plot thickness over time, fig3
        plot_singleValue('thickness', allData, num_bins, bins, bin_centers, ax3, '-', cmap.colors[idx*2], cmap.colors[idx*2+1], label=f"{chemPot} ", rescale=rescaleTime)
        ax3.set(xlabel='Time', ylabel='Thickness (nm)', title=f"Thickness, avg over {numTraj} traj", ylim=(0, 42))
        if not rescaleTime: 
            ax3.set_xscale('log')
        ax3.legend()

    fig1.tight_layout()
    fig4.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
    fig5.tight_layout()
    fig6.tight_layout()
    if rescaleTime: 
        fig1.savefig(f"{writeToDir}uneven_plot_rescale_nAtoms.pdf")
        fig4.savefig(f"{writeToDir}uneven_plot_rescale_surfaceArea.pdf")
        fig2.savefig(f"{writeToDir}uneven_plot_rescale_energy.pdf")
        fig3.savefig(f"{writeToDir}uneven_plot_rescale_thickness.pdf")
        fig5.savefig(f"{writeToDir}uneven_plot_rescale_binding_energy.pdf")
        fig6.savefig(f"{writeToDir}uneven_plot_rescale_surface_energy.pdf")
    else: 
        fig1.savefig(f"{writeToDir}uneven_plot_nAtoms.pdf")
        fig4.savefig(f"{writeToDir}uneven_plot_surfaceArea.pdf")
        fig2.savefig(f"{writeToDir}uneven_plot_energy.pdf")
        fig3.savefig(f"{writeToDir}uneven_plot_thickness.pdf")
        fig5.savefig(f"{writeToDir}uneven_plot_binding_energy.pdf")
        fig6.savefig(f"{writeToDir}uneven_plot_surface_energy.pdf")
    return


######################################################
if __name__ == "__main__":
    ######################################################
    # Plot avg over multiple trajectories
    ######################################################
    all_chemPot_sets = ['-4.0_-3.5', '-4.0_-3.0', '-4.0_-2.5', '-4.0_-2.0', '-3.0_-2.5', '-3.0_-2.0']
    deltaT_list = [0.00008, 0.0003, 0.0005, 0.001, 0.01, 0.05]
    numTraj = 16
    cmap = plt.get_cmap('tab20')

    for idx, chemPot in enumerate(all_chemPot_sets): 
        uneven_avg_one_dataset(f'CALCS_uneven/{chemPot}', numTraj, f'CALCS_uneven/{chemPot}/', deltaT_list[idx], cmap.colors[idx*2], cmap.colors[idx*2+1])

        shutil.copy(f'CALCS_uneven/{chemPot}/plot_energy.pdf', f"PLOTS_uneven/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_energy.pdf")
        shutil.copy(f'CALCS_uneven/{chemPot}/plot_binding_energy.pdf', f"PLOTS_uneven/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_binding_energy.pdf")
        shutil.copy(f'CALCS_uneven/{chemPot}/plot_surface_energy.pdf', f"PLOTS_uneven/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_surface_energy.pdf")
        shutil.copy(f'CALCS_uneven/{chemPot}/plot_nAtoms.pdf', f"PLOTS_uneven/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_nAtoms.pdf")
        # shutil.copy(f'CALCS_uneven/{chemPot}/plot_site_neigh_count.pdf', f"PLOTS_uneven/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_site_neigh_count.pdf")
        shutil.copy(f'CALCS_uneven/{chemPot}/plot_surfaceArea.pdf', f"PLOTS_uneven/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_surfaceArea.pdf")
        # shutil.copy(f'CALCS_uneven/{chemPot}/plot_thickness.pdf', f"PLOTS_uneven/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_thickness.pdf")
        # shutil.copy(f'CALCS_uneven/{chemPot}/plot_VacXY_neigh_count.pdf', f"PLOTS_uneven/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_VacXY_neigh_count.pdf")
        # shutil.copy(f'CALCS_uneven/{chemPot}/plot_VacXY_theta.pdf', f"PLOTS_uneven/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_VacXY_theta.pdf")
        # shutil.copy(f'CALCS_uneven/{chemPot}/plot_VacXY_theta_reverse.pdf', f"PLOTS_uneven/{chemPot.replace('-', 'n').replace('.', 'p')}_plot_VacXY_theta_reverse.pdf")

    uneven_plot_all_datasets(all_chemPot_sets, deltaT_list, numTraj, 'CALCS_uneven/', 'PLOTS_uneven/', rescaleTime=False)
    uneven_plot_all_datasets(all_chemPot_sets, deltaT_list, numTraj, 'CALCS_uneven/', 'PLOTS_uneven/', rescaleTime=True)


    ######################################################
    # Plot individual trajectories
    ######################################################
    uneven_plot_individual_trajectories(all_chemPot_sets, numTraj, cmap, 'CALCS_uneven/')
