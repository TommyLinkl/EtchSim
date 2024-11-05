import shutil
import json, gzip
import numpy as np
from scipy.spatial import ConvexHull, QhullError
import matplotlib
import matplotlib.pyplot as plt
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
from plot_stats import plot_singleValue

def get_stats_from_contour(contour_filename):
    with gzip.open(contour_filename, 'rt') as file:
        contour_data = np.loadtxt(file)
    unique_simTimes = np.unique(contour_data[:, 0])  # sort

    contour_stats_list = []

    for stepNum, simTime in enumerate(unique_simTimes): 
        contour_coord = contour_data[contour_data[:, 0]==simTime][:, [1,2]]

        centroid = np.mean(contour_coord, axis=0)
        distances = np.linalg.norm(contour_coord - centroid, axis=1)
        VacXY_roughness_perimeter = np.std(distances)
        VacXY_rel_std = np.std(distances) / np.mean(distances)
        
        try:
            hull = ConvexHull(contour_coord)
            hull_perimeter = hull.area
            hull_area = hull.volume

            distances = np.sqrt(np.sum(np.diff(contour_coord, axis=0)**2, axis=1))
            actual_perimeter = np.sum(distances) + np.sqrt(np.sum((contour_coord[0] - contour_coord[-1])**2))
            # print(f"actual_perimeter = {actual_perimeter:.5f}")

            x = contour_coord[:, 0]
            y = contour_coord[:, 1]
            # Shoelace formula
            actual_area = 0.5 * np.abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
            # print(f"actual_area = {actual_area:.5f}")

            VacXY_roughness_hull_perimeter = actual_perimeter / hull_perimeter 
            VacXY_roughness_hull_area = actual_area / hull_area
            VacXY_cir_ratio = 4*np.pi*actual_area / (actual_perimeter)**2

        except QhullError as e:
            print(f"ConvexHull failed. (Expected behavior for end of simulation) We are setting VacXY_xxx to placeholder values. ")
            VacXY_roughness_hull_perimeter = 0.0
            VacXY_roughness_hull_area = 0.0
            VacXY_cir_ratio = 0.0

        contour_stats = {
            "stepNum": stepNum, 
            "simTime": simTime, 
            "VacXY_roughness_perimeter": VacXY_roughness_perimeter, 
            "VacXY_rel_std": VacXY_rel_std, 
            "VacXY_roughness_hull_perimeter": VacXY_roughness_hull_perimeter, 
            "VacXY_roughness_hull_area": VacXY_roughness_hull_area, 
            "VacXY_cir_ratio": VacXY_cir_ratio,
        }
        contour_stats_list.append(contour_stats)

    return contour_stats_list


def plot_individual_traj_contour_roughness(all_chemPot_sets, numTraj, cmap, dirPrefix):
    # This function create '_indivTraj.pdf' plots for each parameter set, where each trajectory is plotted individually. 
    for idx, chemPot in enumerate(all_chemPot_sets): 
        allData = []
        for repeat in range(numTraj):
            contour_stats_list = get_stats_from_contour(f"{dirPrefix}{chemPot}_repeat_{repeat+1}/vacXY_contour_pts.dat.gz")
            allData.append(contour_stats_list)
        cmap = plt.get_cmap('tab20')

        #################################################################
        # plot "VacXY_roughness_perimeter"
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        for trajIdx in range(numTraj): 
            times = np.array([item['simTime'] for item in allData[trajIdx]])
            values = np.array([item["VacXY_roughness_perimeter"] for item in allData[trajIdx]])
            if trajIdx == 0: 
                setAlpha = 1.0
            else: 
                setAlpha = 0.15
            ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[idx*2])
        ax.set(xlabel='Time', ylabel='Perimeter Roughness Index', title=f"Perimeter Roughness Index, {chemPot}, 16 traj", ylim=(0.5, 3.0))
        ax.grid(alpha=0.6)
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/contour_VacXY_roughness_perimeter_indivTraj.pdf")
        plt.close()
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        for trajIdx in range(numTraj): 
            times = np.array([item['simTime'] for item in allData[trajIdx]])
            values = np.array([item["VacXY_rel_std"] for item in allData[trajIdx]])
            if trajIdx == 0: 
                setAlpha = 1.0
            else: 
                setAlpha = 0.15
            ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[idx*2])
        ax.set(xlabel='Time', ylabel='Relative std.', title=f"Relative std. of distance, {chemPot}, 16 traj", ylim=(0, 0.2))
        ax.grid(alpha=0.6)
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/contour_VacXY_rel_std_indivTraj.pdf")
        plt.close()

        #################################################################
        # plot "VacXY_roughness_hull_perimeter"
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        for trajIdx in range(numTraj): 
            times = np.array([item['simTime'] for item in allData[trajIdx]])
            values = np.array([item["VacXY_roughness_hull_perimeter"] for item in allData[trajIdx]])
            if trajIdx == 0: 
                setAlpha = 1.0
            else: 
                setAlpha = 0.15
            ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[idx*2])
        ax.set(xlabel='Time', ylabel='VacXY_roughness_hull_perimeter', title=f"Convex Hull vs. Actual Perimeter, {chemPot}, 16 traj", ylim=(0.999, 1.008))
        ax.grid(alpha=0.6)
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/contour_VacXY_roughness_hull_perimeter_indivTraj.pdf")
        plt.close()

        #################################################################
        # plot "VacXY_roughness_hull_area"
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        for trajIdx in range(numTraj): 
            times = np.array([item['simTime'] for item in allData[trajIdx]])
            values = np.array([item["VacXY_roughness_hull_area"] for item in allData[trajIdx]])
            if trajIdx == 0: 
                setAlpha = 1.0
            else: 
                setAlpha = 0.15
            ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[idx*2])
        ax.set(xlabel='Time', ylabel='VacXY_roughness_hull_area', title=f"Convex Hull vs. Actual Area, {chemPot}, 16 traj", ylim=(0.97, 1.005))
        ax.grid(alpha=0.6)
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/contour_VacXY_roughness_hull_area_indivTraj.pdf")
        plt.close()

        #################################################################
        # plot "VacXY_cir_ratio"
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        for trajIdx in range(numTraj): 
            times = np.array([item['simTime'] for item in allData[trajIdx]])
            values = np.array([item["VacXY_cir_ratio"] for item in allData[trajIdx]])
            if trajIdx == 0: 
                setAlpha = 1.0
            else: 
                setAlpha = 0.15
            ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[idx*2])
        ax.set(xlabel='Time', ylabel='VacXY_cir_ratio', title=f"Circularity Ratio, {chemPot}, 16 traj", ylim=(0.9, 1.0))
        ax.grid(alpha=0.6)
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/contour_VacXY_cir_ratio_indivTraj.pdf")
        plt.close()

        shutil.copy(f'{dirPrefix}{chemPot}/contour_VacXY_roughness_perimeter_indivTraj.pdf', f"PLOTS_roughness/{chemPot.replace('-', 'n').replace('.', 'p')}_contour_VacXY_roughness_perimeter_indivTraj.pdf")
        shutil.copy(f'{dirPrefix}{chemPot}/contour_VacXY_rel_std_indivTraj.pdf', f"PLOTS_roughness/{chemPot.replace('-', 'n').replace('.', 'p')}_contour_VacXY_rel_std_indivTraj.pdf")
        shutil.copy(f'{dirPrefix}{chemPot}/contour_VacXY_roughness_hull_perimeter_indivTraj.pdf', f"PLOTS_roughness/{chemPot.replace('-', 'n').replace('.', 'p')}_contour_VacXY_roughness_hull_perimeter_indivTraj.pdf")
        shutil.copy(f'{dirPrefix}{chemPot}/contour_VacXY_roughness_hull_area_indivTraj.pdf', f"PLOTS_roughness/{chemPot.replace('-', 'n').replace('.', 'p')}_contour_VacXY_roughness_hull_area_indivTraj.pdf")
        shutil.copy(f'{dirPrefix}{chemPot}/contour_VacXY_cir_ratio_indivTraj.pdf', f"PLOTS_roughness/{chemPot.replace('-', 'n').replace('.', 'p')}_contour_VacXY_cir_ratio_indivTraj.pdf")


def plot_all_datasets_contour_roughness(all_chemPot_sets, deltaT_list, numTraj, dirPrefix, writeToDir, rescaleTime=False):
    # This function is for creating 1 plot for all parameter sets. 
    fig1, ax1 = plt.subplots(1,1, figsize=(5,3))
    fig2, ax2 = plt.subplots(1,1, figsize=(5,3))
    fig3, ax3 = plt.subplots(1,1, figsize=(5,3))
    fig4, ax4 = plt.subplots(1,1, figsize=(5,3))
    fig9, ax9 = plt.subplots(1,1, figsize=(5,3))
    print("Summarizing all plots... ")

    cmap = plt.get_cmap('tab20')
    
    for idx, chemPot in enumerate(all_chemPot_sets): 
        min_simTime = np.inf
        allData = []
        simTimes = []
        for repeat in range(numTraj):
            contour_stats_list = get_stats_from_contour(f"{dirPrefix}{chemPot}_repeat_{repeat+1}/vacXY_contour_pts.dat.gz")
            allData.append(contour_stats_list)
            simTimes.append(contour_stats_list[-1]['simTime'])

            if contour_stats_list[-1]['simTime'] < min_simTime: 
                min_simTime = contour_stats_list[-1]['simTime']

        num_bins = int(np.floor((min_simTime) / deltaT_list[idx]))
        bins = np.linspace(0.0, min_simTime, num_bins + 1)
        bin_centers = (bins[:-1] + bins[1:]) / 2

        # simTime_std = np.std(simTimes)
        # simTime_mean = np.mean(simTimes)
        # print(f"{idx}, {chemPot} simTime mean and std: {simTime_mean:.6f}, {simTime_std:.6f}. Rel std = {(simTime_std/simTime_mean):.6f}")

        plot_singleValue("VacXY_roughness_perimeter", allData, num_bins, bins, bin_centers, ax1, '-', cmap.colors[idx*2], cmap.colors[idx*2+1], label=f"{chemPot} ", rescale=rescaleTime)
        ax1.set(xlabel='Time', ylabel=r"Perimeter Roughness Index", title=f"Perimeter Roughness Index, avg over {numTraj} traj", ylim=(0.5, 3.0))
        if not rescaleTime: 
            ax1.set_xscale('log')
        ax1.legend(fontsize=7)
        plot_singleValue("VacXY_rel_std", allData, num_bins, bins, bin_centers, ax9, '-', cmap.colors[idx*2], cmap.colors[idx*2+1], label=f"{chemPot} ", rescale=rescaleTime)
        ax9.set(xlabel='Time', ylabel='Relative std.', title=f"Relative std. of distance, avg over {numTraj} traj", ylim=(0, 0.2))
        if not rescaleTime: 
            ax9.set_xscale('log')
        ax9.legend(fontsize=7)

        plot_singleValue("VacXY_roughness_hull_perimeter", allData, num_bins, bins, bin_centers, ax2, '-', cmap.colors[idx*2], cmap.colors[idx*2+1], label=f"{chemPot} ", rescale=rescaleTime)
        ax2.set(xlabel='Time', ylabel='VacXY_roughness_hull_perimeter', title=f"Convex Hull vs. Actual Perimeter, avg over {numTraj} traj", ylim=(0.999, 1.008))
        if not rescaleTime: 
            ax2.set_xscale('log')
        ax2.legend(fontsize=7)
        
        plot_singleValue('VacXY_roughness_hull_area', allData, num_bins, bins, bin_centers, ax3, '-', cmap.colors[idx*2], cmap.colors[idx*2+1], label=f"{chemPot} ", rescale=rescaleTime)
        ax3.set(xlabel='Time', ylabel='VacXY_roughness_hull_area', title=f"Convex Hull vs. Actual Area, avg over {numTraj} traj", ylim=(0.97, 1.005))
        if not rescaleTime: 
            ax3.set_xscale('log')
        ax3.legend(fontsize=7)

        plot_singleValue('VacXY_cir_ratio', allData, num_bins, bins, bin_centers, ax4, '-', cmap.colors[idx*2], cmap.colors[idx*2+1], label=f"{chemPot} ", rescale=rescaleTime)
        ax4.set(xlabel='Time', ylabel='VacXY_cir_ratio', title=f"Circularity Ratio, avg over {numTraj} traj", ylim=(0.9, 1.0))
        if not rescaleTime: 
            ax4.set_xscale('log')
        ax4.legend(fontsize=7)

    fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
    fig4.tight_layout()
    fig9.tight_layout()
    if rescaleTime: 
        fig1.savefig(f"{writeToDir}plot_contour_rescale_VacXY_roughness_perimeter.pdf")
        fig2.savefig(f"{writeToDir}plot_contour_rescale_VacXY_roughness_hull_perimeter.pdf")
        fig3.savefig(f"{writeToDir}plot_contour_rescale_VacXY_roughness_hull_area.pdf")
        fig4.savefig(f"{writeToDir}plot_contour_rescale_cir_ratio.pdf")
        fig9.savefig(f"{writeToDir}plot_contour_rescale_VacXY_rel_std.pdf")
    else: 
        fig1.savefig(f"{writeToDir}plot_contour_VacXY_roughness_perimeter.pdf")
        fig2.savefig(f"{writeToDir}plot_contour_VacXY_roughness_hull_perimeter.pdf")
        fig3.savefig(f"{writeToDir}plot_contour_VacXY_roughness_hull_area.pdf")
        fig4.savefig(f"{writeToDir}plot_contour_cir_ratio.pdf")
        fig9.savefig(f"{writeToDir}plot_contour_VacXY_rel_std.pdf")
    return


######################################################
if __name__ == "__main__":
    # all_chemPot_sets = ['-4.0_-4.0', '-3.5_-3.5', '-3.0_-3.0', '-2.5_-2.5', '-2.0_-2.0'] 
    # deltaT_list = [0.00003, 0.0003, 0.003, 0.05, 2.8]

    all_chemPot_sets = ['-2.5_-2.0', '-3.0_-2.0', '-3.5_-2.0', '-4.0_-2.0'] 
    deltaT_list = [0.4, 0.1, 0.015, 0.002]

    numTraj = 16
    cmap = plt.get_cmap('tab20')

    plot_all_datasets_contour_roughness(all_chemPot_sets, deltaT_list, numTraj, 'CALCS_uneven/', 'PLOTS_roughness_uneven/', rescaleTime=False)
    plot_all_datasets_contour_roughness(all_chemPot_sets, deltaT_list, numTraj, 'CALCS_uneven/', 'PLOTS_roughness_uneven/', rescaleTime=True)

    ######################################################
    # Plot individual trajectories
    ######################################################
    plot_individual_traj_contour_roughness(all_chemPot_sets, numTraj, cmap, 'CALCS_uneven/')