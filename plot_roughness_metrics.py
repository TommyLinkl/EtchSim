import shutil
import json
import numpy as np
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
from plot_stats import plot_singleValue

def plot_individual_traj_roughness(all_chemPot_sets, numTraj, cmap, dirPrefix):
    # This function create '_indivTraj.pdf' plots for each parameter set, where each trajectory is plotted individually. 
    for idx, chemPot in enumerate(all_chemPot_sets): 
        allData = []
        for repeat in range(numTraj):
            json_file_path = f"{dirPrefix}{chemPot}_repeat_{repeat+1}/stats.json"
            with open(json_file_path, 'r') as file:
                stats_list = json.load(file)
            allData.append(stats_list)
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
            ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[(idx*2) % 20])
        ax.set(xlabel='Time', ylabel='Perimeter Roughness Index', title=f"Perimeter Roughness Index, {chemPot}, 16 traj", ylim=(1, 4))
        ax.grid(alpha=0.6)
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/VacXY_roughness_perimeter_indivTraj.pdf")
        plt.close()
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        for trajIdx in range(numTraj): 
            times = np.array([item['simTime'] for item in allData[trajIdx]])
            values = np.array([item["VacXY_rel_std"] for item in allData[trajIdx]])
            if trajIdx == 0: 
                setAlpha = 1.0
            else: 
                setAlpha = 0.15
            ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[(idx*2) % 20])
        ax.set(xlabel='Time', ylabel='Relative std.', title=f"Relative std. of distance, {chemPot}, 16 traj", ylim=(0, 0.4))
        ax.grid(alpha=0.6)
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/VacXY_rel_std_indivTraj.pdf")
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
            ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[(idx*2) % 20])
        ax.set(xlabel='Time', ylabel='VacXY_roughness_hull_perimeter', title=f"Convex Hull vs. Actual Perimeter, {chemPot}, 16 traj", ylim=(0.95, 1.2))
        ax.grid(alpha=0.6)
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/VacXY_roughness_hull_perimeter_indivTraj.pdf")
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
            ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[(idx*2) % 20])
        ax.set(xlabel='Time', ylabel='VacXY_roughness_hull_area', title=f"Convex Hull vs. Actual Area, {chemPot}, 16 traj", ylim=(0.7, 1.05))
        ax.grid(alpha=0.6)
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/VacXY_roughness_hull_area_indivTraj.pdf")
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
            ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[(idx*2) % 20])
        ax.set(xlabel='Time', ylabel='VacXY_cir_ratio', title=f"Circularity Ratio, {chemPot}, 16 traj", ylim=(0.5, 1.1))
        ax.grid(alpha=0.6)
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/VacXY_cir_ratio_indivTraj.pdf")
        plt.close()

        #################################################################
        # plot "XZ_bb_AR"
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        for trajIdx in range(numTraj): 
            times = np.array([item['simTime'] for item in allData[trajIdx]])
            values = np.array([item["XZ_bb_AR"] for item in allData[trajIdx]])
            if trajIdx == 0: 
                setAlpha = 1.0
            else: 
                setAlpha = 0.15
            ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[(idx*2) % 20])
        ax.set(xlabel='Time', ylabel='XZ_bb_AR', title=f"XZ proj, bb AR, {chemPot}, 16 traj", ylim=(1.0, 4.5))
        ax.grid(alpha=0.6)
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/XZ_bb_AR_indivTraj.pdf")
        plt.close()

        #################################################################
        # plot "XZ_area_ratio"
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        for trajIdx in range(numTraj): 
            times = np.array([item['simTime'] for item in allData[trajIdx]])
            values = np.array([item["XZ_area_ratio"] for item in allData[trajIdx]])
            if trajIdx == 0: 
                setAlpha = 1.0
            else: 
                setAlpha = 0.15
            ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[(idx*2) % 20])
        ax.set(xlabel='Time', ylabel='XZ_area_ratio', title=f"XZ proj, hull area ratio, {chemPot}, 16 traj", ylim=(0.70, 1.05))
        ax.grid(alpha=0.6)
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/XZ_area_ratio_indivTraj.pdf")
        plt.close()

        #################################################################
        # plot "YZ_bb_AR"
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        for trajIdx in range(numTraj): 
            times = np.array([item['simTime'] for item in allData[trajIdx]])
            values = np.array([item["YZ_bb_AR"] for item in allData[trajIdx]])
            if trajIdx == 0: 
                setAlpha = 1.0
            else: 
                setAlpha = 0.15
            ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[(idx*2) % 20])
        ax.set(xlabel='Time', ylabel='YZ_bb_AR', title=f"YZ proj, bb AR, {chemPot}, 16 traj", ylim=(1.0, 4.5))
        ax.grid(alpha=0.6)
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/YZ_bb_AR_indivTraj.pdf")
        plt.close()

        #################################################################
        # plot "YZ_area_ratio"
        fig, ax = plt.subplots(1,1, figsize=(4,4))
        for trajIdx in range(numTraj): 
            times = np.array([item['simTime'] for item in allData[trajIdx]])
            values = np.array([item["YZ_area_ratio"] for item in allData[trajIdx]])
            if trajIdx == 0: 
                setAlpha = 1.0
            else: 
                setAlpha = 0.15
            ax.plot(times, values, '-', alpha=setAlpha, color=cmap.colors[(idx*2) % 20])
        ax.set(xlabel='Time', ylabel='YZ_area_ratio', title=f"YZ proj, hull area ratio, {chemPot}, 16 traj", ylim=(0.70, 1.05))
        ax.grid(alpha=0.6)
        fig.tight_layout()
        fig.savefig(f"{dirPrefix}{chemPot}/YZ_area_ratio_indivTraj.pdf")
        plt.close()

        shutil.copy(f'{dirPrefix}{chemPot}/VacXY_roughness_perimeter_indivTraj.pdf', f"PLOTS_roughness/{chemPot.replace('-', 'n').replace('.', 'p')}_VacXY_roughness_perimeter_indivTraj.pdf")
        shutil.copy(f'{dirPrefix}{chemPot}/VacXY_rel_std_indivTraj.pdf', f"PLOTS_roughness/{chemPot.replace('-', 'n').replace('.', 'p')}_VacXY_rel_std_indivTraj.pdf")
        shutil.copy(f'{dirPrefix}{chemPot}/VacXY_roughness_hull_perimeter_indivTraj.pdf', f"PLOTS_roughness/{chemPot.replace('-', 'n').replace('.', 'p')}_VacXY_roughness_hull_perimeter_indivTraj.pdf")
        shutil.copy(f'{dirPrefix}{chemPot}/VacXY_roughness_hull_area_indivTraj.pdf', f"PLOTS_roughness/{chemPot.replace('-', 'n').replace('.', 'p')}_VacXY_roughness_hull_area_indivTraj.pdf")
        shutil.copy(f'{dirPrefix}{chemPot}/VacXY_cir_ratio_indivTraj.pdf', f"PLOTS_roughness/{chemPot.replace('-', 'n').replace('.', 'p')}_VacXY_cir_ratio_indivTraj.pdf")
        shutil.copy(f'{dirPrefix}{chemPot}/XZ_bb_AR_indivTraj.pdf', f"PLOTS_roughness/{chemPot.replace('-', 'n').replace('.', 'p')}_XZ_bb_AR_indivTraj.pdf")
        shutil.copy(f'{dirPrefix}{chemPot}/XZ_area_ratio_indivTraj.pdf', f"PLOTS_roughness/{chemPot.replace('-', 'n').replace('.', 'p')}_XZ_area_ratio_indivTraj.pdf")
        shutil.copy(f'{dirPrefix}{chemPot}/YZ_bb_AR_indivTraj.pdf', f"PLOTS_roughness/{chemPot.replace('-', 'n').replace('.', 'p')}_YZ_bb_AR_indivTraj.pdf")
        shutil.copy(f'{dirPrefix}{chemPot}/YZ_area_ratio_indivTraj.pdf', f"PLOTS_roughness/{chemPot.replace('-', 'n').replace('.', 'p')}_YZ_area_ratio_indivTraj.pdf")


def plot_all_datasets_roughness(all_chemPot_sets, deltaT_list, numTraj, dirPrefix, writeToDir, rescaleTime=False):
    # This function is for creating 1 plot for all parameter sets. 
    fig1, ax1 = plt.subplots(1,1, figsize=(5,3))
    fig2, ax2 = plt.subplots(1,1, figsize=(5,3))
    fig3, ax3 = plt.subplots(1,1, figsize=(5,3))
    fig4, ax4 = plt.subplots(1,1, figsize=(5,3))
    fig5, ax5 = plt.subplots(1,1, figsize=(5,3))
    fig6, ax6 = plt.subplots(1,1, figsize=(5,3))
    fig7, ax7 = plt.subplots(1,1, figsize=(5,3))
    fig8, ax8 = plt.subplots(1,1, figsize=(5,3))
    fig9, ax9 = plt.subplots(1,1, figsize=(5,3))
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
        print(f"{dirPrefix}\nThe shortest simulation time is {min_simTime:.3e}. And we have split it into {num_bins} bins for statistics plotting.")

        simTime_std = np.std(simTimes)
        simTime_mean = np.mean(simTimes)
        print(f"{idx}, {chemPot} simTime mean and std: {simTime_mean:.1e}, {simTime_std:.1e}. Rel std = {(simTime_std/simTime_mean):.1e}")

        plot_singleValue("VacXY_roughness_perimeter", allData, num_bins, bins, bin_centers, ax1, '-', cmap.colors[(idx*2) % 20], cmap.colors[(idx*2+1) % 20], label=f"{chemPot} ", rescale=rescaleTime)
        ax1.set(xlabel='Time', ylabel=r"Perimeter Roughness Index", title=f"Perimeter Roughness Index, avg over {numTraj} traj", ylim=(1, 4))
        if not rescaleTime: 
            ax1.set_xscale('log')
        ax1.legend(fontsize=7)
        plot_singleValue("VacXY_rel_std", allData, num_bins, bins, bin_centers, ax9, '-', cmap.colors[(idx*2) % 20], cmap.colors[(idx*2+1) % 20], label=f"{chemPot} ", rescale=rescaleTime)
        ax9.set(xlabel='Time', ylabel='Relative std.', title=f"Relative std. of distance, avg over {numTraj} traj", ylim=(0, 0.4))
        if not rescaleTime: 
            ax9.set_xscale('log')
        ax9.legend(fontsize=7)

        plot_singleValue("VacXY_roughness_hull_perimeter", allData, num_bins, bins, bin_centers, ax2, '-', cmap.colors[(idx*2) % 20], cmap.colors[(idx*2+1) % 20], label=f"{chemPot} ", rescale=rescaleTime)
        ax2.set(xlabel='Time', ylabel='VacXY_roughness_hull_perimeter', title=f"Convex Hull vs. Actual Perimeter, avg over {numTraj} traj", ylim=(0.95, 1.2))
        if not rescaleTime: 
            ax2.set_xscale('log')
        ax2.legend(fontsize=7)
        
        plot_singleValue('VacXY_roughness_hull_area', allData, num_bins, bins, bin_centers, ax3, '-', cmap.colors[(idx*2) % 20], cmap.colors[(idx*2+1) % 20], label=f"{chemPot} ", rescale=rescaleTime)
        ax3.set(xlabel='Time', ylabel='VacXY_roughness_hull_area', title=f"Convex Hull vs. Actual Area, avg over {numTraj} traj", ylim=(0.7, 1.05))
        if not rescaleTime: 
            ax3.set_xscale('log')
        ax3.legend(fontsize=7)

        plot_singleValue('VacXY_cir_ratio', allData, num_bins, bins, bin_centers, ax4, '-', cmap.colors[(idx*2) % 20], cmap.colors[(idx*2+1) % 20], label=f"{chemPot} ", rescale=rescaleTime)
        ax4.set(xlabel='Time', ylabel='VacXY_cir_ratio', title=f"Circularity Ratio, avg over {numTraj} traj", ylim=(0.5, 1.1))
        if not rescaleTime: 
            ax4.set_xscale('log')
        ax4.legend(fontsize=7)

        plot_singleValue('XZ_bb_AR', allData, num_bins, bins, bin_centers, ax5, '-', cmap.colors[(idx*2) % 20], cmap.colors[(idx*2+1) % 20], label=f"{chemPot} ", rescale=rescaleTime)
        ax5.set(xlabel='Time', ylabel='XZ_bb_AR', title=f"XZ proj, bb AR, avg over {numTraj} traj", ylim=(1.0, 4.5))
        if not rescaleTime: 
            ax5.set_xscale('log')
        ax5.legend(fontsize=7)

        plot_singleValue('XZ_area_ratio', allData, num_bins, bins, bin_centers, ax6, '-', cmap.colors[(idx*2) % 20], cmap.colors[(idx*2+1) % 20], label=f"{chemPot} ", rescale=rescaleTime)
        ax6.set(xlabel='Time', ylabel='XZ_area_ratio', title=f"XZ proj, hull area ratio, avg over {numTraj} traj", ylim=(0.70, 1.05))
        if not rescaleTime: 
            ax6.set_xscale('log')
        ax6.legend(fontsize=7)

        plot_singleValue('YZ_bb_AR', allData, num_bins, bins, bin_centers, ax7, '-', cmap.colors[(idx*2) % 20], cmap.colors[(idx*2+1) % 20], label=f"{chemPot} ", rescale=rescaleTime)
        ax7.set(xlabel='Time', ylabel='YZ_bb_AR', title=f"YZ proj, bb AR, avg over {numTraj} traj", ylim=(1.0, 4.5))
        if not rescaleTime: 
            ax7.set_xscale('log')
        ax7.legend(fontsize=7)

        plot_singleValue('YZ_area_ratio', allData, num_bins, bins, bin_centers, ax8, '-', cmap.colors[(idx*2) % 20], cmap.colors[(idx*2+1) % 20], label=f"{chemPot} ", rescale=rescaleTime)
        ax8.set(xlabel='Time', ylabel='YZ_area_ratio', title=f"YZ proj, hull area ratio, avg over {numTraj} traj", ylim=(0.70, 1.05))
        if not rescaleTime: 
            ax8.set_xscale('log')
        ax8.legend(fontsize=7)

    fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
    fig4.tight_layout()
    fig5.tight_layout()
    fig6.tight_layout()
    fig7.tight_layout()
    fig8.tight_layout()
    fig9.tight_layout()
    if rescaleTime: 
        fig1.savefig(f"{writeToDir}plot_rescale_VacXY_roughness_perimeter.pdf")
        fig2.savefig(f"{writeToDir}plot_rescale_VacXY_roughness_hull_perimeter.pdf")
        fig3.savefig(f"{writeToDir}plot_rescale_VacXY_roughness_hull_area.pdf")
        fig4.savefig(f"{writeToDir}plot_rescale_cir_ratio.pdf")
        fig5.savefig(f"{writeToDir}plot_rescale_XZ_bb_AR.pdf")
        fig6.savefig(f"{writeToDir}plot_rescale_XZ_area_ratio.pdf")
        fig7.savefig(f"{writeToDir}plot_rescale_YZ_bb_AR.pdf")
        fig8.savefig(f"{writeToDir}plot_rescale_YZ_area_ratio.pdf")
        fig9.savefig(f"{writeToDir}plot_rescale_VacXY_rel_std.pdf")
    else: 
        fig1.savefig(f"{writeToDir}plot_VacXY_roughness_perimeter.pdf")
        fig2.savefig(f"{writeToDir}plot_VacXY_roughness_hull_perimeter.pdf")
        fig3.savefig(f"{writeToDir}plot_VacXY_roughness_hull_area.pdf")
        fig4.savefig(f"{writeToDir}plot_cir_ratio.pdf")
        fig5.savefig(f"{writeToDir}plot_XZ_bb_AR.pdf")
        fig6.savefig(f"{writeToDir}plot_XZ_area_ratio.pdf")
        fig7.savefig(f"{writeToDir}plot_YZ_bb_AR.pdf")
        fig8.savefig(f"{writeToDir}plot_YZ_area_ratio.pdf")
        fig9.savefig(f"{writeToDir}plot_VacXY_rel_std.pdf")
    return


######################################################
if __name__ == "__main__":
    # all_chemPot_sets = ['-4.0_-4.0', '-3.5_-3.5', '-3.0_-3.0', '-2.5_-2.5', '-2.0_-2.0'] 
    # deltaT_list = [0.00003, 0.0003, 0.003, 0.05, 2.8]

    all_chemPot_sets = ['-4.0_-2.0', '-4.0_-2.5', '-4.0_-3.0', '-4.0_-3.5'] 
    deltaT_list = [0.002, 0.001, 0.0005, 0.00015]

    all_chemPot_sets = ['-2.5_-2.0', '-3.0_-2.0', '-3.5_-2.0', '-4.0_-2.0'] 
    deltaT_list = [0.4, 0.1, 0.015, 0.002]

    # all_chemPot_sets = ['-2.0_-2.0', '-2.2_-2.2', '-2.4_-2.4', '-2.6_-2.6', '-2.8_-2.8', '-3.0_-3.0'] 
    # deltaT_list = [4.0, 3e-5, 1.5e-9, 7e-14, 4e-18, 1e-22]

    numTraj = 16
    cmap = plt.get_cmap('tab20')

    plot_all_datasets_roughness(all_chemPot_sets, deltaT_list, numTraj, 'CALCS_uneven/', 'PLOTS_roughness_uneven/', rescaleTime=False)
    plot_all_datasets_roughness(all_chemPot_sets, deltaT_list, numTraj, 'CALCS_uneven/', 'PLOTS_roughness_uneven/', rescaleTime=True)

    ######################################################
    # Plot individual trajectories
    ######################################################
    plot_individual_traj_roughness(all_chemPot_sets, numTraj, cmap, 'CALCS_uneven/')