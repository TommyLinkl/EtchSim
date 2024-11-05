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

def get_etching_time_data(chemPot1_list, chemPot2_list, numTraj, dirPrefix):
    plot_data_mean = []
    plot_data_std = []
    
    for (chemPot1, chemPot2) in zip(chemPot1_list, chemPot2_list): 
        chemPot_prefix = f"{chemPot1:.1f}_{chemPot2:.1f}"
        print(chemPot_prefix)

        simTimes = []
        for repeat in range(numTraj):
            json_file_path = f"{dirPrefix}{chemPot_prefix}_repeat_{repeat+1}/stats.json"

            with open(json_file_path, 'r') as file:
                stats_list = json.load(file)

            simTimes.append(stats_list[-1]['simTime'])

        simTime_mean = np.mean(simTimes)
        simTime_std = np.std(simTimes)
        plot_data_mean.append(simTime_mean)
        plot_data_std.append(simTime_std)
        # print(f"{chemPot_prefix} simTime mean and std: {simTime_mean:.1e}, {simTime_std:.1e}. Rel std = {(simTime_std/simTime_mean):.1e}")

    return plot_data_mean, plot_data_std


######################################################
if __name__ == "__main__":
    fig, axs = plt.subplots(1,2, figsize=(8,4))

    chemPot1_list = [-4.0, -3.5, -3.0, -2.5, -2.0]
    chemPot2_list = [-4.0, -3.5, -3.0, -2.5, -2.0]
    numTraj = 16
    eq_mean, eq_std = get_etching_time_data(chemPot1_list, chemPot2_list, numTraj, 'CALCS_diam15.4_thick3.5/')
    axs[0].errorbar(chemPot1_list, eq_mean, yerr=eq_std, fmt='o-', capsize=5, capthick=2, markersize=6, label=r"Equal $\mu$")
    eq_mean, eq_std = get_etching_time_data([-4.0], [-2.0], numTraj, 'CALCS_uneven/')
    axs[0].errorbar([-3.0], eq_mean, yerr=eq_std, fmt='o-', capsize=5, capthick=2, markersize=6, alpha=0.7, label="-4.0_-2.0")
    eq_mean, eq_std = get_etching_time_data([-4.0], [-2.5], numTraj, 'CALCS_uneven/')
    axs[0].errorbar([-3.25], eq_mean, yerr=eq_std, fmt='o-', capsize=5, capthick=2, markersize=6, alpha=0.7, label="-4.0_-2.5")
    eq_mean, eq_std = get_etching_time_data([-4.0], [-3.0], numTraj, 'CALCS_uneven/')
    axs[0].errorbar([-3.5], eq_mean, yerr=eq_std, fmt='o-', capsize=5, capthick=2, markersize=6, alpha=0.7, label="-4.0_-3.0")
    eq_mean, eq_std = get_etching_time_data([-4.0], [-3.5], numTraj, 'CALCS_uneven/')
    axs[0].errorbar([-3.75], eq_mean, yerr=eq_std, fmt='o-', capsize=5, capthick=2, markersize=6, alpha=0.7, label="-4.0_-3.5")
    
    eq_mean, eq_std = get_etching_time_data([-3.5], [-2.0], numTraj, 'CALCS_uneven/')
    axs[0].errorbar([-2.75], eq_mean, yerr=eq_std, fmt='o-', capsize=5, capthick=2, markersize=6, alpha=0.7, label="-3.5_-2.0")
    eq_mean, eq_std = get_etching_time_data([-3.5], [-2.5], numTraj, 'CALCS_uneven/')
    axs[0].errorbar([-3], eq_mean, yerr=eq_std, fmt='o-', capsize=5, capthick=2, markersize=6, alpha=0.7, label="-3.5_-2.5")
    eq_mean, eq_std = get_etching_time_data([-3.5], [-3.0], numTraj, 'CALCS_uneven/')
    axs[0].errorbar([-3.25], eq_mean, yerr=eq_std, fmt='o-', capsize=5, capthick=2, markersize=6, alpha=0.7, label="-3.5_-3.0")

    eq_mean, eq_std = get_etching_time_data([-3.0], [-2.0], numTraj, 'CALCS_uneven/')
    axs[0].errorbar([-2.5], eq_mean, yerr=eq_std, fmt='o-', capsize=5, capthick=2, markersize=6, alpha=0.7, label="-3.0_-2.0")
    eq_mean, eq_std = get_etching_time_data([-3.0], [-2.5], numTraj, 'CALCS_uneven/')
    axs[0].errorbar([-2.75], eq_mean, yerr=eq_std, fmt='o-', capsize=5, capthick=2, markersize=6, alpha=0.7, label="-3.0_-2.5")

    eq_mean, eq_std = get_etching_time_data([-2.5], [-2.0], numTraj, 'CALCS_uneven/')
    axs[0].errorbar([-2.25], eq_mean, yerr=eq_std, fmt='o-', capsize=5, capthick=2, markersize=6, alpha=0.7, label="-2.5_-2.0")

    axs[0].set_yscale("log")
    axs[0].grid(alpha=0.3)
    axs[0].legend(fontsize=8)
    axs[0].set(ylabel="Etching time (arb. unit)", xlabel="Average chemical potential of In and P")

    chemPot1_list = [-3.0, -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2.0]
    chemPot2_list = [-3.0, -2.9, -2.8, -2.7, -2.6, -2.5, -2.4, -2.3, -2.2, -2.1, -2.0]
    numTraj = 8
    eq_mean, eq_std = get_etching_time_data(chemPot1_list, chemPot2_list, numTraj, 'CALCS_lowTemp/')
    axs[1].errorbar(chemPot1_list, eq_mean, yerr=eq_std, fmt='o-', capsize=5, capthick=2, markersize=6, label=r"Equal $\mu$")
    axs[1].set_yscale("log")
    axs[1].grid(alpha=0.3)
    axs[1].legend()
    axs[1].set(ylabel="Etching time (arb. unit)", xlabel="Chemical potentials", title="Lower Temp, tests")

    fig.tight_layout()
    fig.savefig(f"plot_etching_time.pdf")