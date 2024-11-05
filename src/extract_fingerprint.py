import numpy as np
import shutil
import gzip
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.interpolate import CubicSpline
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

from .extract_contour import reconstruct_boundary_pts_over_time


def clockwise_angle_from_x_axis(vector):
    norm = np.linalg.norm(vector)
    
    if norm == 0:
        raise ValueError("Zero vector has no direction.")
    
    angle = np.arctan2(vector[1], vector[0])
    clockwise_angle = -angle % (2 * np.pi)  # Ensure it is in [0, 2π)
    
    return clockwise_angle


def calc_fingerprint_function(calc_dir):
    # Read in file "vacXY_boundary_pts.dat.gz"
    boundary_pts_list, unique_times = reconstruct_boundary_pts_over_time(f"{calc_dir}vacXY_boundary_pts.dat.gz")

    # calculate the fingerprint function and write to compressed file
    fingerprint_zip = f"{calc_dir}fingerprint.dat.gz"
    with gzip.open(fingerprint_zip, 'wt') as f:
        f.write("# time      theta       r(theta, t)\n")
    
    for frame in boundary_pts_list:
        simTime = np.unique(frame[:, 0])[0]
        boundary_coordXY = frame[:, 1:]

        if (boundary_coordXY.size > 4) and (boundary_coordXY.ndim > 1):
            centroid = np.mean(boundary_coordXY, axis=0)
        else:
            continue

        # Record r(theta, t)
        with gzip.open(fingerprint_zip, 'at') as file:
            for point in boundary_coordXY: 
                theta = clockwise_angle_from_x_axis(point-centroid)
                r = np.linalg.norm(point - centroid)
                file.write(f"{simTime}     {theta}      {r}\n")


def plot_fingerprint_function_and_FT(fp_filename, fp_fig_name, fp_FT_fig_name, every_frame=20, end_before=-40):
    fig_fp, ax_fp = plt.subplots(1, 1, figsize=(7, 4))
    fig_fpFT, axs_fpFT = plt.subplots(4, 1, figsize=(6, 8))

    with gzip.open(fp_filename, 'rt') as f: 
        fingerprint_data = np.loadtxt(f)

    unique_times = np.unique(fingerprint_data[:, 0])
    norm = plt.Normalize(vmin=fingerprint_data[:, 0].min(), vmax=fingerprint_data[:, 0].max())
    colors = plt.cm.winter(norm(unique_times))

    for time, color in zip(unique_times[:end_before:every_frame], colors[:end_before:every_frame]):
        mask = fingerprint_data[:, 0] == time
        theta_values = fingerprint_data[mask, 1]
        r_values = fingerprint_data[mask, 2]
        
        # Sort the theta and corresponding r values
        sorted_indices = np.argsort(theta_values)
        theta_values_sorted = theta_values[sorted_indices]
        r_values_sorted = r_values[sorted_indices]
        # ax.plot(theta_values_sorted, r_values_sorted, '-', color=color, label=f'Time: {time}')

        cubic_spline = CubicSpline(theta_values_sorted, r_values_sorted)
        theta_fine = np.linspace(theta_values_sorted.min(), theta_values_sorted.max(), 500)
        r_fine = cubic_spline(theta_fine)
        ax_fp.plot(theta_fine, r_fine, '-', color=color, label=f'Time: {time}')


        # FT
        N = len(theta_fine)
        c_n = np.fft.fft(r_fine) / N
        frequencies = np.fft.fftfreq(N, (theta_fine[1] - theta_fine[0]) / (2 * np.pi))
        frequencies_shifted = np.fft.fftshift(frequencies)  # Shift frequencies
        c_n_shifted = np.fft.fftshift(c_n)  # Shift Fourier coefficients

        axs_fpFT[0].plot(frequencies_shifted, c_n_shifted.real, color=color, label=f'Time: {time}')
        axs_fpFT[0].set(title="Real part of Fourier Coefficients", xlabel="Frequency (n)", ylabel="Re(c_n)", xlim=(0, 20), ylim=(0, 4))
        axs_fpFT[1].plot(frequencies_shifted, c_n_shifted.imag, color=color, label=f'Time: {time}')
        axs_fpFT[1].set(title="Imaginary part of Fourier Coefficients", xlabel="Frequency (n)", ylabel="Im(c_n)", xlim=(0, 20))
        axs_fpFT[2].plot(frequencies_shifted, np.abs(c_n_shifted), color=color, label=f'Time: {time}')
        axs_fpFT[2].set(title="Magnitude of Fourier Coefficients", xlabel="Frequency (n)", ylabel="|c_n|", xlim=(0, 20), ylim=(0, 4))
        axs_fpFT[3].plot(frequencies_shifted, np.angle(c_n_shifted), color=color, label=f'Time: {time}')
        axs_fpFT[3].set(title="Phase of Fourier Coefficients", xlabel="Frequency (n)", ylabel="Phase(c_n)", xlim=(0, 20))
    
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='winter'), ax=ax_fp)
    cbar.set_label('time')
    ax_fp.set(xlabel=r'$\theta$', ylabel=r'$r(\theta, t)$', xlim=(0, 2*np.pi))
    angles = np.linspace(0, 2 * np.pi, 13)  # 0, pi/6, pi/3, pi/2, ..., 2pi
    labels = ["0", r"1/6 $\pi$", r"1/3 $\pi$", r"1/2 $\pi$", r"2/3 $\pi$", r"5/6 $\pi$", r"$\pi$", 
              r"7/6 $\pi$", r"4/3 $\pi$", r"3/2 $\pi$", r"5/3 $\pi$", r"11/6 $\pi$", r"2 $\pi$"]
    ax_fp.set_xticks(angles)
    ax_fp.set_xticklabels(labels)
    fig_fp.tight_layout()
    fig_fp.savefig(fp_fig_name)

    cbar2 = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='winter'), ax=axs_fpFT[0])
    cbar2.set_label('time')
    fig_fpFT.tight_layout()
    fig_fpFT.savefig(fp_FT_fig_name)
    return


######################################################
if __name__ == "__main__": 
    dir_list = ["CALCS_diam15.4_thick3.5/-4.0_-4.0_repeat_1/", 
                "CALCS_diam15.4_thick3.5/-3.5_-3.5_repeat_1/", 
                "CALCS_diam15.4_thick3.5/-3.0_-3.0_repeat_1/", 
                "CALCS_diam15.4_thick3.5/-2.5_-2.5_repeat_1/", 
                "CALCS_diam15.4_thick3.5/-2.0_-2.0_repeat_1/"]
    # dir_list = ["CALCS_test_small/test_repeat_1/"]

    for dir in dir_list:
        print(dir)
        calc_fingerprint_function(f"{dir}")
        plot_fingerprint_function_and_FT(f"{dir}fingerprint.dat.gz", f"{dir}fingerprint.pdf", f"{dir}fingerprint_FT.pdf")

        shutil.copy(f'{dir}/fingerprint.pdf', f"PLOTS_roughness/{dir.split('/')[-2].replace('-', 'n').replace('.', 'p')}_fingerprint.pdf")
        shutil.copy(f'{dir}/fingerprint_FT.pdf', f"PLOTS_roughness/{dir.split('/')[-2].replace('-', 'n').replace('.', 'p')}_fingerprint_FT.pdf")