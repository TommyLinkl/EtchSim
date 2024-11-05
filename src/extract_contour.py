import numpy as np
import shutil
import gzip
from scipy.interpolate import splprep, splev
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.collections import LineCollection
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


def order_points_by_angle(x, y):
    '''
    Function to order points by polar angle around the centroid
    '''

    # Centroid
    centroid_x = np.mean(x)
    centroid_y = np.mean(y)

    angles = np.arctan2(y - centroid_y, x - centroid_x)
    sorted_indices = np.argsort(angles)
    
    return x[sorted_indices], y[sorted_indices]


def smooth_contour_spline(x, y, smoothness=1.0):
    # Ensure the boundary is closed by repeating the first point at the end
    x = np.append(x, x[0])
    y = np.append(y, y[0])

    # Fit a parametric spline to the boundary points
    tck, u = splprep([x, y], s=smoothness, per=1)  # 'per=1' ensures it's periodic (closed contour)

    # Evaluate the spline fit to get smooth contour points
    smooth_points = splev(np.linspace(0, 1, 500), tck)
    smooth_x, smooth_y = smooth_points[0], smooth_points[1]
    smoothed_contour = np.column_stack((smooth_x, smooth_y))
    return smoothed_contour


def wrap_around_list(index, maximum, n_points):
    """
    Generate a list of neighboring indices with circular indexing.
    """
    half_window = n_points // 2
    naive_list = np.arange(index - half_window, index + half_window + 1)

    wrap_list = np.mod(naive_list, maximum)
    
    return wrap_list


def smooth_contour_interpolate(x, y, n_points=4.0):
    """
    Smooth the edge of a curve using a moving average of neighboring points.
    """
    x = np.append(x, x[0])
    y = np.append(y, y[0])
    smooth_x = np.zeros_like(x)
    smooth_y = np.zeros_like(y)
    
    # Loop through each point and calculate the average of its neighbors
    for count in range(len(x)):
        list_temp = wrap_around_list(count, len(x), n_points)
        x_temp = x[list_temp]
        y_temp = y[list_temp]
        smooth_x[count] = np.mean(x_temp)
        smooth_y[count] = np.mean(y_temp)
    smoothed_contour = np.column_stack((smooth_x, smooth_y))
    return smoothed_contour


def reconstruct_boundary_pts_over_time(read_filename): 
    with gzip.open(read_filename, 'rt') as f: 
        data = np.loadtxt(f)

    times = np.unique(data[:, 0])
    blocks = []
    for time in times:
        block = data[data[:, 0] == time]
        blocks.append(block)

    return blocks, times


def calculate_curvature(contour, pad=100):
    # Pad the contour with the first and last 'pad' points to handle boundary conditions
    padded_contour = np.vstack((contour[-pad:], contour, contour[:pad]))

    # First derivatives with respect to the padded contour
    dx = np.gradient(padded_contour[:, 0], edge_order=2)
    dy = np.gradient(padded_contour[:, 1], edge_order=2)

    # Second derivatives
    ddx = np.gradient(dx, edge_order=2)
    ddy = np.gradient(dy, edge_order=2)

    # Calculate curvature using the curvature formula
    curvature = np.abs(dx * ddy - dy * ddx) / ((dx**2 + dy**2)**1.5 + 1e-8)

    # print(curvature[pad:-pad])
    curvature[-pad-1] = (curvature[-pad-2] + curvature[0])/2.0
    # Return only the curvature corresponding to the original contour, excluding the padding
    return curvature[pad:-pad]


def approximate_curvature(contour):
    # Calculate vectors between consecutive points
    vectors = np.diff(contour, axis=0)

    # Include the wrap-around vectors
    wrap_around_vector = contour[0] - contour[-1]
    vectors = np.vstack([vectors, wrap_around_vector])
    
    # Normalize vectors
    norms = np.linalg.norm(vectors, axis=1, keepdims=True)
    unit_vectors = vectors / norms
    
    # Calculate angles between consecutive vectors
    dot_products = np.einsum('ij,ij->i', unit_vectors[:-1], unit_vectors[1:])
    angles = np.arccos(np.clip(dot_products, -1.0, 1.0))  # Clip for numerical stability

    # Approximate curvature as the angle divided by the distance between points
    # Wrap the distance calculation
    distances = np.linalg.norm(contour[np.arange(len(contour)), None] - contour[np.arange(len(contour)), None], axis=1)
    
    curvatures = angles / np.linalg.norm(contour[2:] - contour[:-2], axis=1)
    
    # Pad the curvature array to match the contour length
    curvatures = np.pad(curvatures, (1, 1), mode='edge')

    return curvatures


def get_contour(calc_dir, mode='interpolate', smoothness=1.0, every_frame=1, gen_plot=True, plotOriginalPts=False):
    # Read in file "vacXY_boundary_pts.dat.gz"
    boundary_pts_list, unique_times = reconstruct_boundary_pts_over_time(f"{calc_dir}vacXY_boundary_pts.dat.gz")

    contourPts_zip = f"{calc_dir}vacXY_contour_pts.dat.gz"
    with gzip.open(contourPts_zip, 'wt') as f:
        f.write("# time      x       y\n")
    if gen_plot: 
        fig, ax = plt.subplots(1,1, figsize=(7,5))
    norm = plt.Normalize(vmin=np.unique(boundary_pts_list[0][:, 0]), vmax=np.unique(boundary_pts_list[-1][:, 0]))
    colors = plt.cm.spring(norm(unique_times))
    
    print(len(boundary_pts_list[every_frame::every_frame]))

    for frame, color in zip(boundary_pts_list[every_frame::every_frame], colors[every_frame::every_frame]):
        simTime = np.unique(frame[:, 0])[0]
        boundary_coordXY = frame[:, 1:]

        if len(boundary_coordXY) <= 5: 
            continue

        # Order the points correctly
        x_ordered, y_ordered = order_points_by_angle(boundary_coordXY[:,0], boundary_coordXY[:,1])

        # Smooth the contour
        if mode=='spline':
            smoothed_contour_pts = smooth_contour_spline(x_ordered, y_ordered, smoothness)
        elif mode=='interpolate': 
            smoothed_contour_pts = smooth_contour_interpolate(x_ordered, y_ordered, smoothness)
            smoothed_contour_pts = smooth_contour_spline(smoothed_contour_pts[:, 0], smoothed_contour_pts[:, 1], 0.0)

        # Record contour points
        with gzip.open(contourPts_zip, 'at') as file:
            for point in smoothed_contour_pts: 
                file.write(f"{simTime}     {point[0]}      {point[1]}\n")

        # Plot the contour
        if gen_plot: 
            if plotOriginalPts: 
                ax.scatter(x_ordered, y_ordered, color=color, s=8, label='Original Boundary Points')
            '''
            smoothed_contour_pts = np.vstack([smoothed_contour_pts, smoothed_contour_pts[0]])
            ax.plot(smoothed_contour_pts[:,0], smoothed_contour_pts[:,1], color=color, label='Smooth Contour', lw=2)
            '''

            # smoothed_contour_pts = np.vstack([smoothed_contour_pts, smoothed_contour_pts[0]])
            curvature = calculate_curvature(smoothed_contour_pts)
            curvature[curvature <= 0.12] = curvature.min() 
            curvature[curvature >= 0.12] = curvature.max() 
            sc = ax.scatter(smoothed_contour_pts[:, 0], smoothed_contour_pts[:, 1], c=curvature, cmap='bwr', s=8, lw=0)
            
            '''
            points = smoothed_contour_pts.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)

            lc = LineCollection(segments, cmap='viridis', norm=plt.Normalize(0, curvature.max()))
            lc.set_array(curvature)
            lc.set_linewidth(2)
            ax.add_collection(lc)
            ax.set_xlim(smoothed_contour_pts[:, 0].min(), smoothed_contour_pts[:, 0].max())
            ax.set_ylim(smoothed_contour_pts[:, 1].min(), smoothed_contour_pts[:, 1].max())
            ax.set_aspect('equal', 'box')
            '''


    if gen_plot: 
        # cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='spring'), ax=ax)
        # cbar.set_label('time')

        # cbar = plt.colorbar(lc, ax=ax)
        # cbar.set_label('Curvature')

        cbar = plt.colorbar(sc)
        cbar.set_label('Curvature')

        fig.tight_layout()
        fig.savefig(f"{calc_dir}vacXY_contour.pdf")
    return


######################################################
if __name__ == "__main__": 
    '''
    # dir_list = ["CALCS_diam15.4_thick3.5/-4.0_-4.0_repeat_1/", "CALCS_diam15.4_thick3.5/-3.5_-3.5_repeat_1/", "CALCS_diam15.4_thick3.5/-3.0_-3.0_repeat_1/", "CALCS_diam15.4_thick3.5/-2.5_-2.5_repeat_1/", "CALCS_diam15.4_thick3.5/-2.0_-2.0_repeat_1/"]
    # every_frame = [20, 20, 20, 20, 185]

    dir_list = ["CALCS_uneven/-3.0_-2.0_repeat_1/", "CALCS_uneven/-3.0_-2.5_repeat_1/", "CALCS_uneven/-4.0_-2.0_repeat_1/", "CALCS_uneven/-4.0_-2.5_repeat_1/", "CALCS_uneven/-4.0_-3.0_repeat_1/", "CALCS_uneven/-4.0_-3.5_repeat_1/"]
    every_frame = [20, 20, 20, 20, 20, 20]

    for i, dir in enumerate(dir_list):
        print(dir)
        get_contour(f"{dir}", mode='interpolate', smoothness=2, every_frame=every_frame[i], gen_plot=True, plotOriginalPts=False)
        # get_contour(f"{dir}", mode='spline', smoothness=5.0, every_frame=every_frame[i], plotOriginalPts=False)
        
        shutil.copy(f'{dir}/vacXY_contour.pdf', f"PLOTS_roughness/{dir.split('/')[-2].replace('-', 'n').replace('.', 'p')}_vacXY_contour.pdf")
    '''


    dir_list = []
    every_frame = []

    chemPots = ["-2.5_-2.0", "-3.0_-2.0", "-3.0_-2.5", "-3.5_-2.0", "-3.5_-2.5", "-3.5_-3.0", "-4.0_-2.0", "-4.0_-2.5", "-4.0_-3.0", "-4.0_-3.5"]
    repeats = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
    for r, frame_val in zip(chemPots, repeats):
        dir_list.extend([f"CALCS_uneven/{r}_repeat_{k}/" for k in range(1, 17)])
        every_frame.extend([frame_val] * 16)

    for i, dir in enumerate(dir_list):
        print(dir)
        get_contour(f"{dir}", mode='interpolate', smoothness=2, every_frame=every_frame[i], gen_plot=False)