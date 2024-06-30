import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Given points
points = np.array([
    [2.085, 1.20377531, -2.9925],
    [0.0, 2.40755062, 4.275e-01],
    [2.085, 1.20377531, -4.275e-01],
    [0.0, 2.40755062, 2.9925]
])

# Given unit cell vectors
a = 4.17  # AA
c = 6.84  # AA
unitCellVectors = np.array([
    [a, 0, 0],
    [-0.5*a, np.sqrt(3)/2*a, 0],
    [0, 0, c]
])

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the points
ax.scatter(points[0:2, 0], points[0:2, 1], points[0:2, 2], c='r', marker='o', label='First two points')
ax.scatter(points[2:4, 0], points[2:4, 1], points[2:4, 2], c='b', marker='o', label='Last two points')

origin = np.array([0, 0, 0])  # origin point
for vec in unitCellVectors:
    ax.quiver(origin[0], origin[1], origin[2], vec[0], vec[1], vec[2], 
              color='b', length=1.0, normalize=False)

# Labels and plot settings
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Set the same scale for all axes
ax.set_box_aspect([1,1,1])
max_range = 16
mid_x = (points[:, 0].max()+points[:, 0].min()) * 0.5
mid_y = (points[:, 1].max()+points[:, 1].min()) * 0.5
mid_z = (points[:, 2].max()+points[:, 2].min()) * 0.5
ax.set_xlim(mid_x - max_range/2, mid_x + max_range/2)
ax.set_ylim(mid_y - max_range/2, mid_y + max_range/2)
ax.set_zlim(mid_z - max_range/2, mid_z + max_range/2)

plt.show()
