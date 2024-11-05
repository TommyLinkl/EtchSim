import pickle, gzip, json
import numpy as np
from math import floor
from .constants import *

def xyz_to_lammps_dump(xyz_file, dump_file):
    with open(xyz_file, 'r') as infile:
        # Read the first two lines: atom count and comment
        num_atoms = int(infile.readline().strip())
        infile.readline()  # Skip comment line

        # Prepare the output dump file
        with open(dump_file, 'w') as outfile:
            # Write LAMMPS dump format headers
            outfile.write("ITEM: TIMESTEP\n0\n")
            outfile.write(f"ITEM: NUMBER OF ATOMS\n{num_atoms}\n")
            outfile.write("ITEM: BOX BOUNDS pp pp pp\n-100000 100000\n-100000 100000\n-100000 100000\n")
            outfile.write("ITEM: ATOMS id type x y z user\n")

            # Iterate through the atoms and write data to the dump file
            for i, line in enumerate(infile, start=1):
                parts = line.split()
                atom_type = 1 if parts[0] == 'P' else 2  # Assign a type based on element
                x, y, z = parts[1:4]
                user = parts[4] if len(parts) > 4 else 0  # Use 0 if no fifth column
                outfile.write(f"{i} {atom_type} {x} {y} {z} {user}\n")

    print(f"Conversion complete! LAMMPS dump saved as '{dump_file}'.")


def xyz_to_pdb(xyz_file, pdb_file):
    with open(xyz_file, 'r') as infile, open(pdb_file, 'w') as outfile:
        # Read the number of atoms and skip the comment line
        num_atoms = int(infile.readline().strip())
        infile.readline()  # Skip comment line

        # Process each line of atom data
        for atom_id, line in enumerate(infile, start=1):
            parts = line.split()
            element = parts[0]  # P or In
            x, y, z = map(float, parts[1:4])
            beta = float(parts[4]) if len(parts) > 4 else 0.00  # Use 0.00 if missing

            # Write the atom data in PDB format
            outfile.write(
                f"ATOM  {atom_id:5d}  {element:<2}  MOL     1    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 {beta:6.2f}\n"
            )

    print(f"Conversion complete! PDB file saved as '{pdb_file}'.")


def load_xyz(filename):
    """Load an .xyz file and return atom types and coordinates."""
    with open(filename, 'r') as f:
        num_atoms = int(f.readline().strip())
        f.readline()  # Skip comment line
        atoms = []
        for _ in range(num_atoms):
            parts = f.readline().split()
            atom_type = parts[0]
            coords = list(map(float, parts[1:4]))
            atoms.append((atom_type, np.array(coords)))
    return atoms


def save_xyz_with_neighbors(atoms, neighbors, filename):
    """Save atoms with their neighbor counts as the fifth column."""
    with open(filename, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write("Atoms with neighbor counts as the fifth column\n")
        for (atom_type, coords), neighbor_count in zip(atoms, neighbors):
            f.write(f"{atom_type} {coords[0]:.6f} {coords[1]:.6f} {coords[2]:.6f} {neighbor_count}\n")


def build_grid(atoms, cell_size):
    """Partition the atoms into 3D grids for efficient neighbor search."""
    grid = {}
    for idx, (_, pos) in enumerate(atoms):
        # Compute the cell index (ix, iy, iz) for each atom
        cell = tuple(floor(pos[i] / cell_size[i]) for i in range(3))
        if cell not in grid:
            grid[cell] = []
        grid[cell].append(idx)
    return grid


def neighbors_in_grid(atoms, grid, b_min, b_max):
    """Calculate the number of neighbors within the given distance range."""
    neighbors = [0] * len(atoms)  # Initialize neighbor counts
    b_min_sq = b_min ** 2
    b_max_sq = b_max ** 2

    # Iterate over all grid cells
    for cell, atom_indices in grid.items():
        # Check the current cell and its neighbors
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                for dz in [-1, 0, 1]:
                    neighbor_cell = (cell[0] + dx, cell[1] + dy, cell[2] + dz)
                    if neighbor_cell not in grid:
                        continue

                    # Compare all atoms in the current cell with its neighbor cell
                    for i in atom_indices:
                        for j in grid[neighbor_cell]:
                            if i >= j:  # Avoid double counting
                                continue
                            dist_sq = np.sum((atoms[i][1] - atoms[j][1]) ** 2)
                            if b_min_sq <= dist_sq <= b_max_sq:
                                neighbors[i] += 1
                                neighbors[j] += 1
    return neighbors


def main(xyz_file, output_xyz_file, output_pdb_file, b_min, b_max):
    # Load atoms from the XYZ file
    atoms = load_xyz(xyz_file)

    # Define cell size for the grid (usually b_max works well)
    cell_size = [b_max] * 3

    # Build the 3D grid for spatial partitioning
    grid = build_grid(atoms, cell_size)

    # Calculate neighbors for each atom within distance bounds
    neighbors = neighbors_in_grid(atoms, grid, b_min, b_max)

    # Save the modified XYZ file with neighbor counts as the fifth column
    save_xyz_with_neighbors(atoms, neighbors, output_xyz_file)

    print(f"Neighbor calculation complete! Results saved to '{output_xyz_file}'.")

    xyz_to_pdb(output_xyz_file, output_pdb_file)


if __name__ == "__main__":
    xyz_file = "CALCS_rotated_hex/init_atoms.xyz"
    output_xyz_file = "CALCS_rotated_hex/init_atoms_neighbors.xyz"
    output_pdb_file = "CALCS_rotated_hex/init_atoms_neighbors.pdb"
    b_min = bond_length_max/2.0
    b_max = bond_length_max

    main(xyz_file, output_xyz_file, output_pdb_file, b_min, b_max)