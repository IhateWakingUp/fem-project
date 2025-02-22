from Solver import Solver
import numpy as np
# Example usage
solve = Solver(E=200e9, nu=0.3, thickness=0.01)

# Define coordinates and node mappings
coordinates_list = [
    [(0,0), (1,0), (0.5,0.5)],    # Element 1
    [(1,1), (1,0), (0.5,0.5)],    # Element 2
    [(1,1), (0,1), (0.5,0.5)],    # Element 3
    [(0,0), (0,1), (0.5,0.5)]     # Element 4
]

node_mappings = [
    [0, 1, 4],  # Element 1 nodes
    [1, 2, 4],  # Element 2 nodes
    [2, 3, 4],  # Element 3 nodes
    [0, 3, 4]   # Element 4 nodes
]

# Add elements
solve.add_elements(coordinates_list, node_mappings)

# Define boundary conditions and loads
fixed_dofs = [0, 1, 6]
forces = np.zeros(10)
forces[1] = 0.5
forces[2] = 0.5

# Perform analysis
d, strains, vm_stresses, safety_factor = solve.analyze(forces, fixed_dofs, n_nodes=5)

# Calculate strain and stress for specific element
triangle_1 = CST(coordinates_list[0])
eps_1 = solve.strain(triangle_1, node_mappings[0], d)  
sigma_1 = solve.stress(eps_1)  
