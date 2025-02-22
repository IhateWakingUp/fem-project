import numpy as np
from CST import CST

class Assembly:
    def __init__(self, E, nu, thickness):
        self.E = E
        self.nu = nu
        self.thickness = thickness
        self.elements = []
        self.coordinates_list = []
        
    def add_element(self, coordinates):
        """Add element coordinates and create CST element"""
        self.coordinates_list.append(coordinates)
        element = CST(coordinates)
        element.Stiff(self.E, self.nu, self.thickness)
        self.elements.append(element)
            
    def assemble_global_stiffness(self, n_nodes, node_mappings):
        #n_nodes: Total number of nodes in the system
        #node_mappings: List of node numbers for each element
        
        K_global = np.zeros((2*n_nodes, 2*n_nodes))
        
        for elem, nodes in zip(self.elements, node_mappings):
            K_local = elem.K
            dofs = []
            for node in nodes:
                dofs.extend([2*node, 2*node+1])
            
            for i in range(6):
                for j in range(6):
                    K_global[dofs[i], dofs[j]] += K_local[i,j]
                    
        # Using symmetry to fill rest of the matrix
        K_global = K_global + K_global.T - np.diag(np.diag(K_global))
        return K_global
 
