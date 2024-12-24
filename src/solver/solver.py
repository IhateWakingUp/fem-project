import numpy as np
from scipy.optimize import fsolve
from CST import CST
from Assembly import Assembly

class Solver:
    def __init__(self, E=200e9, nu=0.3, thickness=0.01):
        """Initialize solver with material properties"""
        self.E = E
        self.nu = nu
        self.thickness = thickness
        self.assembly = Assembly(E, nu, thickness)
        self.elements = []
        self.node_mappings = []
        self.coordinates_list = []
        
    def add_elements(self, coordinates_list, node_mappings):
        """Add multiple elements at once"""
        self.coordinates_list = coordinates_list
        self.node_mappings = node_mappings
        
        # Create elements and add to assembly
        for coords in coordinates_list:
            self.elements.append(CST(coords))
            self.assembly.add_element(coords)
            
    def D(self):
        """Calculate D matrix for plane stress"""
        return self.E * np.array([[1, self.nu, 0],
                                [self.nu, 1, 0],
                                [0, 0, (1-self.nu)/2]]) / (1-(self.nu)**2)
    
    def disp(self, K, forces, fixed_dofs):
        """Calculate displacements"""
        m = K.shape[0]
        f = np.array(forces).reshape(-1, 1)
        free_dofs = list(set(range(m)) - set(fixed_dofs))
        K_reduced = K[np.ix_(free_dofs, free_dofs)]
        f_reduced = f[free_dofs]
        
        d_reduced = np.linalg.solve(K_reduced, f_reduced)
        d = np.zeros((m, 1))
        d[free_dofs] = d_reduced
        return d
    
    def strain(self, element, node_mapping, d):
        """Calculate strain for a single element"""
        if not hasattr(element, 'B'):
            element.B_matrix()
        elem_dofs = []
        for node in node_mapping:
            elem_dofs.extend([2*int(node), 2*int(node)+1])
        elem_disp = d[elem_dofs]
        return element.B @ elem_disp
    
    def stress(self, eps):
        """Calculate stress"""
        D_matrix = self.D()
        return D_matrix @ eps
    
    def Von_Mises(self, sigma):
        """Calculate von Mises stress"""
        return np.sqrt(sigma[0]**2 + sigma[1]**2 - sigma[0]*sigma[1] + 3*sigma[2]**2)
    
    def analyze(self, forces, fixed_dofs, n_nodes, sigma_yield=250e6):
        """Perform complete structural analysis"""
        # Get global stiffness matrix
        K = self.assembly.assemble_global_stiffness(n_nodes, self.node_mappings)
        
        # Calculate displacements
        d = self.disp(K, forces, fixed_dofs)
        
        # Calculate stresses and strains for each element
        vm_stresses = []
        strains = []
        
        for element, nodes in zip(self.elements, self.node_mappings):
            eps = self.strain(element, nodes, d)
            strains.append(eps)
            sigma = self.stress(eps)
            vm_stress = float(self.Von_Mises(sigma))
            vm_stresses.append(vm_stress)
        
        # Safety analysis
        max_vm_stress = max(vm_stresses)
        safety_factor = float(sigma_yield / max_vm_stress)
        
        if max_vm_stress > sigma_yield:
            print(f"Structure will fail. Maximum stress {max_vm_stress/1e6:.2f} MPa exceeds yield strength")
        else:
            print(f"Structure is safe. Factor of safety: {safety_factor:.2f}")
        
        return d, strains, vm_stresses, safety_factor 
