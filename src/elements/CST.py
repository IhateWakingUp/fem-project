import numpy as np

class CST:
    def __init__(self, coordinates):
        self.coordinates = np.array(coordinates)
        if self.coordinates.shape != (3, 2):
            raise ValueError("CST requires exactly 3 coordinates in 2D")
            
        self.x1 = self.coordinates[0][0]
        self.y1 = self.coordinates[0][1]
        self.x2 = self.coordinates[1][0] 
        self.y2 = self.coordinates[1][1]
        self.x3 = self.coordinates[2][0] 
        self.y3 = self.coordinates[2][1]
        
        self.calc_area()
        self.calc_coeff()
    
    def calc_area(self):  
        self.X = np.array([
            [1, self.x1, self.y1],
            [1, self.x2, self.y2],
            [1, self.x3, self.y3]
        ])
        self.area = abs(np.linalg.det(self.X)/2)  
    
    def calc_coeff(self):
        # Calculate b and c coefficients directly
        self.a = np.array([
            self.x2*self.y3 - self.x3*self.y2,
            self.x3*self.y1 - self.x1*self.y3,
            self.x1*self.y2 - self.y1*self.x2,
        ])
        self.b = np.array([
            self.y2 - self.y3,
            self.y3 - self.y1,
            self.y1 - self.y2
        ])
        
        self.c = np.array([
            self.x3 - self.x2,
            self.x1 - self.x3,
            self.x2 - self.x1
        ])
        self.Coeff = np.vstack((self.a, self.b, self.c)).T
    def B_matrix(self):
        # Form B matrix using b and c coefficients
        self.B = (1/(2*self.area)) * np.array([
            [self.b[0], 0, self.b[1], 0, self.b[2], 0],
            [0, self.c[0], 0, self.c[1], 0, self.c[2]],
            [self.c[0], self.b[0], self.c[1], self.b[1], self.c[2], self.b[2]]
        ])
        return self.B
    
    def D_matrix(self, E, nu):
        return E/(1 - nu**2) * np.array([
            [1, nu, 0],
            [nu, 1, 0],
            [0, 0, (1-nu)/2]
        ])
    
    def Stiff(self, E, nu, t):
        if not hasattr(self, 'B'):
            self.B_matrix()
            
        D = self.D_matrix(E, nu)
        self.K = t * self.area * self.B.T @ D @ self.B
        return self.K

    def shape_func(self, x, y):
        N = np.zeros(3)
        for i in range(3):
            N[i] = (1/(2*self.area)) * (
                self.a[i] + 
                self.b[i]*x + 
                self.c[i]*y
            )
        return tuple(N)

 
