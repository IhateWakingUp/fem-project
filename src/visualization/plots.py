import matplotlib.pyplot as plt
def plot_shape_functions(triangle, num_points=50):
    # Create grid of points
    x = np.linspace(min(triangle.coordinates[:,0]), max(triangle.coordinates[:,0]), num_points)
    y = np.linspace(min(triangle.coordinates[:,1]), max(triangle.coordinates[:,1]), num_points)
    X, Y = np.meshgrid(x, y)
    
    # Calculate shape functions at each point
    Z1 = np.zeros_like(X)
    Z2 = np.zeros_like(X)
    Z3 = np.zeros_like(X)
    
    for i in range(num_points):
        for j in range(num_points):
            N1, N2, N3 = triangle.shape_func(X[i,j], Y[i,j])
            Z1[i,j] = N1
            Z2[i,j] = N2
            Z3[i,j] = N3
    
    # Create 3D plots
    fig = plt.figure(figsize=(15, 5))
    
    # Get triangle vertices with the first point repeated to close the triangle
    x_coords = np.append(triangle.coordinates[:,0], triangle.coordinates[0,0])
    y_coords = np.append(triangle.coordinates[:,1], triangle.coordinates[0,1])
    z_coords = np.zeros_like(x_coords)
    
    # Generate equations using coefficients
    equations = []
    for i in range(3):
        a0 = triangle.Coeff[i,0]/(2*triangle.area)
        a1 = triangle.Coeff[i,1]/(2*triangle.area)
        a2 = triangle.Coeff[i,2]/(2*triangle.area)
        
        # Format coefficients for better readability
        terms = []
        if abs(a0) > 1e-10:  # Add constant term if not approximately zero
            terms.append(f"{a0:.3f}")
        if abs(a1) > 1e-10:  # Add x term if not approximately zero
                terms.append(f"{a1:.3f}x")
        if abs(a2) > 1e-10:  # Add y term if not approximately zero
                terms.append(f"{a2:.3f}y")
                
        equation = r'$N_{' + str(i+1) + r'}(x,y) = ' + ' + '.join(terms).replace('+ -', '- ') + '$'
        equations.append(equation)
    
    titles = ['Shape Function N1', 'Shape Function N2', 'Shape Function N3']
    surfaces = [Z1, Z2, Z3]
    
    for idx, (title, Z, eq) in enumerate(zip(titles, surfaces, equations), 1):
        ax = fig.add_subplot(1, 3, idx, projection='3d')
        
        # Plot shape function surface
        surf = ax.plot_surface(X, Y, Z, cmap='viridis', alpha=0.8)
        
        # Plot triangle outline (always at z=0)
        ax.plot(x_coords, y_coords, z_coords, 'k-', linewidth=2, label='Triangle')
        # Plot triangle vertices (always at z=0)
        ax.scatter(x_coords[:-1], y_coords[:-1], z_coords[:-1], 
                  color='red', s=100, label='Vertices')
        
        # Add title with equation
        ax.set_title(f'{title}\n{eq}', pad=20)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('N')
        ax.set_zlim(0, 1.0)
        fig.colorbar(surf, ax=ax)
        ax.legend()
    
    plt.tight_layout()
    plt.show()
coordinates = [(0,0), (1,0), (0.5,0.5)]
triangle = CST(coordinates)
plot_shape_functions(triangle)

 
