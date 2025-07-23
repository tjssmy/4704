"""
Simplified 2D Cylindrical Wave Equation Solver
==============================================

A more direct approach using proper finite difference discretization
for the cylindrical Laplacian operator.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv, jn_zeros
from scipy.linalg import eigh

def solve_cylindrical_modes():
    """
    Solve the 2D wave equation in cylindrical coordinates using 
    a direct finite difference approach.
    """
    # Parameters
    R = 1.0  # Domain radius
    nr = 100  # Number of radial points
    
    # Create radial grid (excluding r=0 and r=R for interior problem)
    r = np.linspace(R/(nr+1), R*nr/(nr+1), nr)
    dr = r[1] - r[0]
    
    def solve_mode_m(m):
        """Solve for azimuthal mode number m"""
        # Build the matrix for (1/r)(d/dr)(r du/dr) - (m²/r²)u = -k²u
        A = np.zeros((nr, nr))
        
        for i in range(nr):
            ri = r[i]
            
            if i == 0:
                # Forward difference for first point
                A[i, i] = -2/(dr**2) - m**2/ri**2 + 1/(ri*dr)
                A[i, i+1] = 2/(dr**2) - 1/(ri*dr)
            elif i == nr-1:
                # Backward difference for last point  
                A[i, i-1] = 2/(dr**2) + 1/(ri*dr)
                A[i, i] = -2/(dr**2) - m**2/ri**2 - 1/(ri*dr)
            else:
                # Central difference for interior points
                A[i, i-1] = 1/(dr**2) - 1/(2*ri*dr)
                A[i, i] = -2/(dr**2) - m**2/ri**2
                A[i, i+1] = 1/(dr**2) + 1/(2*ri*dr)
        
        # Solve eigenvalue problem
        eigenvals, eigenvecs = eigh(-A)
        
        # Sort and take positive eigenvalues
        idx = eigenvals > 0
        eigenvals = eigenvals[idx]
        eigenvecs = eigenvecs[:, idx]
        
        # Sort by eigenvalue
        sort_idx = np.argsort(eigenvals)
        eigenvals = eigenvals[sort_idx]
        eigenvecs = eigenvecs[:, sort_idx]
        
        # Convert to k values
        k_values = np.sqrt(eigenvals)
        
        return k_values, eigenvecs, r
    
    # Solve for different m values
    results = {}
    for m in [0, 1, 2]:
        k_vals, eigenvecs, r_grid = solve_mode_m(m)
        results[m] = {'k': k_vals[:4], 'eigenvecs': eigenvecs[:, :4], 'r': r_grid}
    
    return results

def plot_comparison_and_modes():
    """
    Create plots comparing numerical and analytical solutions
    """
    R = 1.0
    
    # Get analytical eigenvalues
    analytical = {
        0: jn_zeros(0, 4) / R,  # J0 zeros
        1: jn_zeros(1, 4) / R,  # J1 zeros  
        2: jn_zeros(2, 4) / R   # J2 zeros
    }
    
    # Solve numerically
    results = solve_cylindrical_modes()
    
    # Print comparison
    print("Eigenvalue Comparison (k_mn values):")
    print("=" * 50)
    print("Mode (m,n)  Numerical   Analytical   Error (%)")
    print("-" * 50)
    
    total_error = 0
    count = 0
    
    for m in [0, 1, 2]:
        for n in range(3):  # First 3 modes for each m
            k_num = results[m]['k'][n]
            k_ana = analytical[m][n]
            error = abs(k_num - k_ana) / k_ana * 100
            print(f"({m},{n+1})      {k_num:.6f}   {k_ana:.6f}    {error:.3f}")
            total_error += error
            count += 1
    
    print("-" * 50)
    print(f"Average error: {total_error/count:.3f}%")
    print()
    
    # Create visualization
    fig = plt.figure(figsize=(16, 10))
    
    # Create coordinate grids for plotting
    nr_plot = 80
    nphi_plot = 100
    r_plot = np.linspace(0, R, nr_plot)
    phi_plot = np.linspace(0, 2*np.pi, nphi_plot)
    R_mesh, PHI_mesh = np.meshgrid(r_plot, phi_plot)
    X_mesh = R_mesh * np.cos(PHI_mesh)
    Y_mesh = R_mesh * np.sin(PHI_mesh)
    
    # Plot first 8 modes
    modes_to_plot = [
        (0, 1), (0, 2), (0, 3), (1, 1),
        (1, 2), (1, 3), (2, 1), (2, 2)
    ]
    
    for idx, (m, n) in enumerate(modes_to_plot):
        ax = plt.subplot(2, 4, idx + 1)
        
        # Get analytical mode
        k_mn = analytical[m][n-1]
        
        # Create 2D mode shape
        U = np.zeros_like(R_mesh)
        for i in range(nr_plot):
            for j in range(nphi_plot):
                r_val = r_plot[i]
                phi_val = phi_plot[j]
                
                if r_val == 0:
                    if m == 0:
                        radial = 1.0
                    else:
                        radial = 0.0
                else:
                    radial = jv(m, k_mn * r_val)
                
                if m == 0:
                    angular = 1.0
                else:
                    angular = np.cos(m * phi_val)
                
                U[j, i] = radial * angular
        
        # Plot
        levels = np.linspace(-1, 1, 21)
        im = ax.contourf(X_mesh, Y_mesh, U, levels=levels, cmap='RdBu_r', extend='both')
        ax.contour(X_mesh, Y_mesh, U, levels=10, colors='black', alpha=0.3, linewidths=0.5)
        
        # Add boundary
        theta = np.linspace(0, 2*np.pi, 100)
        ax.plot(np.cos(theta), np.sin(theta), 'k-', linewidth=2)
        
        ax.set_aspect('equal')
        ax.set_xlim(-1.1, 1.1)
        ax.set_ylim(-1.1, 1.1)
        ax.set_xticks([])
        ax.set_yticks([])
        
        # Title
        k_num = results[m]['k'][n-1]
        k_ana = analytical[m][n-1]
        error = abs(k_num - k_ana) / k_ana * 100
        
        title_str = f'Mode ({m},{n})\nk = {k_ana:.3f}\nFD Error: {error:.2f}%'
        ax.set_title(title_str, fontsize=10)
        
        # Add colorbar
        plt.colorbar(im, ax=ax, shrink=0.8)
    
    plt.suptitle('Cylindrical Wave Modes: Analytical Solutions\n∇²u = -k²u, u(R,φ) = 0', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('eigenvalues_cylindrical.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    plot_comparison_and_modes()
