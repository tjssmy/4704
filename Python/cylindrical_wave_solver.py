"""
2D Harmonic Wave Equation Solution in Cylindrical Coordinates
============================================================

Solves the 2D harmonic wave equation in cylindrical coordinates (r, φ):
∇²u = (1/r)(∂/∂r)(r ∂u/∂r) + (1/r²)(∂²u/∂φ²) = -k²u

with boundary condition:
- u(r=R, φ) = 0 (Dirichlet BC at outer boundary)

This leads to eigenvalue problems with solutions of the form:
u(r,φ) = J_m(k_mn * r) * [A*cos(mφ) + B*sin(mφ)]

where J_m are Bessel functions of the first kind and k_mn are eigenvalues.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import eigs
from scipy.special import jv, jn_zeros
import warnings
warnings.filterwarnings('ignore')

class CylindricalWaveSolver:
    def __init__(self, R=1.0, nr=50, nphi=64):
        """
        Initialize the cylindrical wave equation solver
        
        Parameters:
        -----------
        R : float
            Radius of the cylindrical domain
        nr : int
            Number of radial grid points
        nphi : int
            Number of angular grid points
        """
        self.R = R
        self.nr = nr
        self.nphi = nphi
        
        # Create grids
        self.dr = R / (nr - 1)
        self.dphi = 2 * np.pi / nphi
        
        self.r = np.linspace(0, R, nr)
        self.phi = np.linspace(0, 2*np.pi, nphi)
        
        # Create meshgrid for plotting
        self.R_mesh, self.PHI_mesh = np.meshgrid(self.r, self.phi)
        self.X_mesh = self.R_mesh * np.cos(self.PHI_mesh)
        self.Y_mesh = self.R_mesh * np.sin(self.PHI_mesh)
        
    def build_operators(self, m=0):
        """
        Build finite difference operators for given azimuthal mode number m
        
        Parameters:
        -----------
        m : int
            Azimuthal mode number
            
        Returns:
        --------
        A : sparse matrix
            Discrete Laplacian operator
        """
        # For azimuthal mode m, the equation becomes:
        # (1/r)(d/dr)(r du/dr) - (m²/r²)u = -k²u
        
        # Build radial operator using finite differences
        # We exclude r=0 and r=R to avoid singularity and apply BC
        nr_inner = self.nr - 2  # Interior points only
        r_inner = self.r[1:-1]  # r from dr to R-dr
        
        # For the operator (1/r)(d/dr)(r du/dr), use finite differences
        # d/dr(r du/dr) ≈ [r_{i+1/2}(u_{i+1} - u_i)/dr - r_{i-1/2}(u_i - u_{i-1})/dr]/dr
        # where r_{i±1/2} = r_i ± dr/2
        
        main_diag = np.zeros(nr_inner)
        lower_diag = np.zeros(nr_inner-1)
        upper_diag = np.zeros(nr_inner-1)
        
        for i in range(nr_inner):
            r_i = r_inner[i]
            
            if i == 0:  # First interior point
                # Special treatment near r=0
                r_plus = r_i + self.dr/2
                r_minus = max(r_i - self.dr/2, self.dr/4)  # Avoid r=0
                
                main_diag[i] = -(r_plus + r_minus) / (r_i * self.dr**2) - m**2/r_i**2
                upper_diag[i] = r_plus / (r_i * self.dr**2)
                
            elif i == nr_inner-1:  # Last interior point
                r_plus = r_i + self.dr/2
                r_minus = r_i - self.dr/2
                
                main_diag[i] = -(r_plus + r_minus) / (r_i * self.dr**2) - m**2/r_i**2
                lower_diag[i-1] = r_minus / (r_i * self.dr**2)
                
            else:  # General interior points
                r_plus = r_i + self.dr/2
                r_minus = r_i - self.dr/2
                
                main_diag[i] = -(r_plus + r_minus) / (r_i * self.dr**2) - m**2/r_i**2
                lower_diag[i-1] = r_minus / (r_i * self.dr**2)
                upper_diag[i] = r_plus / (r_i * self.dr**2)
        
        # Create sparse matrix
        A = diags([lower_diag, main_diag, upper_diag], 
                 [-1, 0, 1], shape=(nr_inner, nr_inner))
        
        return A, r_inner
    
    def solve_eigenvalue_problem(self, m=0, n_modes=10):
        """
        Solve eigenvalue problem for given azimuthal mode number
        
        Parameters:
        -----------
        m : int
            Azimuthal mode number
        n_modes : int
            Number of eigenvalues to compute
            
        Returns:
        --------
        eigenvals : array
            Eigenvalues (k²)
        eigenvecs : array
            Eigenvectors (radial modes)
        r_inner : array
            Radial grid for interior points
        """
        A, r_inner = self.build_operators(m)
        
        # Solve eigenvalue problem: A*u = lambda*u where lambda = -k²
        eigenvals, eigenvecs = eigs(-A, k=n_modes, which='SM', sigma=0)
        
        # Sort by eigenvalue magnitude
        idx = np.argsort(np.real(eigenvals))
        eigenvals = np.real(eigenvals[idx])
        eigenvecs = np.real(eigenvecs[:, idx])
        
        # Normalize eigenvectors
        for i in range(n_modes):
            eigenvecs[:, i] = eigenvecs[:, i] / np.max(np.abs(eigenvecs[:, i]))
        
        return eigenvals, eigenvecs, r_inner
    
    def get_analytical_solution(self, m, n):
        """
        Get analytical eigenvalues using Bessel function zeros
        
        Parameters:
        -----------
        m : int
            Azimuthal mode number
        n : int
            Radial mode number (1-indexed)
            
        Returns:
        --------
        k_analytical : float
            Analytical eigenvalue k_mn
        """
        # Get the nth zero of Bessel function J_m
        zeros = jn_zeros(m, n)
        k_analytical = zeros[-1] / self.R
        return k_analytical
    
    def create_2d_mode(self, eigenvals, eigenvecs, r_inner, m, mode_idx, mode_type='cos'):
        """
        Create 2D mode from radial eigenfunction
        
        Parameters:
        -----------
        eigenvals : array
            Eigenvalues
        eigenvecs : array
            Radial eigenfunctions
        r_inner : array
            Radial grid for interior points
        m : int
            Azimuthal mode number
        mode_idx : int
            Index of the mode to visualize
        mode_type : str
            'cos' or 'sin' for azimuthal variation
            
        Returns:
        --------
        u_2d : array
            2D mode shape
        """
        # Interpolate radial eigenfunction to full grid
        radial_func = np.zeros(self.nr)
        radial_func[1:-1] = eigenvecs[:, mode_idx]
        radial_func[0] = 0  # Boundary condition at r=0 (regularity)
        radial_func[-1] = 0  # Boundary condition at r=R
        
        # Create 2D mode
        u_2d = np.zeros((self.nphi, self.nr))
        
        for i in range(self.nphi):
            phi_val = self.phi[i]
            if mode_type == 'cos':
                angular_factor = np.cos(m * phi_val)
            else:
                angular_factor = np.sin(m * phi_val)
            
            u_2d[i, :] = radial_func * angular_factor
        
        return u_2d

def plot_geometry():
    """
    Plot the cylindrical geometry and discretization
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot 1: Geometry
    theta = np.linspace(0, 2*np.pi, 100)
    R = 1.0
    x_boundary = R * np.cos(theta)
    y_boundary = R * np.sin(theta)
    
    ax1.plot(x_boundary, y_boundary, 'k-', linewidth=2, label='Boundary u=0')
    ax1.fill(x_boundary, y_boundary, alpha=0.1, color='lightblue')
    ax1.set_xlim(-1.2, 1.2)
    ax1.set_ylim(-1.2, 1.2)
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_title('Cylindrical Domain\n∇²u = -k²u, u(R,φ) = 0')
    ax1.legend()
    
    # Add coordinate system
    ax1.arrow(0, 0, 0.8, 0, head_width=0.05, head_length=0.05, fc='red', ec='red')
    ax1.arrow(0, 0, 0, 0.8, head_width=0.05, head_length=0.05, fc='red', ec='red')
    ax1.text(0.85, -0.1, 'x', fontsize=12, color='red')
    ax1.text(-0.1, 0.85, 'y', fontsize=12, color='red')
    
    # Plot 2: Grid discretization
    solver = CylindricalWaveSolver(R=1.0, nr=15, nphi=24)
    
    # Plot radial lines
    for i in range(0, solver.nphi, 3):
        phi_val = solver.phi[i]
        r_line = np.linspace(0, R, 20)
        x_line = r_line * np.cos(phi_val)
        y_line = r_line * np.sin(phi_val)
        ax2.plot(x_line, y_line, 'b-', alpha=0.5, linewidth=0.8)
    
    # Plot circular lines
    for i in range(1, solver.nr, 2):
        r_val = solver.r[i]
        theta_circle = np.linspace(0, 2*np.pi, 100)
        x_circle = r_val * np.cos(theta_circle)
        y_circle = r_val * np.sin(theta_circle)
        ax2.plot(x_circle, y_circle, 'g-', alpha=0.5, linewidth=0.8)
    
    # Boundary
    ax2.plot(x_boundary, y_boundary, 'k-', linewidth=2)
    
    ax2.set_xlim(-1.2, 1.2)
    ax2.set_ylim(-1.2, 1.2)
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_title('Finite Difference Grid\n(r,φ) coordinates')
    
    plt.tight_layout()
    plt.savefig('geometry_cylindrical.png', dpi=300, bbox_inches='tight')
    plt.show()

def solve_and_visualize_modes():
    """
    Solve for eigenvalues and visualize the first 8 modes
    """
    solver = CylindricalWaveSolver(R=1.0, nr=40, nphi=64)
    
    # Collect modes from different m values
    modes_data = []
    
    # m = 0 modes (axisymmetric)
    eigenvals_0, eigenvecs_0, r_inner_0 = solver.solve_eigenvalue_problem(m=0, n_modes=4)
    for i in range(3):
        k_numerical = np.sqrt(eigenvals_0[i])
        k_analytical = solver.get_analytical_solution(0, i+1)
        u_2d = solver.create_2d_mode(eigenvals_0, eigenvecs_0, r_inner_0, 0, i, 'cos')
        modes_data.append({
            'mode': f'(0,{i+1})',
            'k_num': k_numerical,
            'k_ana': k_analytical,
            'u_2d': u_2d,
            'm': 0,
            'n': i+1
        })
    
    # m = 1 modes
    eigenvals_1, eigenvecs_1, r_inner_1 = solver.solve_eigenvalue_problem(m=1, n_modes=4)
    for i in range(3):
        k_numerical = np.sqrt(eigenvals_1[i])
        k_analytical = solver.get_analytical_solution(1, i+1)
        u_2d = solver.create_2d_mode(eigenvals_1, eigenvecs_1, r_inner_1, 1, i, 'cos')
        modes_data.append({
            'mode': f'(1,{i+1})',
            'k_num': k_numerical,
            'k_ana': k_analytical,
            'u_2d': u_2d,
            'm': 1,
            'n': i+1
        })
    
    # m = 2 modes
    eigenvals_2, eigenvecs_2, r_inner_2 = solver.solve_eigenvalue_problem(m=2, n_modes=3)
    for i in range(2):
        k_numerical = np.sqrt(eigenvals_2[i])
        k_analytical = solver.get_analytical_solution(2, i+1)
        u_2d = solver.create_2d_mode(eigenvals_2, eigenvecs_2, r_inner_2, 2, i, 'cos')
        modes_data.append({
            'mode': f'(2,{i+1})',
            'k_num': k_numerical,
            'k_ana': k_analytical,
            'u_2d': u_2d,
            'm': 2,
            'n': i+1
        })
    
    # Plot the first 8 modes
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.ravel()
    
    for idx in range(8):
        mode_data = modes_data[idx]
        u_2d = mode_data['u_2d']
        
        # Convert to Cartesian for plotting
        X = solver.X_mesh
        Y = solver.Y_mesh
        
        # Plot
        im = axes[idx].contourf(X, Y, u_2d, levels=20, cmap='RdBu_r')
        axes[idx].contour(X, Y, u_2d, levels=10, colors='black', alpha=0.3, linewidths=0.5)
        
        # Add boundary circle
        theta = np.linspace(0, 2*np.pi, 100)
        x_boundary = np.cos(theta)
        y_boundary = np.sin(theta)
        axes[idx].plot(x_boundary, y_boundary, 'k-', linewidth=2)
        
        axes[idx].set_aspect('equal')
        axes[idx].set_xlim(-1.1, 1.1)
        axes[idx].set_ylim(-1.1, 1.1)
        
        # Title with mode info
        m, n = mode_data['m'], mode_data['n']
        k_num = mode_data['k_num']
        k_ana = mode_data['k_ana']
        error = abs(k_num - k_ana) / k_ana * 100
        
        axes[idx].set_title(f'Mode ({m},{n})\nk = {k_num:.3f}\nError: {error:.2f}%', fontsize=10)
        axes[idx].set_xticks([])
        axes[idx].set_yticks([])
        
        # Add colorbar
        plt.colorbar(im, ax=axes[idx], shrink=0.8)
    
    plt.suptitle('Cylindrical Wave Equation Modes\n∇²u = -k²u, u(R,φ) = 0', fontsize=14)
    plt.tight_layout()
    plt.savefig('modes_cylindrical.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print eigenvalue comparison
    print("Eigenvalue Comparison (k values):")
    print("Mode\tNumerical\tAnalytical\tError (%)")
    print("-" * 45)
    for mode_data in modes_data[:8]:
        m, n = mode_data['m'], mode_data['n']
        k_num = mode_data['k_num']
        k_ana = mode_data['k_ana']
        error = abs(k_num - k_ana) / k_ana * 100
        print(f"({m},{n})\t{k_num:.6f}\t{k_ana:.6f}\t{error:.3f}")

def main():
    """
    Main function to demonstrate the cylindrical wave equation solver
    """
    print("2D Cylindrical Wave Equation Solver")
    print("=" * 40)
    print("Equation: ∇²u = -k²u")
    print("Boundary: u(R,φ) = 0")
    print("Domain: 0 ≤ r ≤ R, 0 ≤ φ ≤ 2π")
    print()
    
    # Plot geometry
    print("Generating geometry plots...")
    plot_geometry()
    
    # Solve and visualize modes
    print("Solving eigenvalue problems and generating mode plots...")
    solve_and_visualize_modes()
    
    print("\nSolution complete!")
    print("Generated files:")
    print("- geometry_cylindrical.png")
    print("- modes_cylindrical.png")

if __name__ == "__main__":
    main()
