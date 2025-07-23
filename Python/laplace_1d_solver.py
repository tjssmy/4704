"""
1D Finite Difference Solution for Laplace's Equation
====================================================

Solves the 1D Laplace equation: d²u/dx² = 0
with boundary conditions:
- du/dx(0) = 2 (Neumann BC at x=0)
- u(L) = 0 (Dirichlet BC at x=L)

The analytical solution for this problem is: u(x) = 2*(x-L)
"""

import numpy as np
import matplotlib.pyplot as plt

def solve_laplace_1d(L=1.0, n_points=51):
    """
    Solve 1D Laplace equation using finite differences
    
    Parameters:
    -----------
    L : float
        Length of the domain (default: 1.0)
    n_points : int
        Number of grid points (default: 51)
    
    Returns:
    --------
    x : numpy array
        Grid points
    u : numpy array
        Solution values
    """
    
    # Create grid
    dx = L / (n_points - 1)
    x = np.linspace(0, L, n_points)
    
    # Initialize coefficient matrix A and right-hand side vector b
    A = np.zeros((n_points, n_points))
    b = np.zeros(n_points)
    
    # For interior points (i = 1 to n-2), use central difference:
    # (u[i+1] - 2*u[i] + u[i-1]) / dx² = 0
    # This gives: u[i-1] - 2*u[i] + u[i+1] = 0
    for i in range(1, n_points - 1):
        A[i, i-1] = 1.0
        A[i, i] = -2.0
        A[i, i+1] = 1.0
        b[i] = 0.0
    
    # Boundary condition at x = 0: du/dx(0) = 2
    # Using forward difference: (u[1] - u[0]) / dx = 2
    # This gives: u[1] - u[0] = 2*dx
    # Rearranging: -u[0] + u[1] = 2*dx
    A[0, 0] = -1.0
    A[0, 1] = 1.0
    b[0] = 2.0 * dx
    
    # Boundary condition at x = L: u(L) = 0
    A[n_points-1, n_points-1] = 1.0
    b[n_points-1] = 0.0
    
    # Solve the linear system
    u = np.linalg.solve(A, b)
    
    return x, u

def analytical_solution(x, L):
    """
    Analytical solution: u(x) = 2*(x-L)
    """
    return 2.0 * (x - L)

def plot_results(x, u_numerical, u_analytical):
    """
    Plot numerical and analytical solutions
    """
    plt.figure(figsize=(10, 6))
    
    plt.subplot(1, 2, 1)
    plt.plot(x, u_numerical, 'bo-', label='Numerical', markersize=4)
    plt.plot(x, u_analytical, 'r-', label='Analytical', linewidth=2)
    plt.xlabel('x')
    plt.ylabel('u(x)')
    plt.title('1D Laplace Equation Solution')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.subplot(1, 2, 2)
    error = np.abs(u_numerical - u_analytical)
    plt.plot(x, error, 'g-', linewidth=2)
    plt.xlabel('x')
    plt.ylabel('|Error|')
    plt.title('Absolute Error')
    plt.grid(True, alpha=0.3)
    plt.yscale('log')
    
    plt.tight_layout()
    
    # Save the figure
    plt.savefig('laplace_1d_solution.png', dpi=300, bbox_inches='tight')
    print("Figure saved as 'laplace_1d_solution.png'")
    
    plt.show()

def main():
    """
    Main function to solve and visualize the problem
    """
    # Problem parameters
    L = 1.0  # Domain length
    n_points = 51  # Number of grid points
    
    print("1D Finite Difference Solution for Laplace's Equation")
    print("=" * 55)
    print(f"Domain: [0, {L}]")
    print(f"Number of grid points: {n_points}")
    print(f"Grid spacing (dx): {L/(n_points-1):.6f}")
    print("\nBoundary Conditions:")
    print("- du/dx(0) = 2 (Neumann)")
    print("- u(L) = 0 (Dirichlet)")
    print()
    
    # Solve numerically
    x, u_numerical = solve_laplace_1d(L, n_points)
    
    # Calculate analytical solution
    u_analytical = analytical_solution(x, L)
    
    # Calculate errors
    max_error = np.max(np.abs(u_numerical - u_analytical))
    rms_error = np.sqrt(np.mean((u_numerical - u_analytical)**2))
    
    print("Results:")
    print(f"Maximum absolute error: {max_error:.2e}")
    print(f"RMS error: {rms_error:.2e}")
    print()
    
    # Print some values for verification
    print("Sample values:")
    print("x\t\tNumerical\tAnalytical\tError")
    print("-" * 50)
    for i in range(0, len(x), len(x)//5):
        print(f"{x[i]:.3f}\t\t{u_numerical[i]:.6f}\t{u_analytical[i]:.6f}\t\t{abs(u_numerical[i] - u_analytical[i]):.2e}")
    
    # Plot results
    plot_results(x, u_numerical, u_analytical)
    
    # Verify boundary conditions
    print(f"\nBoundary condition verification:")
    print(f"du/dx(0) ≈ {(u_numerical[1] - u_numerical[0])/(x[1] - x[0]):.6f} (should be 2.0)")
    print(f"u(L) = {u_numerical[-1]:.6f} (should be 0.0)")

if __name__ == "__main__":
    main()
