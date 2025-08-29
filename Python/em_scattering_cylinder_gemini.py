"""
2D Electromagnetic Scattering from a Perfectly Conducting Cylinder
=================================================================

A clear and modern implementation for solving the scattering of a plane
wave by a perfectly conducting cylinder using eigenfunction expansion.

Problem Setup:
- Incident plane wave with Ez polarization (TE mode)
- Perfectly conducting cylinder of radius 'b' at the origin
- Boundary condition: Total Ez field is zero at r = b

Mathematical Formulation:
- Incident field: Ez_inc = E0 * exp(ik*x)
- This is expanded using the Jacobi-Anger identity:
  exp(ikr*cos(φ)) = Σ [i^n * J_n(kr) * exp(inφ)]
- Scattered field: Represented by a series of outgoing cylindrical waves:
  Ez_scat = Σ [a_n * H_n^(1)(kr) * exp(inφ)]
- The coefficients 'a_n' are found by enforcing the boundary condition.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv, hankel1

class CylinderScatteringSolver:
    """
    Solves the EM scattering problem for a plane wave incident on a PEC cylinder.
    """
    def __init__(self, radius, wavelength, field_amplitude=1.0):
        """
        Initializes the solver with the physical parameters of the problem.

        Args:
            radius (float): The radius 'b' of the conducting cylinder.
            wavelength (float): The wavelength 'λ' of the incident plane wave.
            field_amplitude (float): The amplitude 'E0' of the incident field.
        """
        self.b = radius
        self.lambda_ = wavelength
        self.E0 = field_amplitude
        self.k = 2 * np.pi / self.lambda_
        self.kb = self.k * self.b

        # Determine the number of terms needed for the series expansion
        # This is a well-known rule of thumb for convergence.
        self.n_terms = int(2 * self.kb + 10)

        # Pre-calculate the scattering coefficients
        self.coefficients = self._calculate_coefficients()

    def _calculate_coefficients(self):
        """
        Calculates the scattering coefficients 'a_n' based on the boundary condition.
        For a PEC cylinder, the condition Ez_total(r=b) = 0 leads to:
        a_n = -E0 * i^n * J_n(kb) / H_n^(1)(kb)
        """
        coeffs = {}
        for n in range(-self.n_terms, self.n_terms + 1):
            jn_kb = jv(n, self.kb)
            hn_kb = hankel1(n, self.kb)

            # Avoid division by zero for numerically unstable Hankel functions
            if np.abs(hn_kb) > 1e-30:
                coeffs[n] = -self.E0 * (1j)**n * jn_kb / hn_kb
            else:
                coeffs[n] = 0
        return coeffs

    def get_fields(self, r, phi):
        """
        Calculates the incident, scattered, and total fields at given points.

        Args:
            r (np.ndarray): Radial coordinates of the grid points.
            phi (np.ndarray): Angular coordinates of the grid points.

        Returns:
            tuple: A tuple containing (Ez_inc, Ez_scat, Ez_total).
        """
        Ez_inc = np.zeros_like(r, dtype=complex)
        Ez_scat = np.zeros_like(r, dtype=complex)

        for n, an in self.coefficients.items():
            # Incident field using Jacobi-Anger expansion
            jn_kr = jv(n, self.k * r)
            Ez_inc += self.E0 * (1j)**n * jn_kr * np.exp(1j * n * phi)

            # Scattered field using Hankel functions for outgoing waves
            hn_kr = hankel1(n, self.k * r)
            Ez_scat += an * hn_kr * np.exp(1j * n * phi)

        # Total field is the sum, but zero inside the cylinder
        Ez_total = Ez_inc + Ez_scat
        Ez_total[r <= self.b] = 0
        
        # The scattered field should also be zero inside
        Ez_scat[r <= self.b] = 0

        return Ez_inc, Ez_scat, Ez_total

def visualize_fields(solver):
    """
    Creates a three-panel plot of the incident, scattered, and total fields.
    """
    # Set up the computational grid
    grid_size = 3 * solver.b
    n_points = 250
    x = np.linspace(-grid_size, grid_size, n_points)
    y = np.linspace(-grid_size, grid_size, n_points)
    X, Y = np.meshgrid(x, y)

    # Convert grid to cylindrical coordinates
    R = np.sqrt(X**2 + Y**2)
    PHI = np.arctan2(Y, X)

    print("Calculating fields on the grid...")
    Ez_inc, Ez_scat, Ez_total = solver.get_fields(R, PHI)
    print("Field calculation complete.")

    # Create the plots
    fig, axes = plt.subplots(1, 3, figsize=(20, 6.5))
    field_data = {
        'Incident Field': Ez_inc,
        'Scattered Field': Ez_scat,
        'Total Field': Ez_total
    }

    # Determine a common color scale for visual consistency
    max_val = np.max(np.abs(np.real(Ez_total)))
    norm = plt.Normalize(vmin=-max_val, vmax=max_val)
    cmap = 'RdBu_r'

    for ax, (title, data) in zip(axes, field_data.items()):
        im = ax.imshow(np.real(data), extent=[-grid_size, grid_size, -grid_size, grid_size],
                       origin='lower', cmap=cmap, norm=norm)
        
        # Add the cylinder boundary
        cylinder_boundary = plt.Circle((0, 0), solver.b, color='black', fill=False, linewidth=2)
        ax.add_artist(cylinder_boundary)
        
        # Fill the cylinder to represent the conductor
        cylinder_fill = plt.Circle((0, 0), solver.b, color='gray', alpha=0.6)
        ax.add_artist(cylinder_fill)

        ax.set_title(title, fontsize=16, fontweight='bold')
        ax.set_xlabel('x / b', fontsize=12)
        ax.set_ylabel('y / b', fontsize=12)
        ax.set_aspect('equal')
        ax.grid(True, linestyle='--', alpha=0.5)

    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.8, label='Re(Ez)')
    fig.suptitle(f'EM Scattering from a PEC Cylinder (Radius b = {solver.b}, kb = {solver.kb:.2f})',
                 fontsize=18, fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    # Save the figure
    filename = 'em_scattering_cylinder_gemini.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Figure saved as '{filename}'")
    
    plt.show()

def main():
    """
    Main function to set up and run the scattering simulation.
    """
    # --- Parameters ---
    cylinder_radius = 1.5  # meters
    incident_wavelength = 1.0  # meters
    # ------------------

    print("Initializing simulation...")
    solver = CylinderScatteringSolver(radius=cylinder_radius, wavelength=incident_wavelength)
    
    print(f"Cylinder Radius (b): {solver.b} m")
    print(f"Wavelength (λ): {solver.lambda_} m")
    print(f"Wave number (k): {solver.k:.2f} rad/m")
    print(f"Size parameter (kb): {solver.kb:.2f}")
    print("-" * 30)

    visualize_fields(solver)

if __name__ == "__main__":
    main()
