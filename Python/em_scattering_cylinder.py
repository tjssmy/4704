"""
2D Electromagnetic Scattering from a Perfectly Conducting Cylinder
=================================================================

Solves the scattering of a plane wave by a perfectly conducting cylinder
using eigenfunction expansion in cylindrical coordinates.

Problem Setup:
- Incident plane wave with Ez polarization
- Perfectly conducting cylinder of radius b at origin
- 2D problem (no z-variation, TE mode)
- Boundary condition: Ez = 0 at r = b

Mathematical Formulation:
- Incident field: Ez_inc = E0 * exp(ik*x) = E0 * exp(ik*r*cos(φ))
- Scattered field: Eigenfunction expansion using Hankel functions
- Total field: Ez_total = Ez_inc + Ez_scattered
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv, yv, hankel1, iv, kv
import warnings
warnings.filterwarnings('ignore')

class EMScattering:
    def __init__(self, k=2*np.pi, b=1.0, E0=1.0):
        """
        Initialize electromagnetic scattering problem
        
        Parameters:
        -----------
        k : float
            Wave number (2π/λ)
        b : float
            Cylinder radius
        E0 : float
            Incident field amplitude
        """
        self.k = k
        self.b = b
        self.E0 = E0
        self.kb = k * b  # Normalized frequency parameter
        
        # Calculate scattering coefficients
        self.coefficients = self.calculate_scattering_coefficients()
        
    def calculate_scattering_coefficients(self):
        """
        Calculate the scattering coefficients a_n using boundary conditions
        
        For a perfectly conducting cylinder: Ez = 0 at r = b
        This gives: a_n = -J_n(kb) / H_n^(1)(kb)
        """
        # Number of terms in expansion (converged for kb)
        n_terms = max(20, int(2 * self.kb + 10))
        
        coefficients = np.zeros(2*n_terms + 1, dtype=complex)
        
        for n in range(-n_terms, n_terms + 1):
            n_idx = n + n_terms  # Array index
            
            # Bessel function of first kind
            jn_kb = jv(n, self.kb)
            
            # Hankel function of first kind
            hn_kb = hankel1(n, self.kb)
            
            # Scattering coefficient
            if abs(hn_kb) > 1e-15:  # Avoid division by zero
                coefficients[n_idx] = -jn_kb / hn_kb
            else:
                coefficients[n_idx] = 0.0
                
        return coefficients
    
    def incident_field(self, r, phi):
        """
        Calculate incident field: Ez_inc = E0 * exp(ik*r*cos(φ))
        Using Jacobi-Anger expansion: exp(ikr*cos(φ)) = Σ i^n * J_n(kr) * exp(inφ)
        """
        # Grid dimensions
        if np.isscalar(r):
            r = np.array([r])
            phi = np.array([phi])
        
        Ez_inc = np.zeros_like(r, dtype=complex)
        n_terms = max(20, int(2 * self.kb + 10))
        
        for n in range(-n_terms, n_terms + 1):
            # Jacobi-Anger expansion coefficient
            coeff = (1j)**n
            
            # Bessel function
            jn_kr = jv(n, self.k * r)
            
            # Add to expansion
            Ez_inc += coeff * jn_kr * np.exp(1j * n * phi)
            
        return self.E0 * Ez_inc
    
    def scattered_field(self, r, phi):
        """
        Calculate scattered field using eigenfunction expansion
        Ez_scat = Σ a_n * H_n^(1)(kr) * exp(inφ)
        """
        # Only calculate for r > b (outside cylinder)
        Ez_scat = np.zeros_like(r, dtype=complex)
        
        # Mask for points outside cylinder
        outside_mask = r > self.b
        
        if np.any(outside_mask):
            n_terms = max(20, int(2 * self.kb + 10))
            
            for n in range(-n_terms, n_terms + 1):
                n_idx = n + n_terms
                an = self.coefficients[n_idx]
                
                # Hankel function of first kind
                hn_kr = hankel1(n, self.k * r[outside_mask])
                
                # Add to expansion
                Ez_scat[outside_mask] += an * hn_kr * np.exp(1j * n * phi[outside_mask])
                
        return Ez_scat
    
    def total_field(self, r, phi):
        """
        Calculate total field: Ez_total = Ez_inc + Ez_scat (outside cylinder)
        Ez_total = 0 (inside cylinder - perfect conductor)
        """
        Ez_total = np.zeros_like(r, dtype=complex)
        
        # Outside cylinder: incident + scattered
        outside_mask = r > self.b
        if np.any(outside_mask):
            Ez_inc = self.incident_field(r[outside_mask], phi[outside_mask])
            Ez_scat = self.scattered_field(r[outside_mask], phi[outside_mask])
            Ez_total[outside_mask] = Ez_inc + Ez_scat
            
        # Inside cylinder: Ez = 0 (perfect conductor)
        # Already initialized to zero
        
        return Ez_total

def create_scattering_plots():
    """
    Create three subplot figure showing incident, scattered, and total fields
    """
    # Problem parameters
    wavelength = 1.0
    k = 2 * np.pi / wavelength
    b = 0.8  # Cylinder radius
    E0 = 1.0
    
    # Create scattering solver
    em_scatter = EMScattering(k=k, b=b, E0=E0)
    
    print(f"EM Scattering Analysis:")
    print(f"Wavelength λ = {wavelength:.2f}")
    print(f"Wave number k = {k:.2f}")
    print(f"Cylinder radius b = {b:.2f}")
    print(f"Normalized frequency kb = {em_scatter.kb:.2f}")
    print()
    
    # Create computational grid
    x_max = 3.0
    y_max = 3.0
    nx, ny = 200, 200
    
    x = np.linspace(-x_max, x_max, nx)
    y = np.linspace(-y_max, y_max, ny)
    X, Y = np.meshgrid(x, y)
    
    # Convert to cylindrical coordinates
    R = np.sqrt(X**2 + Y**2)
    PHI = np.arctan2(Y, X)
    
    # Avoid singularity at origin
    R = np.maximum(R, 1e-10)
    
    # Calculate fields
    print("Calculating incident field...")
    Ez_inc = em_scatter.incident_field(R, PHI)
    
    print("Calculating scattered field...")
    Ez_scat = em_scatter.scattered_field(R, PHI)
    
    print("Calculating total field...")
    Ez_total = em_scatter.total_field(R, PHI)
    
    # Create figure with three subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Common plotting parameters
    levels = 30
    extent = [-x_max, x_max, -y_max, y_max]
    
    # Determine common color scale
    max_field = max(np.max(np.abs(np.real(Ez_inc))), 
                   np.max(np.abs(np.real(Ez_scat))), 
                   np.max(np.abs(np.real(Ez_total))))
    vmin, vmax = -max_field*0.8, max_field*0.8
    
    # Plot 1: Incident field
    ax1 = axes[0]
    im1 = ax1.contourf(X, Y, np.real(Ez_inc), levels=levels, 
                       cmap='RdBu_r', vmin=vmin, vmax=vmax, extent=extent)
    ax1.contour(X, Y, np.real(Ez_inc), levels=15, colors='black', alpha=0.3, linewidths=0.5)
    
    # Add cylinder boundary
    theta = np.linspace(0, 2*np.pi, 100)
    x_circle = b * np.cos(theta)
    y_circle = b * np.sin(theta)
    ax1.plot(x_circle, y_circle, 'k-', linewidth=3, label='Conducting cylinder')
    ax1.fill(x_circle, y_circle, color='gray', alpha=0.7)
    
    ax1.set_xlim(-x_max, x_max)
    ax1.set_ylim(-y_max, y_max)
    ax1.set_aspect('equal')
    ax1.set_xlabel('x/λ', fontsize=12)
    ax1.set_ylabel('y/λ', fontsize=12)
    ax1.set_title('Incident Field\nRe[Ez_inc]', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Add wave direction arrow
    ax1.arrow(-2.5, 2.5, 0.8, 0, head_width=0.15, head_length=0.2, 
              fc='red', ec='red', linewidth=2)
    ax1.text(-2.0, 2.8, 'k', fontsize=12, color='red', fontweight='bold')
    
    plt.colorbar(im1, ax=ax1, shrink=0.8, label='Ez field')
    
    # Plot 2: Scattered field
    ax2 = axes[1]
    im2 = ax2.contourf(X, Y, np.real(Ez_scat), levels=levels, 
                       cmap='RdBu_r', vmin=vmin, vmax=vmax, extent=extent)
    ax2.contour(X, Y, np.real(Ez_scat), levels=15, colors='black', alpha=0.3, linewidths=0.5)
    
    # Add cylinder boundary
    ax2.plot(x_circle, y_circle, 'k-', linewidth=3)
    ax2.fill(x_circle, y_circle, color='gray', alpha=0.7)
    
    ax2.set_xlim(-x_max, x_max)
    ax2.set_ylim(-y_max, y_max)
    ax2.set_aspect('equal')
    ax2.set_xlabel('x/λ', fontsize=12)
    ax2.set_ylabel('y/λ', fontsize=12)
    ax2.set_title('Scattered Field\nRe[Ez_scat]', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    plt.colorbar(im2, ax=ax2, shrink=0.8, label='Ez field')
    
    # Plot 3: Total field
    ax3 = axes[2]
    im3 = ax3.contourf(X, Y, np.real(Ez_total), levels=levels, 
                       cmap='RdBu_r', vmin=vmin, vmax=vmax, extent=extent)
    ax3.contour(X, Y, np.real(Ez_total), levels=15, colors='black', alpha=0.3, linewidths=0.5)
    
    # Add cylinder boundary
    ax3.plot(x_circle, y_circle, 'k-', linewidth=3)
    ax3.fill(x_circle, y_circle, color='gray', alpha=0.7)
    
    ax3.set_xlim(-x_max, x_max)
    ax3.set_ylim(-y_max, y_max)
    ax3.set_aspect('equal')
    ax3.set_xlabel('x/λ', fontsize=12)
    ax3.set_ylabel('y/λ', fontsize=12)
    ax3.set_title('Total Field\nRe[Ez_total] = Re[Ez_inc + Ez_scat]', fontsize=14, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    
    plt.colorbar(im3, ax=ax3, shrink=0.8, label='Ez field')
    
    # Main title
    plt.suptitle(f'EM Scattering: Plane Wave on Conducting Cylinder (kb = {em_scatter.kb:.2f})', 
                 fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    
    # Save figure
    plt.savefig('em_scattering_fields.png', dpi=300, bbox_inches='tight')
    print("Saved figure: em_scattering_fields.png")
    plt.show()
    
    # Print some analysis
    print("\nField Analysis:")
    print("-" * 40)
    
    # Calculate scattering cross section (simplified)
    sigma_normalized = 0
    n_terms = min(10, len(em_scatter.coefficients)//2)
    for n in range(-n_terms, n_terms+1):
        n_idx = n + len(em_scatter.coefficients)//2
        an = em_scatter.coefficients[n_idx]
        sigma_normalized += np.abs(an)**2
    
    sigma_normalized *= 4 / k
    
    print(f"Geometric cross section: {2*b:.3f}λ")
    print(f"Scattering cross section: {sigma_normalized:.3f}λ")
    print(f"Efficiency factor: {sigma_normalized/(2*b):.3f}")
    
    # Resonance information
    if em_scatter.kb > 1:
        print(f"\nThis is in the resonance/optical regime (kb = {em_scatter.kb:.2f} > 1)")
        print("Multiple scattering lobes and interference patterns are visible")
    else:
        print(f"\nThis is in the Rayleigh regime (kb = {em_scatter.kb:.2f} < 1)")
        print("Scattering is predominantly forward/backward")

def main():
    """
    Main function to run the electromagnetic scattering analysis
    """
    print("2D Electromagnetic Scattering Analysis")
    print("=" * 50)
    print("Problem: Plane wave scattering by perfectly conducting cylinder")
    print("Method: Eigenfunction expansion with Hankel functions")
    print("Polarization: Ez (TE mode)")
    print()
    
    # Create and display plots
    create_scattering_plots()
    
    print("\nAnalysis complete!")

if __name__ == "__main__":
    main()
