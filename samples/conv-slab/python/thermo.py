#!/usr/bin/python
# -*- encoding: utf8 -*-
"""
Calculate all relevant thermodynamical variables from lnrho and ss for box convection runs.

This is a Python translation of the IDL script: ../idl/thermo.pro

The script computes thermodynamic quantities from simulation variables:
- Pressure (pp)
- Sound speed squared (cs2)
- Temperature (TT)
- Initial profiles of temperature, density, and entropy for stratified convection

These are derived from the logarithmic density (lnrho) and specific entropy (ss)
using the equation of state and polytropic stratification profiles.
"""

import pencil as pc
import numpy as np
import matplotlib.pyplot as plt
import os

# Default sample directory, change the path to your working directory
#sdir = os.path.join(os.environ['PENCIL_HOME'],'sample/conv-slab')

# Create simulation object
sdir = 'illa/test/conv-slab'
sim = pc.sim.simulation(sdir)

# Read parameters from the simulation
# This reads both start.in and run.in parameters
print("Reading simulation parameters...")
#par = pc.read.param(datadir=sim.datadir, quiet=False)
par = sim.param

# Extract key thermodynamic parameters

cp = par.get('cp',1.0)
rho0 = par.get('rho0',1.0)
cs20 = par.get('cs0',1.0)**2
gamma = par.get('gamma',1.6666666)
gamma_m1 = gamma - 1.0
lnrho0 = np.log(rho0)

# Get z-coordinate information
z1 = par.get('z1', -0.5)
z2 = par.get('z2', 0.5)
zref = par.get('zref', 1.32)
z0 = par.get('xyz0', [-0.68])[-1]
Lz = par.get('Lz', 2.0)
ztop = z0 + Lz
gravz = par.get('gravz', -1.0)

# Get stratification parameters (polytropic indices)

mpoly0 = par.get('mpoly0', 1.0)
mpoly1 = par.get('mpoly1', 3.0)
mpoly2 = par.get('mpoly2', 0.0)
isothtop = par.get('isothtop', 1)

print(f"Thermodynamic parameters loaded:")
print(f"  cp={cp}, rho0={rho0}, cs20={cs20}")
print(f"  gamma={gamma}, z1={z1}, z2={z2}, zref={zref}")
print(f"  Stratification: mpoly0={mpoly0}, mpoly1={mpoly1}, mpoly2={mpoly2}")
print(f"  Isothermal top: {isothtop}")


# Function to initialize temperature, density, and entropy profiles
# This mimics the IDL thermo.pro routine
def initialize_stratification(z, z0=z0, z1=z1, z2=z2, ztop=ztop, zref=zref):
    """
    Calculate initial profiles of temperature, density, and entropy.
    
    Based on the piecewise polytropic stratification with:
    - Unstable layer (z < z1): polytropic index mpoly1
    - Stable layer (z1 < z < z2): polytropic index mpoly0
    - Top layer (z > z2): polytropic index mpoly2 (or isothermal if isothtop=1)
    
    Parameters
    ----------
    z : array
        Height array
    z0, z1, z2, ztop, zref : float
        Height boundaries and reference height
        
    Returns
    -------
    Tinit : array
        Initial temperature profile
    lnrhoinit : array
        Initial log-density profile
    ssinit : array
        Initial specific entropy profile
    """
    
    Tinit = np.zeros_like(z)
    lnrhoinit = np.zeros_like(z)
    ssinit = np.zeros_like(z)
    
    # Reference state values at zref (top of stable layer)
    Tref = cs20 / gamma_m1 / cp
    lnrhoref = np.log(rho0)
    ssref = 0.0
    
    # Upper layer (z >= z2): polytropic or isothermal
    top_mask = z >= z2
    if np.any(top_mask):
        zint = min(z2, ztop)  # interface height
        
        if isothtop == 0:  # polytropic top layer
            beta = gamma / gamma_m1 / cp * gravz / (mpoly2 + 1)
            Tinit[top_mask] = np.maximum(Tref + beta * (z[top_mask] - zref), 1.e-5 * Tref)
            lnrhoinit[top_mask] = lnrhoref + mpoly2 * np.log(Tinit[top_mask] / Tref)
            ssinit[top_mask] = (ssref + 
                               (1 - mpoly2 * (gamma - 1)) / gamma * 
                               np.log(Tinit[top_mask] / Tref))
            
            # Update reference state for transition to stable layer
            temp_ratio = 1.0 + beta * (zint - zref) / Tref
            lnrhoref = lnrhoref + mpoly2 * np.log(temp_ratio)
            ssref = ssref + (1 - mpoly2 * (gamma - 1)) / gamma * np.log(temp_ratio)
            Tref = Tref + beta * (zint - zref)
            
        else:  # isothermal top layer (isothtop == 1)
            beta = 0.0
            Tinit[top_mask] = Tref
            lnrhoinit[top_mask] = lnrhoref + gamma / gamma_m1 * gravz * (z[top_mask] - zref) / Tref
            ssinit[top_mask] = ssref - gravz * (z[top_mask] - zref) / Tref
            
            # Update reference state for transition
            lnrhoref = lnrhoref + gamma / gamma_m1 * gravz * (z2 - zref) / Tref
            ssref = ssref - gravz * (zint - zref) / Tref
            # Tref remains the same
    
    # Stable layer (z1 <= z < z2): polytropic with mpoly0
    stab_mask = (z <= z2) & (z >= z1)
    if np.any(stab_mask):
        beta = gamma / gamma_m1 / cp * gravz / (mpoly0 + 1)
        Tinit[stab_mask] = Tref + beta * (z[stab_mask] - zint)
        lnrhoinit[stab_mask] = lnrhoref + mpoly0 * np.log(Tinit[stab_mask] / Tref)
        ssinit[stab_mask] = (ssref + 
                            (1 - mpoly0 * (gamma - 1)) / gamma * 
                            np.log(Tinit[stab_mask] / Tref))
    
    # Unstable layer (z < z1): polytropic with mpoly1
    unstab_mask = z <= z1
    if np.any(unstab_mask):
        # Update reference state from stable layer
        lnrhoref = lnrhoref + mpoly0 * np.log(1.0 + beta * (z1 - z2) / Tref)
        ssref = (ssref + 
                (1 - mpoly0 * (gamma - 1)) / gamma * 
                np.log(1.0 + beta * (z1 - z2) / Tref))
        Tref = Tref + beta * (z1 - z2)
        
        # Now apply mpoly1 for unstable layer
        beta = gamma / gamma_m1 / cp * gravz / (mpoly1 + 1)
        Tinit[unstab_mask] = Tref + beta * (z[unstab_mask] - z1)
        lnrhoinit[unstab_mask] = lnrhoref + mpoly1 * np.log(Tinit[unstab_mask] / Tref)
        ssinit[unstab_mask] = (ssref + 
                              (1 - mpoly1 * (gamma - 1)) / gamma * 
                              np.log(Tinit[unstab_mask] / Tref))
    
    return Tinit, lnrhoinit, ssinit

def plot_profiles(temperature, lnrho, ss, title,savefig=False):
    # Plot the stratification profiles
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)

    axes[0].plot(temperature, z, 'b-', linewidth=2)
    axes[0].axhline(y=z1, color='r', linestyle='--', alpha=0.5, label='z1 (layer boundary)')
    axes[0].axhline(y=z2, color='g', linestyle='--', alpha=0.5, label='z2 (layer boundary)')
    axes[0].set_xlabel('Temperature')
    axes[0].set_ylabel('Height (z)')
    axes[0].set_title('Initial Temperature Profile')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(np.exp(lnrho), z, 'g-', linewidth=2)
    axes[1].axhline(y=z1, color='r', linestyle='--', alpha=0.5, label='z1')
    axes[1].axhline(y=z2, color='g', linestyle='--', alpha=0.5, label='z2')
    axes[1].set_xlabel('Density (ρ)')
    axes[1].set_ylabel('Height (z)')
    axes[1].set_title('Initial Density Profile')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(ss, z, 'r-', linewidth=2)
    axes[2].axhline(y=z1, color='r', linestyle='--', alpha=0.5, label='z1')
    axes[2].axhline(y=z2, color='g', linestyle='--', alpha=0.5, label='z2')
    axes[2].set_xlabel('Specific Entropy (s)')
    axes[2].set_ylabel('Height (z)')
    axes[2].set_title('Initial Entropy Profile')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)

    plt.suptitle(title, fontsize=14, fontweight='bold')

    

    if savefig:
            # Directory to save figures. Set to 'fig' by default
            figdir = os.path.join(sdir,'fig')
            # Check that the figs directory exist and create it if it does not
            if not os.path.isdir(figdir):
                    os.makedirs(figdir)
                    print(f"Created directory 'fig' inside {sdir}")


            savefig = os.path.join(figdir,'stratification_profiles.png')
            plt.savefig(savefig, format='png', bbox_inches='tight')
            print(f"Figure saved as {savefig}")

    plt.show()


# Read simulation data (optional: only if you want to compute thermodynamics from actual snapshots)
# Uncomment the following to read actual simulation data:
# 
# print("\nReading simulation data...")
# var = pc.read.var(sim=sim)  # Read full 3D snapshot
# lnrho = var.lnrho
# ss = var.ss
# 
# # Compute thermodynamic variables from actual snapshot data
# pp = cs20 * rho0 / gamma * np.exp(gamma * (ss / cp + lnrho - lnrho0))
# cs2 = cs20 * np.exp(gamma * ss / cp + gamma_m1 * (lnrho - lnrho0))
# TT = cs2 / gamma_m1 / cp

# For now, compute initial stratification profiles only
print("\nComputing initial stratification profiles...")
# Create z-array from grid information
dim = sim.dim
z = np.linspace(z0, ztop, dim.nz)

Tinit, lnrhoinit, ssinit = initialize_stratification(z)

print(f"  Temperature range: {Tinit.min():.4f} to {Tinit.max():.4f}")
print(f"  Log-density range: {lnrhoinit.min():.4f} to {lnrhoinit.max():.4f}")
print(f"  Entropy range: {ssinit.min():.4f} to {ssinit.max():.4f}")

title = 'Initial Stratification Profiles'
# Savefig controls if you want to save the result as a figure
plot_profiles(Tinit, lnrhoinit, ssinit, title, savefig=True)

print("\nThermodynamic analysis complete!")