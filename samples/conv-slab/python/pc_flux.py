#!/usr/bin/python
# -*- encoding: utf8 -*-
"""
Analyze and visualize 1-D energy fluxes from convection simulations (conv-slab).

This is a Python translation of the IDL script: ../idl/pc_fluxz.pro

The script reads vertical flux data from xy-averaged files and plots various
energy flux components:
- Kinetic energy flux (fkin)
- Radiative flux (frad)
- Convective flux (fconv)
- Cooling flux (fcool) - computed from cooling rate by integration
- Turbulent flux (fturb)

These fluxes are analyzed over time to understand energy transport in the
convective layer.
"""

import pencil as pc
import numpy as np
import matplotlib.pyplot as plt
import os

# Default sample directory, change the path to your working directory
#sdir = os.path.join(os.environ['PENCIL_HOME'],'sample/conv-slab')
workdir = 'illa/test/conv-slab'
sdir = os.path.join(os.environ['PENCIL_HOME'],workdir)
# Create simulation object
#sdir = 'illa/test/conv-slab'
sim = pc.sim.simulation(sdir)

# Read parameters from the simulation
print("Reading simulation parameters...")
par = sim.param
dim = sim.dim
param = pc.read.param(datadir=os.path.join(sdir,'data/'))

# Extract grid information
print(f"Grid dimensions: nx={dim.nx}, ny={dim.ny}, nz={dim.nz}")

# Read xy-averaged data which contains flux information
print("\nReading xy-averaged data...")
aver = pc.read.aver(simdir=sdir, datadir=os.path.join(sdir,'data/'),param=param,quiet=False)

# Get time and z-coordinate information
t = aver.t
z0 = par.get('xyz0', [-0.68])[-1]
Lz = par.get('Lz', 2.0)
ztop = z0 + Lz
z = np.linspace(z0, ztop, dim.nz)

# Get the number of timesteps and z-grid points
nt = len(t)
nz = dim.nzgrid

print(f"  Time range: {t[0]:.2f} to {t[-1]:.2f}")
print(f"  Vertical range: {z[0]:.2f} to {z[-1]:.2f}")
print(f"  Number of timesteps: {nt}")
print(f"  Number of z-levels: {nz}")

# Extract flux components if they exist
# Keys in the averages file are stored as attributes
fluxes = {}

# Try to access various flux components
if hasattr(aver.xy, 'fmasszmz'):
    fluxes['fmass'] = aver.xy.fmasszmz  # Mass flux
    print("  ✓ fmassz (mass flux)")

if hasattr(aver.xy, 'fkinzmz'):
    fluxes['fkin'] = aver.xy.fkinzmz   # Kinetic energy flux
    print("  ✓ fkinz (kinetic energy flux)")

if hasattr(aver.xy, 'fradz'):
    fluxes['frad'] = aver.xy.fradz   # Radiative flux
    print("  ✓ fradz (radiative flux)")

if hasattr(aver.xy, 'fconvz'):
    fluxes['fconv'] = aver.xy.fconvz # Convective flux
    print("  ✓ fconvz (convective flux)")

if hasattr(aver.xy, 'fturbz'):
    fluxes['fturb'] = aver.xy.fturbz # Turbulent flux
    print("  ✓ fturbz (turbulent flux)")

if hasattr(aver.xy, 'dcoolz'):
    dcool = aver.xy.dcoolz            # Cooling rate
    print("  ✓ dcoolz (cooling rate)")
    
    # Compute cooling flux by vertical integration
    # fcool = -∫ dcool dz (negative because cooling removes energy)
    fcool = np.zeros_like(dcool)
    for it in range(1, nt):
        # Integrate cooling rate vertically using trapezoid rule
        fcool[it,:] = -np.trapezoid(dcool[it,:], z)
    fluxes['fcool'] = fcool
    print("  ✓ Computed fcool (cooling flux via integration)")


def plot_fluxes_final_time(z, fluxes, t, savefig=False):
    """
    Plot all available flux components at the final timestep.
    
    Parameters
    ----------
    z : array
        Vertical coordinate array
    fluxes : dict
        Dictionary of flux arrays with shape (nz, nt)
    t : array
        Time array
    savefig : bool
        If True, save figure to fig/flux_final_time.png
    """
    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
    
    # Get color mapping for different fluxes
    colors = {
        'fkin': '#FF6B6B',    # red
        'frad': '#4ECDC4',    # teal
        'fconv': '#FFE66D',   # yellow
        'fcool': '#95E1D3',   # mint
        'fturb': '#A8DADC',   # light blue
        'fmass': '#C7CEEA'    # lavender
    }
    
    #nt = len(t)
    final_time = t[-1]
    
    # Plot each available flux at final time
    for flux_name, flux_data in fluxes.items():
        if flux_data is not None:
            ax.plot(z, flux_data[-1], label=flux_name, 
                   color=colors.get(flux_name, None), linewidth=2.5)
    
    ax.axhline(y=0, color='k', linestyle='-', alpha=0.2, linewidth=0.8)
    ax.set_xlabel('Height (z)', fontsize=12)
    ax.set_ylabel('Flux', fontsize=12)
    ax.set_title(f'Energy Fluxes at Final Time\n(t = {final_time:.2f})', fontsize=13, fontweight='bold')
    ax.legend(loc='best', framealpha=0.9)
    ax.grid(True, alpha=0.3)
    
    if savefig:
        figdir = os.path.join(sdir, 'fig')
        if not os.path.isdir(figdir):
            os.makedirs(figdir)
            print(f"Created directory 'fig' inside {sdir}")
        
        figpath = os.path.join(figdir, 'flux_final_time.png')
        fig.savefig(figpath, format='png', dpi=150, bbox_inches='tight')
        print(f"Figure saved as {figpath}")
    
    plt.show()


def plot_fluxes_time_evolution(z, fluxes, t, time_indices=None, ncols=3, savefig=False):
    """
    Create a grid of plots showing flux evolution at different times.
    
    Parameters
    ----------
    z : array
        Vertical coordinate array
    fluxes : dict
        Dictionary of flux arrays with shape (nt, nz)
    t : array
        Time array
    time_indices : array, optional
        Specific time indices to plot. If None, evenly spaced times are chosen.
    ncols : int
        Number of columns in subplot grid
    savefig : bool
        If True, save figure to fig/flux_evolution.png
    """
    #nt = len(t)
    
    # Select time indices to plot
    if time_indices is None:
        if nt <= 6:
            time_indices = np.arange(nt)
        else:
            time_indices = np.linspace(0, nt-1, 6, dtype=int)
    
    nrows = int(np.ceil(len(time_indices) / ncols))
    
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, 
                            figsize=(15, 4*nrows), constrained_layout=True)
    axs = axs.flatten()
    
    colors = {
        'fkin': '#FF6B6B',
        'frad': '#4ECDC4',
        'fconv': '#FFE66D',
        'fcool': '#95E1D3',
        'fturb': '#A8DADC',
        'fmass': '#C7CEEA'
    }
    
    for plot_idx, time_idx in enumerate(time_indices):
        ax = axs[plot_idx]
        
        # Plot all flux components at this time
        for flux_name, flux_data in fluxes.items():
            if flux_data is not None:
                ax.plot(z, flux_data[time_idx], label=flux_name,
                       color=colors.get(flux_name, None), linewidth=2)
        
        ax.axhline(y=0, color='k', linestyle='-', alpha=0.2, linewidth=0.8)
        ax.set_xlabel('Height (z)', fontsize=10)
        ax.set_ylabel('Flux', fontsize=10)
        ax.set_title(f't = {t[time_idx]:.2f}', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best', fontsize=9)
    
    # Hide unused subplots
    for plot_idx in range(len(time_indices), len(axs)):
        axs[plot_idx].set_visible(False)
    
    fig.suptitle('Energy Flux Evolution Over Time', fontsize=14, fontweight='bold')

    if savefig:
        figdir = os.path.join(sdir, 'fig')
        if not os.path.isdir(figdir):
            os.makedirs(figdir)
            print(f"Created directory 'fig' inside {sdir}")
        
        figpath = os.path.join(figdir, 'flux_evolution.png')
        fig.savefig(figpath, format='png', dpi=150, bbox_inches='tight')
        print(f"Figure saved as {figpath}")
    
    plt.show()
    
    


def plot_total_flux_time_series(z, fluxes, t, savefig=False):
    """
    Plot the time evolution of total flux at different heights.
    
    Parameters
    ----------
    z : array
        Vertical coordinate array
    fluxes : dict
        Dictionary of flux arrays with shape (nt, nz)
    t : array
        Time array
    savefig : bool
        If True, save figure to fig/total_flux_timeseries.png
    """
    # Compute total flux at each time and height
    ftot = np.zeros((len(t), len(z)))
    for flux_data in fluxes.values():
        if flux_data is not None:
            ftot += flux_data
    
    fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)
    
    # Plot total flux at different heights
    #nz = len(z)
    z_indices = np.linspace(0, nz-1, 5, dtype=int)  # Sample 5 heights
    
    for z_idx in z_indices:
        ax.plot(t, ftot[:,z_idx], label=f'z = {z[z_idx]:.2f}', linewidth=2)
    
    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel('Total Flux', fontsize=12)
    ax.set_title('Total Energy Flux Time Series at Different Heights', fontsize=13, fontweight='bold')
    ax.legend(loc='best', framealpha=0.9)
    ax.grid(True, alpha=0.3)
    
    if savefig:
        figdir = os.path.join(sdir, 'fig')
        if not os.path.isdir(figdir):
            os.makedirs(figdir)
            print(f"Created directory 'fig' inside {sdir}")
        
        figpath = os.path.join(figdir, 'total_flux_timeseries.png')
        fig.savefig(figpath, format='png', dpi=150, bbox_inches='tight')
        print(f"Figure saved as {figpath}")
    plt.show()
    


# Generate plots
print("\n" + "="*60)
print("Generating flux visualizations...")
print("="*60)

# Plot 1: All fluxes at final time
plot_fluxes_final_time(z, fluxes, t, savefig=True)

# Plot 2: Time evolution at multiple times
plot_fluxes_time_evolution(z, fluxes, t, savefig=True)

# Plot 3: Total flux time series at different heights
plot_total_flux_time_series(z, fluxes, t, savefig=True)

# Print summary statistics
print("\n" + "="*60)
print("Flux Statistics")
print("="*60)

for flux_name, flux_data in fluxes.items():
    if flux_data is not None:
        print(f"\n{flux_name}:")
        print(f"  Min: {flux_data.min():.6e}")
        print(f"  Max: {flux_data.max():.6e}")
        print(f"  Mean: {flux_data.mean():.6e}")
        print(f"  At final time - Min: {flux_data[:, -1].min():.6e}, Max: {flux_data[:, -1].max():.6e}")

print("\n" + "="*60)
print("Analysis complete!")
print("="*60)

plt.show()
