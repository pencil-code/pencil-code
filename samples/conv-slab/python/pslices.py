#!/usr/bin/python
# -*- encoding: utf8 -*-
# Monitoring numerical stability



import pencil as pc
import numpy as np
import pylab as plt
import os.path as op

plt.rc('text', usetex=True)
plt.rc('font', family='arial')
plt.rcParams['font.size'] = 20

plt.xticks(fontsize=20, family='serif')
plt.yticks(fontsize=20, family='serif')

sdir = 'illa/test/conv-slab/'

sim = pc.sim.simulation(sdir)


# Read time series
ts = pc.read.ts(sim=sim)

# Plot the urms

rows = 1
cols = 1

fig, axs = plt.subplots(nrows=rows, ncols=cols, figsize=(12, 6), constrained_layout = True)

# Read slices
# Before reading the slices, one needs to ensemble them:
# make read_videofiles
# ./src/read_videofiles.x
# select uu1,uu2,uu3

# Read 3D field data at a specific time

slices = pc.read.slices()

    

# Access velocity field (3D array with shape [nt, nx, ny])
uux_xy = slices.xy.uu1
  
    

# Access density field
rho = np.exp(slices.xz.lnrho)


# Plot vertical slice of vertical velocity to visualize convective plumes

plt.contourf(uu[:, :, 64, 2])  # Vertical component at mid-height
plt.colorbar(label='Vertical velocity')
plt.title('Convective Plumes: Vertical Velocity')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
    

# Plot horizontal slice of density
plt.figure(figsize=(10, 8))
plt.contourf(rho[:, :, 64])
plt.colorbar(label='Density')
plt.title('Density Structure at Mid-height')
plt.xlabel('x')
plt.ylabel('y')
plt.show()


savefig = False

if savefig:
        savefig = op.join(sdir, 'fig','pvar.png')
        plt.savefig(savefig, format='png', bbox_inches='tight')
        print(f"Figure saved as {savefig}")



plt.show()