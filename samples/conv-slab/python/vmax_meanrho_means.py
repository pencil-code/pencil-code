#!/usr/bin/python
# -*- encoding: utf8 -*-
# Tracking maximum velocity and mean quantities



import pencil as pc
import numpy as np
import pylab as plt
import os

plt.rc('text', usetex=True)
plt.rc('font', family='arial')
plt.rcParams['font.size'] = 20

plt.xticks(fontsize=20, family='serif')
plt.yticks(fontsize=20, family='serif')

# Default sample directory, change the path to your working directory
#sdir = os.path.join(os.environ['PENCIL_HOME'],'sample/conv-slab')

sdir = 'illa/test/conv-slab/'

# Create simulation object
sim = pc.sim.simulation(sdir)

# Read time series
ts = pc.read.ts(sim=sim)

# Select number of rows and columns for figure
rows = 1
cols = 3

# Create figure
fig, axs = plt.subplots(nrows=rows, ncols=cols, figsize=(16, 4),constrained_layout = True)

# Plot maximum velocity to see if the flow is becoming turbulent

axs[0].plot(ts.t, ts.umax, 'r-')
axs[0].set_xlabel('Time')
axs[0].set_ylabel('Max velocity')
axs[0].set_title('Peak Velocity')
axs[0].grid(True)

# Monitor mean density changes

axs[1].plot(ts.t, ts.rhom, 'g-')
axs[1].set_xlabel('Time')
axs[1].set_ylabel('Mean density')
axs[1].set_title('Mean Density Evolution')
axs[1].grid(True)

# Monitor mean entropy (important for stratification)

axs[2].plot(ts.t, ts.ssm, 'b-')
axs[2].set_xlabel('Time')
axs[2].set_ylabel('Mean entropy')
axs[2].set_title('Mean Entropy Evolution')
axs[2].grid(True)


# Chose if you want to save your figures
savefig = True


if savefig:
        # Directory to save figures. Set to 'fig' by default
        figdir = os.path.join(sdir,'fig')
        # Check that the figs directory exist and create it if it does not
        if not os.path.isdir(figdir):
                os.makedirs(figdir)
                print(f"Created directory 'fig' inside {sdir}")


        savefig = os.path.join(figdir,'umax-meanpho-means.png')
        plt.savefig(savefig, format='png', bbox_inches='tight')
        print(f"Figure saved as {savefig}")



plt.show()