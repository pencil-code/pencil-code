#!/usr/bin/python
# -*- encoding: utf8 -*-
# Monitoring numerical stability



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

# Select the number of rows and columns
rows = 1
cols = 1

# Create figure
fig, axs = plt.subplots(nrows=rows, ncols=cols, figsize=(12, 6), constrained_layout = True)


# Plot maximum velocity to see if the flow is becoming turbulent

axs.plot(ts.t, ts.dtc, label='dtc (compressibility)',ls='--', linewidth=2)
axs.plot(ts.t, ts.dtu, label='dtu (velocity)',ls='-.', linewidth=2)
axs.plot(ts.t, ts.dtnu, label='dtnu (viscosity)',ls=':', linewidth=2)
axs.plot(ts.t, ts.dtchi, label='dtchi (heat conduction)', ls='-',linewidth=2)
axs
axs.set_xlabel('Time')
axs.set_ylabel('Timestep constraint')
axs.set_title('Timestep Limiting Factors')
fig.legend(loc='upper right',bbox_to_anchor=(1, 0.8))
axs.grid(True)


# Chose if you want to save your figures
savefig = True


if savefig:
        # Directory to save figures. Set to 'fig' by default
        figdir = os.path.join(sdir,'fig')
        # Check that the figs directory exist and create it if it does not
        if not os.path.isdir(figdir):
                os.makedirs(figdir)
                print(f"Created directory 'fig' inside {sdir}")


        savefig = os.path.join(figdir,'numerical_stability.png')
        plt.savefig(savefig, format='png', bbox_inches='tight')
        print(f"Figure saved as {savefig}")



plt.show()