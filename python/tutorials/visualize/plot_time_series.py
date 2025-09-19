# Import all the necessary libraries.

import numpy as np
import matplotlib.pyplot as plt
import pencil as pc

# Read the time series.
ts = pc.read.ts()

# Plot the kinetic energies: fluid and particle.
plt.plot(ts.t, ts.ekin, linestyle='-', color='r', linewidth=2, label=r'$E_{\rm kin}$')
plt.plot(ts.t, ts.ekinp, linestyle='--', color='green', linewidth=2, label=r'$E_{\rm kin}^{\rm p}$')

# Add lables.
plt.xlabel(r'$t$')
plt.ylabel(r'$\rho$')

# Add a legend.
plt.legend(loc=0, shadow=False, fancybox=False, numpoints=1)
leg = plt.gca().get_legend()
frame = leg.get_frame()
frame.set_facecolor('1.0')
leg.draw_frame(False)

# Save the figure.
plt.savefig('rho.pdf')
