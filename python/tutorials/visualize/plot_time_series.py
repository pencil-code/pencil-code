# Import all the necessary libraries.

import matplotlib.pyplot as plt
import pencil as pc

# Read the time series.
ts = pc.read.ts()

# Plot the kinetic energies: fluid and particle.
plt.plot(ts.t, ts.uxmax, linestyle='-', color='r', linewidth=2,
         label=r'$u_x^{\rm max}$')
plt.plot(ts.t, ts.uymax, linestyle='', color='green', marker='x',
         markersize=12, label=r'$u_y^{\rm max}$')
plt.plot(ts.t, ts.uzmax, linestyle='--', color=(1, 0, 1), linewidth=2,
         label=r'$u_z^{\rm max}$')

# Add lables.
plt.xlabel(r'$t$')
plt.ylabel(r'$u^{\rm max}$')

# Add a legend.
plt.legend(loc=0, shadow=False, fancybox=False, numpoints=1)
leg = plt.gca().get_legend()
frame = leg.get_frame()
frame.set_facecolor('1.0')
leg.draw_frame(False)

# Save the figure.
plt.savefig('u_max.pdf')
