# plot_animation.py
'''
Plots some quantities as animation.
'''

import os
import numpy as np
import pencilnew as pn

run_dir = 'n256_tw2_dz16_oo01_nu1e-5_cl_peri_peri'
datadir = os.path.join(run_dir, 'data')

# Read the data.
#slices = pn.read.slices(datadir=datadir, field=['oo1'], extension=['yz'])

# Prepare the plot.
width = 5
height = 7
#plt.rc('text', usetex=True)
#plt.rc('font', family='arial')
plt.rc("figure.subplot", left=0.15)
plt.rc("figure.subplot", right=0.75)
plt.rc("figure.subplot", bottom=0.2)
plt.rc("figure.subplot", top=0.95)
fig = plt.figure(figsize=(width, height))
ax = fig.add_subplot(111)

# Plot the data.
im = plt.imshow(slices.yz.oo1[0, :, :], extent=[-5.5, 5.5, -16, 16], origin='lower', interpolation='nearest', cmap='bwr', vmin=-0.05, vmax=0.05)

# Make plot pretty.
#plt.xticks([0, 2, 4, 6, 8, 10])
#plt.yticks([0, 2, 4, 6, 8, 10])
plt.tick_params(axis='both', which='major', length=8, labelsize=20)
plt.tick_params(axis='both', which='minor', length=4, labelsize=20)

plt.xlabel(r'$y$', fontsize=25)
plt.ylabel(r'$z$', fontsize=25)
#plt.xlabel(r'$x [{\rm Mm}]$', fontsize=25)
#plt.ylabel(r'$y [{\rm Mm}]$', fontsize=25)
ax.xaxis.label.set_fontname('serif')

#plt.xlim(-0.5, 10.5)
#plt.ylim(-0.5, 10.5)

for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontname('serif')
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontname('serif')

for label in ax.xaxis.get_ticklabels():
    label.set_position((0, -0.03))
for label in ax.yaxis.get_ticklabels():
    label.set_position((-0.03, 0))

# Add a color bar.
cb = plt.colorbar()
cb.set_label(r'$\omega_r$', fontsize = 25)
cbytick_obj = plt.getp(cb.ax.axes, 'yticklabels')
plt.setp(cbytick_obj, fontsize = 15, family = 'serif')

plt.draw()

fig.savefig('animation0000.png', dpi=200)
# Loop over all times.
for t_idx, t in enumerate(slices.t[1:]):
    im.set_data(slices.yz.oo1[t_idx, :, :])
    fig.savefig('animation{0:04}.png'.format(t_idx), dpi=200)

