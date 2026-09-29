#!/usr/bin/env python

import pencil as pc
import matplotlib.pyplot as plt
import numpy as np

sl = pc.read.slices()
sim = pc.sim.Simulation(".", quiet=True)

im = plt.contourf(sim.grid.x, sim.grid.z, sl.xz.pp[-1] - sl.xz.pp[0], cmap='bwr', levels=101)
im.set_clim(np.array([-1,1])*np.max(np.abs(im.get_clim())))
plt.colorbar()

plt.show()
