"""
This script assumes pc_start has already been run. Otherwise, we would need to hardcode the values of nx,ny,nz.
"""

import pencil as pc
import numpy as np
import h5py

dim = pc.read.dim()

a = np.zeros((3,dim.nz,dim.ny,dim.nx))

a[1,:,:,:] = 0.1
a[2,:,:,:] = 0.2
a[0,1,1,1] = 0.05
a[0,0,1,2] = 0.01

with h5py.File("forcing_cont.h5", mode='w') as f:
	f["forcing_cont/x"] = a[0]
	f["forcing_cont/y"] = a[1]
	f["forcing_cont/z"] = a[2]
