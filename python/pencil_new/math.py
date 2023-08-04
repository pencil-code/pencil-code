#!/usr/bin/python3
# Last Modification: $Id$
#=======================================================================
# math.py
#
# Facilities for mathematical operations.
#
# Chao-Chin Yang, 2016-06-29
#=======================================================================
def curl(v1, v2, v3, grid):
    """Takes the curl of a vector field in Cartesian coordinates.

    Positional Arguments
        v1, v2, v3
            Components of the vector field.
        grid
            Grid geometry returned by read.grid(trim=True).
    """
    # Author: Chao-Chin Yang
    # Created: 2016-06-29
    import numpy as np

    ng = 3	# Number of ghost cells.

    c1 = (grid.dy_1[np.newaxis,:,np.newaxis] * _deriv1(v3, axis=1)[ng:-ng,:,ng:-ng] -
          grid.dz_1[np.newaxis,np.newaxis,:] * _deriv1(v2, axis=2)[ng:-ng,ng:-ng,:])
    c2 = (grid.dz_1[np.newaxis,np.newaxis,:] * _deriv1(v1, axis=2)[ng:-ng,ng:-ng,:] -
          grid.dx_1[:,np.newaxis,np.newaxis] * _deriv1(v3, axis=0)[:,ng:-ng,ng:-ng])
    c3 = (grid.dx_1[:,np.newaxis,np.newaxis] * _deriv1(v2, axis=0)[:,ng:-ng,ng:-ng] -
          grid.dy_1[np.newaxis,:,np.newaxis] * _deriv1(v1, axis=1)[ng:-ng,:,ng:-ng])

    return c1, c2, c3
#=======================================================================
def _deriv1(a, axis=-1):
    """Takes the first derivative over the specified axis, assuming the
    grid spacing is unity.

    Positional Argument
        a
            Array to be operated on.

    Keyword Arguments
        axis
            Axis of a to be opearted on.
    """
    # Author: Chao-Chin Yang
    # Created: 2016-06-29
    b = a.view() if axis == -1 else a.swapaxes(-1, axis)
    db = 0.75 * (b[...,4:-2] - b[...,2:-4]) - 0.15 * (b[...,5:-1] - b[...,1:-5]) + (1/60) * (b[...,6:] - b[...,:-6])
    return db if axis == -1 else db.swapaxes(axis,-1)
