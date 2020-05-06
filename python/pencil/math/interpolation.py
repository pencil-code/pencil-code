# interpolation.py
#
# Interpolation routines for scalar and vector fields.
#
# Written by Simon Candelaresi (iomsn1@gmail.com)
"""
Interpolation routines for scalar and vector fields.
"""

def vec_int(xyz, field, dxyz, oxyz, nxyz, interpolation='trilinear'):
    """
    Interpolates the field around position xyz.

    call signature:

        vec_int(xyz, field, dxyz, oxyz, nxyz, interpolation='trilinear')

    Keyword arguments:

    *xyz*:
      Position vector around which will be interpolated.

    *field*:
      Vector field to be interpolated with shape [nz, ny, nx].

    *dxyz*:
      Array with the three deltas.

    *oxyz*:
      Array with the position of the origin.

    *nxyz*:
      Number of grid points in each direction.

    *interpolation*:
      Interpolation method. Can be 'mean' or 'trilinear'.
    """

    import numpy as np

    # Find the adjacent indices.
    i = (xyz[0] - oxyz[0])/dxyz[0]
    if i < 0:
        i = 0
    if i > nxyz[0] - 1:
        i = nxyz[0] - 1
    ii = np.array([int(np.floor(i)), int(np.ceil(i))])

    j = (xyz[1] - oxyz[1])/dxyz[1]
    if j < 0:
        j = 0
    if j > nxyz[1] - 1:
        j = nxyz[1] - 1
    jj = np.array([int(np.floor(j)), int(np.ceil(j))])

    k = (xyz[2] - oxyz[2])/dxyz[2]
    if k < 0:
        k = 0
    if k > nxyz[2] - 1:
        k = nxyz[2] - 1
    kk = np.array([int(np.floor(k)), int(np.ceil(k))])

    # Interpolate the field.
    if interpolation == 'mean':
        return np.mean(field[:, kk[0]:kk[1]+1, jj[0]:jj[1]+1, ii[0]:ii[1]+1],
                       axis=(1, 2, 3))
    elif interpolation == 'trilinear':
        if ii[0] == ii[1]:
            w1 = np.array([1, 1])
        else:
            w1 = i - ii[::-1]
        if jj[0] == jj[1]:
            w2 = np.array([1, 1])
        else:
            w2 = j - jj[::-1]
        if kk[0] == kk[1]:
            w3 = np.array([1, 1])
        else:
            w3 = k - kk[::-1]
        weight = abs(w1.reshape((2, 1, 1))*w2.reshape((1, 2, 1))*\
                 w3.reshape((1, 1, 2)))
        return np.sum(field[:, kk[0]:kk[1]+1, jj[0]:jj[1]+1,
                            ii[0]:ii[1]+1]*weight, axis=(1, 2, 3)) \
                            /np.sum(weight)
    else:
        print('Error: cannot find interpolation method {0}.'.format(interpolation))
        return -1
