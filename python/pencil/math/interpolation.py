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

    if (interpolation == 'mean') or (interpolation == 'trilinear'):
        # Find the adjacent indices.
        i = (xyz[0]-oxyz[0])/dxyz[0]
        ii = np.array([int(np.floor(i)), int(np.ceil(i))])
        if i < 0:
            i = 0
        if i > nxyz[0]-1:
            i = nxyz[0]-1
        ii = np.array([int(np.floor(i)), int(np.ceil(i))])

        j = (xyz[1]-oxyz[1])/dxyz[1]
        jj = np.array([int(np.floor(j)), int(np.ceil(j))])
        if j < 0:
            j = 0
        if j > nxyz[1]-1:
            j = nxyz[1]-1
        jj = np.array([int(np.floor(j)), int(np.ceil(j))])

        k = (xyz[2]-oxyz[2])/dxyz[2]
        kk = np.array([int(np.floor(k)), int(np.ceil(k))])
        if k < 0:
            k = 0
        if k > nxyz[2]-1:
            k = nxyz[2]-1
        kk = np.array([int(np.floor(k)), int(np.ceil(k))])

    # Interpolate the field.
    if interpolation == 'mean':
        sub_field = field[:, :, :, [ii[0], ii[1]]]
        sub_field = sub_field[:, :, [jj[0], jj[1]], :]
        sub_field = sub_field[:, [kk[0], kk[1]], :, :]
        return np.mean(sub_field, axis=(1, 2, 3))

    if interpolation == 'trilinear':
        if ii[0] == ii[1]:
            w1 = np.array([1, 1])
        else:
            w1 = (i-ii[::-1])

        if jj[0] == jj[1]:
            w2 = np.array([1, 1])
        else:
            w2 = np.array([nxyz[1]-j, j-jj[0]])

        if kk[0] == kk[1]:
            w3 = np.array([1, 1])
        else:
            w3 = np.array([nxyz[2]-k, k-kk[0]])

        weight = abs(w3.reshape((2, 1, 1))*w2.reshape((1, 2, 1))*w1.reshape((1, 1, 2)))
        sub_field = field[:, :, :, [ii[0], ii[1]]]
        sub_field = sub_field[:, :, [jj[0], jj[1]], :]
        sub_field = sub_field[:, [kk[0], kk[1]], :, :]
        return np.sum(sub_field*weight, axis=(1, 2, 3))/np.sum(weight)

    # If the point lies outside the domain, return 0.
    if (ii[0] < -1) or (ii[1] > nxyz[0]) or (jj[0] < -1) or (jj[1] > nxyz[1]) \
        or (kk[0] < -1) or (kk[1] > nxyz[2]):
        return np.zeros([0, 0, 0])
