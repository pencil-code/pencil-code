# interpolation.py
# Written by Simon Candelaresi (iomsn1@gmail.com)

"""
Interpolation routines for scalar and vector fields.
"""

import numpy as np


def vec_int(xyz, var, field, interpolation='trilinear'):
    """
    Interpolates the field around this position.

    call signature:

        vecInt(xyz, var, field, interpolation = 'weighted')

    Keyword arguments:

    *xyz*:
      Position vector around which will be interpolated.

    *var*:
        The var object from the read_var routine.

    *field*:
      Vector field to be interpolated.

    *interpolation*:
      Interpolation of the vector field.
      'mean': Takes the mean of the adjacent grid points.
      'trilinear': Weights the adjacent grid points according to their distance.
    """

    Ox = var.x[0]
    Oy = var.y[0]
    Oz = var.z[0]
    nx = len(var.x)
    ny = len(var.y)
    nz = len(var.z)

    # Find the adjacent indices.
    i = (xyz[0]-Ox)/var.dx
    if i < 0:
        i = 0
    if i > nx-1:
        i = nx-1
    ii = np.array([int(np.floor(i)), \
                    int(np.ceil(i))])

    j = (xyz[1]-Oy)/var.dy
    if j < 0:
        j = 0
    if j > ny-1:
        j = ny-1
    jj = np.array([int(np.floor(j)), \
                    int(np.ceil(j))])

    k = (xyz[2]-Oz)/var.dz
    if k < 0:
        k = 0
    if k > nz-1:
        k = nz-1
    kk = np.array([int(np.floor(k)), \
                    int(np.ceil(k))])

    # Interpolate the field.
    field = np.swapaxes(field, 1, 3)
    if interpolation == 'mean':
        return np.mean(field[:, ii[0]:ii[1]+1, jj[0]:jj[1]+1, kk[0]:kk[1]+1],
                       axis=(1, 2, 3))
    if interpolation == 'trilinear':
        if ii[0] == ii[1]:
            w1 = np.array([1, 1])
        else:
            w1 = (i-ii[::-1])
        if jj[0] == jj[1]:
            w2 = np.array([1, 1])
        else:
            w2 = (j-jj[::-1])
        if kk[0] == kk[1]:
            w3 = np.array([1, 1])
        else:
            w3 = (k-kk[::-1])
        weight = abs(w1.reshape((2, 1, 1))*w2.reshape((1, 2, 1))*\
                 w3.reshape((1, 1, 2)))
        return np.sum(field[:, ii[0]:ii[1]+1, jj[0]:jj[1]+1,
                            kk[0]:kk[1]+1]*weight, axis=(1, 2, 3)) \
                            /np.sum(weight)


def vec_int_no_var(xx, field, params, interpolation='weighted'):
    """
    Interpolates the field around this position without the need of
    the full var object.

    call signature:

        vec_int_no_var(xx, field, p, interpolation = 'weighted')

    Keyword arguments:

    *xx*:
      Position vector around which will be interpolated.

    *field*:
      Vector field to be interpolated.

    *params*:
      Parameter object.

    *interpolation*:
      Interpolation of the vector field.
      'mean': takes the mean of the adjacent grid point.
      'weighted': weights the adjacent grid points according to their distance.
    """

    # Find the adjacent indices.
    i = (xx[0]-params.Ox)/params.dx
    if i < 0:
        i = 0
    if i > params.nx-1:
        i = params.nx-1
    ii = np.array([int(np.floor(i)), int(np.ceil(i))])

    j = (xx[1]-params.Oy)/params.dy
    if j < 0:
        j = 0
    if j > params.ny-1:
        j = params.ny-1
    jj = np.array([int(np.floor(j)), int(np.ceil(j))])

    k = (xx[2]-params.Oz)/params.dz
    if k < 0:
        k = 0
    if k > params.nz-1:
        k = params.nz-1
    kk = np.array([int(np.floor(k)), int(np.ceil(k))])

    field = np.swapaxes(field, 1, 3)
    # Interpolate the field.
    if interpolation == 'mean':
        return np.mean(field[:, ii[0]:ii[1]+1, jj[0]:jj[1]+1, kk[0]:kk[1]+1],
                       axis=(1, 2, 3))
    if interpolation == 'weighted':
        if ii[0] == ii[1]:
            w1 = np.array([1, 1])
        else:
            w1 = (i-ii[::-1])
        if jj[0] == jj[1]:
            w2 = np.array([1, 1])
        else:
            w2 = (j-jj[::-1])
        if kk[0] == kk[1]:
            w3 = np.array([1, 1])
        else: w3 = (k-kk[::-1])
        weight = abs(w1.reshape((2, 1, 1))*w2.reshape((1, 2, 1))*\
                 w3.reshape((1, 1, 2)))
        return np.sum(field[:, ii[0]:ii[1]+1, jj[0]:jj[1]+1, kk[0]:kk[1]+1]*weight,
                      axis=(1, 2, 3))/np.sum(weight)
