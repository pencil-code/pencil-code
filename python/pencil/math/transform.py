# transform.py
#
# Routines to transform coordinates and fields.
#
# Author:
# J. Aarnes (jorgenaarnes@gmail.com)
# S. Candelaresi (iomsn1@gmail.com)
"""
Routines to transform coordinates and fields.
"""


def coordinate_transformation(x, y, z, xyz_from='', xyz_to=''):
    """
    Tranform coordinates between coordinate systems.
    Returns a tuple (x, y, z).

    call signature:

    coordinate_transformation(x, y, z, xyz_from='', xyz_to='')

    Keyword arguments:

    *x, y, z*: 1d ndarray
      Input coordinates.

    *xyz_from*: str
      Origin coordinate system: 'cartesian', 'cylindrical' and 'spherical'.

    *xyz_to*: str
      Destination coordinate system: 'cartesian', 'cylindrical' and 'spherical'.
    """

    import numpy as np

    x, y, z = np.meshgrid(x, y, z, indexing='ij')

    if xyz_from == 'cartesian':
        if xyz_to == 'cylindrical':
            r = np.sqrt(x**2 + y**2)
            phi = np.arctan2(y, x)
            z = z
            return (r, phi, z)
        if xyz_to == 'spherical':
            r = np.sqrt(x**2 + y**2 + z**2)
            theta = np.arctan2(z, r)
            phi = np.arctan2(y, x)
            return (r, theta, phi)
    if xyz_from == 'cylindrical':
        if xyz_to == 'cartesian':
            return (x*np.cos(y), x*np.sin(y), z)
        if xyz_to == 'spherical':
            r = np.sqrt(x**2 + z**2)
            theta = np.arctan2(x, z)
            phi = y
            return (r, theta, phi)
    if xyz_from == 'spherical':
        if xyz_to == 'cartesian':
            return (x*np.cos(y)*np.cos(z), x*np.sin(y)*np.cos(z), x*np.sin(z))
        if xyz_to == 'cylindrical':
            r = x*np.cos(z)
            theta = x*np.sin(z)
            phi = y
            return (r, theta, phi)


#def field_transformation():
#    pass


def pospolar2cart(r, th):
    """
    Transform from polar to cartesian coordinates
    """

    import numpy as np

    rsize = r.size
    thsize = th.size

    x = np.empty([thsize + 1, rsize])
    y = np.empty([thsize + 1, rsize])
    for i in range(rsize):
        for j in range(thsize):
            x[j, i] = r[i] * np.cos(th[j])
            y[j, i] = r[i] * np.sin(th[j])
        x[thsize, i] = r[i] * np.cos(th[0])
        y[thsize, i] = r[i] * np.sin(th[0])

    return x, y


def velpolar2cart(vr, vth, r, th, zcoord=0):
    """
    Transform from polar velocities to cartesian velocities

    call signature:

    velpolar2cart(vr, vth, r, th, zcoord=0)

    Keyword arguments:

    *vr*:
      3D array of velocities in radial direction

    *vth*:
      3D array of velocities in theta direction

    *r*
      coordinates in radial direction

    *th*
      coordinates in theta direction

    *zcoord*
      plane in z-dir
    """

    import numpy as np

    rsize = r.size
    thsize = th.size
    vrsize = vr.size
    vthsize = vth.size
    if vrsize != vthsize:
        print("ERROR: vr and vth not equal in size!")
        return [], []
    vx = np.empty([thsize + 1, rsize])
    vy = np.empty([thsize + 1, rsize])
    for i in range(rsize):
        for j in range(thsize):
            vx[j, i] = vr[zcoord, j, i] * np.cos(th[j]) - vth[zcoord, j, i] * np.sin(
                th[j]
            )
            vy[j, i] = vr[zcoord, j, i] * np.sin(th[j]) + vth[zcoord, j, i] * np.cos(
                th[j]
            )
        vx[thsize, i] = vr[zcoord, 0, i] * np.cos(th[0]) - vth[zcoord, 0, i] * np.sin(
            th[0]
        )
        vy[thsize, i] = vr[zcoord, 0, i] * np.sin(th[0]) + vth[zcoord, 0, i] * np.cos(
            th[0]
        )

    return vx, vy
