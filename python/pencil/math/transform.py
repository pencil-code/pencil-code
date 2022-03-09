# transform.py
#
# Routines to transform coordinates and fields.
"""
Routines to transform coordinates and fields.
"""


def coordinate_transformation(x, y, z, xyz_from="", xyz_to=""):
    """
    coordinate_transformation(x, y, z, xyz_from='', xyz_to='')

    Tranform coordinates between coordinate systems.
    Returns a tuple (x, y, z).

    Parameters
    ----------
    x, y, z : 1d ndarray
        Input coordinates.

    xyz_from : string
        Origin coordinate system: 'cartesian', 'cylindrical' and 'spherical'.

    xyz_to : string
        Destination coordinate system: 'cartesian', 'cylindrical' and 'spherical'.

    Returns
    -------
    Tuple of 3 ndarrays containing the transformed coordinates.
    """

    import numpy as np

    x, y, z = np.meshgrid(x, y, z, indexing="ij")

    if xyz_from == "cartesian":
        if xyz_to == "cylindrical":
            r = np.sqrt(x ** 2 + y ** 2)
            phi = np.arctan2(y, x)
            z = z
            return (r, phi, z)
        if xyz_to == "spherical":
            r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
            theta = np.arctan2(z, r)
            phi = np.arctan2(y, x)
            return (r, theta, phi)
    if xyz_from == "cylindrical":
        if xyz_to == "cartesian":
            return (x * np.cos(y), x * np.sin(y), z)
        if xyz_to == "spherical":
            r = np.sqrt(x ** 2 + z ** 2)
            theta = np.arctan2(x, z)
            phi = y
            return (r, theta, phi)
    if xyz_from == "spherical":
        if xyz_to == "cartesian":
            return (x * np.cos(y) * np.cos(z), x * np.sin(y) * np.cos(z), x * np.sin(z))
        if xyz_to == "cylindrical":
            r = x * np.cos(z)
            theta = x * np.sin(z)
            phi = y
            return (r, theta, phi)


def vector_field_transformation(field, x, y, z, xyz_from="", xyz_to=""):
    """
    vector_field_transformation(field, x, y, z, xyz_from='', xyz_to='')

    Transform a vector field from one coordinate system to another.

    Parameters
    ----------
    field : 4d ndarray
        Vector field.

    x, y, z : 1d ndarrays
        Input coordinates of the 'from' system.

    xyz_from : string
        Origin coordinate system: 'cartesian', 'cylindrical' and 'spherical'.

    xyz_to : string
        Destination coordinate system: 'cartesian', 'cylindrical' and 'spherical'.

    Returns
    -------
    ndarray with the transformed vector field.
    """

    import numpy as np

    (u, v, w) = coordinate_transformation(x, y, z, xyz_from=xyz_from, xyz_to=xyz_to)
    u = np.swapaxes(u, 0, 2)
    v = np.swapaxes(v, 0, 2)
    w = np.swapaxes(w, 0, 2)

    if xyz_from == "cartesian":
        if xyz_to == "cylindrical":
            field_r = np.cos(v) * field[0] + np.sin(v) * field[1]
            field_phi = -np.cos(v) * field[0] + np.cos(v) * field[1]
            field_z = field[2]
            return np.array([field_r, field_phi, field_z])
        if xyz_to == "spherical":
            field_r = (
                np.sin(v) * np.cos(w) * field[0]
                + np.sin(v) * np.sin(w) * field[1]
                + np.cos(w) * field[2]
            )
            field_theta = (
                np.cos(v) * np.cos(w) * field[0]
                + np.cos(v) * np.sin(w) * field[1]
                - np.sin(w) * field[2]
            )
            field_phi = -np.sin(w) * field[0] + np.cos(w) * field[1]
            return np.array([field_r, field_theta, field_phi])
    if xyz_from == "cylindrical":
        if xyz_to == "cartesian":
            field_x = np.cos(y) * field[0] - np.sin(y) * field[1]
            field_y = np.cos(y) * field[0] + np.cos(y) * field[1]
            field_x = field[2]
            return np.array([field_x, field_y, field_z])
        if xyz_to == "spherical":
            field_r = np.sin(v) * field[0] + np.cos(v) * field[2]
            field_theta = np.cos(v) * field[0] - np.sin(v) * field[2]
            field_phi = field[1]
            return np.array([field_r, field_theta, field_phi])
    if xyz_from == "spherical":
        if xyz_to == "cartesian":
            field_x = (
                np.sin(y) * np.cos(z) * field[0]
                + np.cos(y) * np.cos(z) * field[1]
                + -np.sin(z) * field[2]
            )
            field_y = (
                np.sin(y) * np.sin(z) * field[0]
                + np.cos(y) * np.sin(z) * field[1]
                + np.cos(z) * field[2]
            )
            field_z = np.cos(z) * field[0] - np.sin(z) * field[1]
            return np.array([field_x, field_y, field_z])
        if xyz_to == "cylindrical":
            field_r = np.sin(y) * field[0] + np.cos(y) * field[1]
            field_phi = field[2]
            field_z = np.cos(y) * field[0] - np.sin(y) * field[1]
            return np.array([field_phi, field_theta, field_z])


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
