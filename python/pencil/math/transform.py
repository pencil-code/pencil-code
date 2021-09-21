# transform.py
#
# Routines to transform coordinates and fields.
#
# Author:
# J. Aarnes (jorgenaarnes@gmail.com)
"""
Routines to transform coordinates and fields.
"""


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
