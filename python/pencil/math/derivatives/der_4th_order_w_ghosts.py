# der_4th_order_w_ghosts.py
#
# Calculate 4th order spatial derivatives with ghost zones assumed to be in
# the arrays.
#
# Author: Simon Candelaresi (iomsn1@gmail.com).
"""
4th Order derivatives. Currently only equidistant grids are supported.
"""


def xder3_4th(f, dx):
    """
    xder3_4th(f, dx)

    Compute the 3rd order derivative, 4th order accurate in x.

    Parameters
    ----------
    f : ndarray
        Array for which to compute the derivative.

    dx : float
        Grid-spacing in x.
    """

    import numpy as np

    if f.ndim != 3 and f.ndim != 4:
        print("{0} dimension arrays not handled.".format(str(f.ndim)))
        raise ValueError

    dx3 = 1.0 / (8.0 * dx ** 3.0)
    dfdx = np.zeros_like(f)
    l1 = 3
    l2 = f.shape[-1] - 3

    if l2 > l1:
        dfdx[..., l1:l2] = dx3 * (
            +13.0 * (f[..., l1 - 1 : l2 - 1] - f[..., l1 + 1 : l2 + 1])
            - 8.0 * (f[..., l1 - 2 : l2 - 2] - f[..., l1 + 2 : l2 + 2])
            + 1.0 * (f[..., l1 - 3 : l2 - 3] - f[..., l1 + 3 : l2 + 3])
        )
    else:
        dfdx = 0.0

    return dfdx


def yder3_4th(f, dy):
    """
    yder3_4th(f, dy)

    Compute the 3rd order derivative, 4th order accurate in y.

    Parameters
    ----------
    f : ndarray
        Array for which to compute the derivative.

    dy : float
        Grid-spacing in y.
    """

    import numpy as np

    if f.ndim != 3 and f.ndim != 4:
        print("{0} dimension arrays not handled.".format(str(f.ndim)))
        raise ValueError

    dy3 = 1.0 / (8.0 * dy ** 3.0)
    dfdy = np.zeros_like(f)
    m1 = 3
    m2 = f.shape[-2] - 3

    if m2 > m1:
        dfdy[..., m1:m2, :] = dy3 * (
            +13.0 * (f[..., m1 - 1 : m2 - 1, :] - f[..., m1 + 1 : m2 + 1, :])
            - 8.0 * (f[..., m1 - 2 : m2 - 2, :] - f[..., m1 + 2 : m2 + 2, :])
            + 1.0 * (f[..., m1 - 3 : m2 - 3, :] - f[..., m1 + 3 : m2 + 3, :])
        )
    else:
        dfdy = 0.0

    return dfdy


def zder3_4th(f, dz):
    """
    zder3_4th(f, dz)

    Compute the 3rd order derivative, 4th order accurate in z.

    Parameters
    ----------
    f : ndarray
        Array for which to compute the derivative.

    dz : float
        Grid-spacing in z.
    """

    import numpy as np

    if f.ndim != 3 and f.ndim != 4:
        print("{0} dimension arrays not handled.".format(str(f.ndim)))
        raise ValueError

    dz3 = 1.0 / (8.0 * dz ** 3.0)
    dfdz = np.zeros_like(f)
    n1 = 3
    n2 = f.shape[-3] - 3
    if n2 > n1:
        dfdz[..., n1:n2, :, :] = dz3 * (
            +13.0 * (f[..., n1 - 1 : n2 - 1, :, :] - f[..., n1 + 1 : n2 + 1, :, :])
            - 8.0 * (f[..., n1 - 2 : n2 - 2, :, :] - f[..., n1 + 2 : n2 + 2, :, :])
            + 1.0 * (f[..., n1 - 3 : n2 - 3, :, :] - f[..., n1 + 3 : n2 + 3, :, :])
        )
    else:
        dfdz = 0.0

    return dfdz
