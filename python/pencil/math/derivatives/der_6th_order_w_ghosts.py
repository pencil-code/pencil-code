# der_6th_order_w_ghosts.py
#
# Calculate 6th order spatial derivatives with ghost zones assumed to be in
# the arrays.
"""
6th Order derivatives. Currently only equidistant grids are supported.

FAG: added 5th derivative and filled derivative ghosts (except corners)
     with periodic boundary values. TODO: include proper boundary conditions
     for non-internal ghosts. Add separate get_ghosts routine to handle this.
"""


def xder_6th(f, dx):
    """
    xder_6th(f, dx)

    Compute the 1st order derivative, 6th order accurate in x.

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

    dx2 = 1.0 / (60.0 * dx)
    dfdx = np.zeros_like(f)
    l1 = 3
    l2 = f.shape[-1] - 3

    if l2 > l1:
        dfdx[..., l1:l2] = dx2 * (
            45.0 * (f[..., l1 + 1 : l2 + 1] - f[..., l1 - 1 : l2 - 1])
            - 9.0 * (f[..., l1 + 2 : l2 + 2] - f[..., l1 - 2 : l2 - 2])
            + (f[..., l1 + 3 : l2 + 3] - f[..., l1 - 3 : l2 - 3])
        )
        dfdx[..., :l1] = dfdx[..., l2 - 3 : l2]
        dfdx[..., l2:] = dfdx[..., l1 : l1 + 3]
    else:
        dfdx = 0.0

    return dfdx


def yder_6th(f, dy):
    """
    yder_6th(f, dy)

    Compute the 1st order derivative, 6th order accurate in y.

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

    dy2 = 1.0 / (60.0 * dy)
    dfdy = np.zeros_like(f)
    m1 = 3
    m2 = f.shape[-2] - 3

    if m2 > m1:
        dfdy[..., m1:m2, :] = dy2 * (
            45.0 * (f[..., m1 + 1 : m2 + 1, :] - f[..., m1 - 1 : m2 - 1, :])
            - 9.0 * (f[..., m1 + 2 : m2 + 2, :] - f[..., m1 - 2 : m2 - 2, :])
            + (f[..., m1 + 3 : m2 + 3, :] - f[..., m1 - 3 : m2 - 3, :])
        )
        dfdy[..., :m1, :] = dfdy[..., m2 - 3 : m2, :]
        dfdy[..., m2:, :] = dfdy[..., m1 : m1 + 3, :]
    else:
        dfdy = 0.0

    return dfdy


def zder_6th(f, dz, run2D=False):
    """
    zder_6th(f, dz)

    Compute the 1st order derivative, 6th order accurate in z.

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

    dz2 = 1.0 / (60.0 * dz)
    dfdz = np.zeros_like(f)
    n1 = 3

    if run2D:
        n2 = f.shape[1] - 3
    else:
        n2 = f.shape[-3] - 3

    if n2 > n1:
        if run2D:
            # Case f[..., z, x] or f[..., z, y].
            dfdz[..., n1:n2, :] = dz2 * (
                45.0 * (f[..., n1 + 1 : n2 + 1, :] - f[..., n1 - 1 : n2 - 1, :])
                - 9.0 * (f[..., n1 + 2 : n2 + 2, :] - f[..., n1 - 2 : n2 - 2, :])
                + (f[..., n1 + 3 : n2 + 3, :] - f[..., n1 - 3 : n2 - 3, :])
            )
            dfdz[..., :n1, :] = dfdz[..., n2 - 3 : n2, :]
            dfdz[..., n2:, :] = dfdz[..., n1 : n1 + 3, :]
        else:
            # Case f[...,z,y,x].
            dfdz[..., n1:n2, :, :] = dz2 * (
                45.0 * (f[..., n1 + 1 : n2 + 1, :, :] - f[..., n1 - 1 : n2 - 1, :, :])
                - 9.0 * (f[..., n1 + 2 : n2 + 2, :, :] - f[..., n1 - 2 : n2 - 2, :, :])
                + (f[..., n1 + 3 : n2 + 3, :, :] - f[..., n1 - 3 : n2 - 3, :, :])
            )
            dfdz[..., :n1, :, :] = dfdz[..., n2 - 3 : n2, :, :]
            dfdz[..., n2:, :, :] = dfdz[..., n1 : n1 + 3, :, :]
    else:
        dfdz = 0

    return dfdz


def xder2_6th(f, dx):
    """
    xder_6th(f, dx)

    Compute the 2nd order derivative, 6th order accurate in x.

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

    dx2 = 1.0 / (180.0 * dx ** 2.0)
    dfdx = np.zeros_like(f)
    l1 = 3
    l2 = f.shape[-1] - 3

    if l2 > l1:
        dfdx[..., l1:l2] = dx2 * (
            -490.0 * f[..., l1:l2]
            + 270.0 * (f[..., l1 - 1 : l2 - 1] + f[..., l1 + 1 : l2 + 1])
            - 27.0 * (f[..., l1 - 2 : l2 - 2] + f[..., l1 + 2 : l2 + 2])
            + 2.0 * (f[..., l1 - 3 : l2 - 3] + f[..., l1 + 3 : l2 + 3])
        )
        dfdx[..., :l1] = dfdx[..., l2 - 3 : l2]
        dfdx[..., l2:] = dfdx[..., l1 : l1 + 3]
    else:
        dfdx = 0.0

    return dfdx


def yder2_6th(f, dy):
    """
    yder2_6th(f, dy)

    Compute the 2nd order derivative, 6th order accurate in y.

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

    dy2 = 1.0 / (180.0 * dy ** 2.0)
    dfdy = np.zeros_like(f)
    m1 = 3
    m2 = f.shape[-2] - 3

    if m2 > m1:
        dfdy[..., m1:m2, :] = dy2 * (
            -490.0 * f[..., m1:m2, :]
            + 270.0 * (f[..., m1 - 1 : m2 - 1, :] + f[..., m1 + 1 : m2 + 1, :])
            - 27.0 * (f[..., m1 - 2 : m2 - 2, :] + f[..., m1 + 2 : m2 + 2, :])
            + 2.0 * (f[..., m1 - 3 : m2 - 3, :] + f[..., m1 + 3 : m2 + 3, :])
        )
        dfdy[..., :m1, :] = dfdy[..., m2 - 3 : m2, :]
        dfdy[..., m2:, :] = dfdy[..., m1 : m1 + 3, :]
    else:
        dfdy = 0.0

    return dfdy


def zder2_6th(f, dz):
    """
    zder2_6th(f, dz)

    Compute the 2nd order derivative, 6th order accurate in z.

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

    dz2 = 1.0 / (180.0 * dz ** 2.0)
    dfdz = np.zeros_like(f)
    n1 = 3
    n2 = f.shape[-3] - 3
    if n2 > n1:
        dfdz[..., n1:n2, :, :] = dz2 * (
            -490.0 * f[..., n1:n2, :, :]
            + 270.0 * (f[..., n1 - 1 : n2 - 1, :, :] + f[..., n1 + 1 : n2 + 1, :, :])
            - 27.0 * (f[..., n1 - 2 : n2 - 2, :, :] + f[..., n1 + 2 : n2 + 2, :, :])
            + 2.0 * (f[..., n1 - 3 : n2 - 3, :, :] + f[..., n1 + 3 : n2 + 3, :, :])
        )
        dfdz[..., :n1, :, :] = dfdz[..., n2 - 3 : n2, :, :]
        dfdz[..., n2:, :, :] = dfdz[..., n1 : n1 + 3, :, :]
    else:
        dfdz = 0.0

    return dfdz


def xder5_6th(f, dx):
    """
    xder5_6th(f, dx)

    Compute the 5th order derivative, 6th order accurate in x.

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

    fac = 1 / dx ** 5
    d5fdx = np.zeros_like(f)
    l1 = 3
    l2 = f.shape[-1] - 3

    if l2 > l1:
        d5fdx[..., l1:l2] = fac * (
            +2.5 * (f[..., l1 + 1 : l2 + 1] + f[..., l1 - 1 : l2 - 1])
            - 2.0 * (f[..., l1 + 2 : l2 + 2] + f[..., l1 - 2 : l2 - 2])
            + 0.5 * (f[..., l1 + 3 : l2 + 3] + f[..., l1 - 3 : l2 - 3])
        )
        d5fdx[..., :l1] = d5fdx[..., l2 - 3 : l2]
        d5fdx[..., l2:] = d5fdx[..., l1 : l1 + 3]

    return d5fdx


def yder5_6th(f, dy):
    """
    yder5_6th(f, dy)

    Compute the 5th order derivative, 6th order accurate in y.

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

    fac = 1 / dy ** 5
    m1 = 3
    m2 = f.shape[-2] - 3
    d5fdy = np.zeros_like(f)

    if m2 > m1:
        d5fdy[..., m1:m2, :] = fac * (
            +2.5 * (f[..., m1 + 1 : m2 + 1, :] + f[..., m1 - 1 : m2 - 1, :])
            - 2.0 * (f[..., m1 + 2 : m2 + 2, :] + f[..., m1 - 2 : m2 - 2, :])
            + 0.5 * (f[..., m1 + 3 : m2 + 3, :] + f[..., m1 - 3 : m2 - 3, :])
        )
        d5fdy[..., :m1, :] = d5fdy[..., m2 - 3 : m2, :]
        d5fdy[..., m2:, :] = d5fdy[..., m1 : m1 + 3, :]

    return d5fdy


def zder5_6th(f, dz):
    """
    zder5_6th(f, dy)

    Compute the 5th order derivative, 6th order accurate in z.

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

    fac = 1 / dz ** 5
    n1 = 3
    n2 = f.shape[-3] - 3
    d5fdz = np.zeros_like(f)

    if n2 > n1:
        d5fdz[..., n1:n2, :, :] = fac * (
            +2.5 * (f[..., n1 + 1 : n2 + 1, :, :] + f[..., n1 - 1 : n2 - 1, :, :])
            - 2.0 * (f[..., n1 + 2 : n2 + 2, :, :] + f[..., n1 - 2 : n2 - 2, :, :])
            + 0.5 * (f[..., n1 + 3 : n2 + 3, :, :] + f[..., n1 - 3 : n2 - 3, :, :])
        )

        d5fdz[..., :n1, :, :] = d5fdz[..., n2 - 3 : n2, :, :]
        d5fdz[..., n2:, :, :] = d5fdz[..., n1 : n1 + 3, :, :]

    return d5fdz


def xder6_6th(f, dx):
    """
    xder6_6th(f, dx)

    Compute the 6th order derivative, 6th order accurate in x.

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

    fac = 1 / dx ** 6
    d6fdx = np.zeros_like(f)
    l1 = 3
    l2 = f.shape[-1] - 3

    if l2 > l1:
        d6fdx[..., l1:l2] = fac * (
            -20.0 * f[..., l1:l2]
            + 15.0 * (f[..., l1 + 1 : l2 + 1] + f[..., l1 - 1 : l2 - 1])
            - 6.0 * (f[..., l1 + 2 : l2 + 2] + f[..., l1 - 2 : l2 - 2])
            + (f[..., l1 + 3 : l2 + 3] + f[..., l1 - 3 : l2 - 3])
        )
        d6fdx[..., :l1] = d6fdx[..., l2 - 3 : l2]
        d6fdx[..., l2:] = d6fdx[..., l1 : l1 + 3]

    return d6fdx


def yder6_6th(f, dy):
    """
    yder6_6th(f, dy)

    Compute the 6th order derivative, 6th order accurate in y.

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

    fac = 1 / dy ** 6
    m1 = 3
    m2 = f.shape[-2] - 3
    d6fdy = np.zeros_like(f)

    if m2 > m1:
        d6fdy[..., m1:m2, :] = fac * (
            -20.0 * f[..., m1:m2, :]
            + 15.0 * (f[..., m1 + 1 : m2 + 1, :] + f[..., m1 - 1 : m2 - 1, :])
            - 6.0 * (f[..., m1 + 2 : m2 + 2, :] + f[..., m1 - 2 : m2 - 2, :])
            + (f[..., m1 + 3 : m2 + 3, :] + f[..., m1 - 3 : m2 - 3, :])
        )

        d6fdy[..., :m1, :] = d6fdy[..., m2 - 3 : m2, :]
        d6fdy[..., m2:, :] = d6fdy[..., m1 : m1 + 3, :]
    return d6fdy


def zder6_6th(f, dz):
    """
    zder6_6th(f, dz)

    Compute the 6th order derivative, 6th order accurate in z.

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

    fac = 1 / dz ** 6
    n1 = 3
    n2 = f.shape[-3] - 3
    d6fdz = np.zeros_like(f)

    if n2 > n1:
        d6fdz[..., n1:n2, :, :] = fac * (
            -20.0 * f[..., n1:n2, :, :]
            + 15.0 * (f[..., n1 + 1 : n2 + 1, :, :] + f[..., n1 - 1 : n2 - 1, :, :])
            - 6.0 * (f[..., n1 + 2 : n2 + 2, :, :] + f[..., n1 - 2 : n2 - 2, :, :])
            + (f[..., n1 + 3 : n2 + 3, :, :] + f[..., n1 - 3 : n2 - 3, :, :])
        )
        d6fdz[..., :n1, :, :] = d6fdz[..., n2 - 3 : n2, :, :]
        d6fdz[..., n2:, :, :] = d6fdz[..., n1 : n1 + 3, :, :]

    return d6fdz
