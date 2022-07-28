# NOTE: der3, der4, der5 not implemented for nonequidistant grid even in deriv.f90, so we do not bother here.

import warnings
import numpy as np


def der_6th(f, dx_1, axis):
    """
    der_6th(f, dx_1, axis)

    Compute the 1st order derivative, 6th order accurate in x. Adapted from der2_main in deriv.f90.
    This supports nonequidistant grids.

    Parameters
    ----------
    f : ndarray
        Array for which to compute the derivative.

    dx_1 : 1D array
        grid.dx_1 (or dy_1 or dz_1), where grid is a Pencil grid object, and the array passed should correspond to the variable that is differentiated. In the case of equidistant grids, all elements are just 1/grid_spacing.

    axis: int
        Axis of f along which the derivative should be taken.
    """

    nghost = 3

    nax = f.ndim
    if axis >= nax:
        try:
            raise np.AxisError(
                axis, ndim=nax
            )  # May not be available in old numpy versions
        except AttributeError:
            raise ValueError(
                msg_prefix=f"Invalid axis={axis}. Given array has only {nax} axes."
            )

    if np.size(dx_1) != np.size(f, axis):
        raise ValueError(
            "Size of dx_1 ({}) != size of axis of f to be differentiated ({})".format(
                np.size(dx_1), np.size(f, axis)
            )
        )

    f_t = np.moveaxis(f, axis, -1)

    l1 = nghost
    l2 = np.size(f_t,-1) - nghost

    fac = dx_1 / 60

    df_t = np.zeros_like(f_t)

    if l2 > l1:
        df_t[..., l1:l2] = fac[l1:l2] * (
            45.0 * (f_t[..., l1 + 1 : l2 + 1] - f_t[..., l1 - 1 : l2 - 1])
            - 9.0 * (f_t[..., l1 + 2 : l2 + 2] - f_t[..., l1 - 2 : l2 - 2])
            + (f_t[..., l1 + 3 : l2 + 3] - f_t[..., l1 - 3 : l2 - 3])
        )
        df_t[..., :l1] = df_t[..., l2 - nghost : l2]
        df_t[..., l2:] = df_t[..., l1 : l1 + nghost]
    else:
        df_t = 0.0

    return np.moveaxis(df_t, -1, axis)


def der2_6th(f, dx_1, dx_tilde, axis):
    """
    der2_6th(f, dx_1, dx_tilde, axis)

    Compute the 2nd order derivative, 6th order accurate in x. Adapted from der2_main in deriv.f90.
    This supports nonequidistant grids.

    Parameters
    ----------
    f : ndarray
        Array for which to compute the derivative.

    dx_1 : 1D array
        grid.dx_1 (or dy_1 or dz_1), where grid is a Pencil grid object, and the array passed should correspond to the variable that is differentiated. In the case of equidistant grids, all elements are just 1/grid_spacing.

    dx_tilde : 1D array
        grid.dx_tilde (or dy_tilde or dz_tilde). Is just zero for equidistant grids.

    axis: int
        Axis of f along which the derivative should be taken.
    """

    nghost = 3

    nax = f.ndim
    if axis >= nax:
        try:
            raise np.AxisError(
                axis, ndim=nax
            )  # May not be available in old numpy versions
        except AttributeError:
            raise ValueError(
                msg_prefix=f"Invalid axis={axis}. Given array has only {nax} axes."
            )

    if np.size(dx_1) != np.size(f, axis):
        raise ValueError(
            "Size of dx_1 ({}) != size of axis of f to be differentiated ({})".format(
                np.size(dx_1), np.size(f, axis)
            )
        )

    f_t = np.moveaxis(f, axis, -1)

    l1 = nghost
    l2 = np.size(f_t,-1) - nghost

    fac = dx_1**2
    der2_coef0 = -490 / 180
    der2_coef1 = 270 / 180
    der2_coef2 = -27 / 180
    der2_coef3 = 2 / 180

    d2f_t = np.zeros_like(f_t)

    if l2 > l1:
        d2f_t[..., l1:l2] = fac[l1:l2] * (
            der2_coef0 * f_t[..., l1:l2]
            + der2_coef1 * (f_t[..., l1 + 1 : l2 + 1] + f_t[..., l1 - 1 : l2 - 1])
            + der2_coef2 * (f_t[..., l1 + 2 : l2 + 2] + f_t[..., l1 - 2 : l2 - 2])
            + der2_coef3 * (f_t[..., l1 + 3 : l2 + 3] + f_t[..., l1 - 3 : l2 - 3])
        )

        d2f_t[..., l1:l2] += dx_tilde[l1:l2] * der_6th(f_t, dx_1, axis=axis)[..., l1:l2]

        d2f_t[..., :l1] = d2f_t[..., l2 - nghost : l2]
        d2f_t[..., l2:] = d2f_t[..., l1 : l1 + nghost]
    else:
        d2f_t = 0.0

    return np.moveaxis(d2f_t, -1, axis)


def der6_6th(f, dx_1, axis):
    """
    der6_6th(f, dx_1, axis)

    Compute the 6th order derivative, 6th order accurate in x. Adapted from der6_main in deriv.f90.
    This supports nonequidistant grids.

    Parameters
    ----------
    f : ndarray
        Array for which to compute the derivative.

    dx_1 : 1D array
        grid.dx_1 (or dy_1 or dz_1), where grid is a Pencil grid object, and the array passed should correspond to the variable that is differentiated. In the case of equidistant grids, all elements are just 1/grid_spacing.

    axis: int
        Axis of f along which the derivative should be taken.
    """

    nghost = 3

    nax = f.ndim
    if axis >= nax:
        try:
            raise np.AxisError(
                axis, ndim=nax
            )  # May not be available in old numpy versions
        except AttributeError:
            raise ValueError(
                msg_prefix=f"Invalid axis={axis}. Given array has only {nax} axes."
            )

    if np.size(dx_1) != np.size(f, axis):
        raise ValueError(
            "Size of dx_1 ({}) != size of axis of f to be differentiated ({})".format(
                np.size(dx_1), np.size(f, axis)
            )
        )

    f_t = np.moveaxis(f, axis, -1)

    l1 = nghost
    l2 = np.size(f_t,-1) - nghost

    fac = dx_1**6

    d6f_t = np.zeros_like(f_t)

    if l2 > l1:
        d6f_t[..., l1:l2] = fac[l1:l2] * (
            -20 * f_t[..., l1:l2]
            + 15 * (f_t[..., l1 + 1 : l2 + 1] + f_t[..., l1 - 1 : l2 - 1])
            - 6 * (f_t[..., l1 + 2 : l2 + 2] + f_t[..., l1 - 2 : l2 - 2])
            + (f_t[..., l1 + 3 : l2 + 3] + f_t[..., l1 - 3 : l2 - 3])
        )

        d6f_t[..., :l1] = d6f_t[..., l2 - nghost : l2]
        d6f_t[..., l2:] = d6f_t[..., l1 : l1 + nghost]
    else:
        d6f_t = 0.0

    return np.moveaxis(d6f_t, -1, axis)


def xder_6th(f, dx=None, dx_1=None):
    """
    xder_6th(f, dx=None, dx_1=None)

    Compute the 1st order derivative, 6th order accurate in x.

    Parameters
    ----------
    f : ndarray
        Array for which to compute the derivative.

    dx : float
        Grid-spacing in x. For nonequidistant grids, leave this as None and specify dx_1 (=grid.dx_1) instead. If this is specified, the dx_1 argument is ignored.

    dx_1 : ndarray
        Inverse grid spacing. Specify this and leave dx=None if your grid is nonequidistant
    """

    if f.ndim != 3 and f.ndim != 4:
        print("{0} dimension arrays not handled.".format(str(f.ndim)))
        raise ValueError

    if dx is not None:
        warnings.warn("Assuming equidistant grid")
        dx_1 = np.ones(np.size(f,-1)) / dx

    return der_6th(f, dx_1, axis=-1)


def yder_6th(f, dy=None, dy_1=None):
    """
    Same as xder_6th, but for y axis
    """

    if f.ndim != 3 and f.ndim != 4:
        print("{0} dimension arrays not handled.".format(str(f.ndim)))
        raise ValueError

    if dy is not None:
        warnings.warn("Assuming equidistant grid")
        dy_1 = np.ones(np.size(f,-2)) / dy

    return der_6th(f, dy_1, axis=-2)


def zder_6th(f, dz=None, dz_1=None):
    """
    Same as xder_6th, but for z axis
    """

    if f.ndim != 3 and f.ndim != 4:
        print("{0} dimension arrays not handled.".format(str(f.ndim)))
        raise ValueError

    if dz is not None:
        warnings.warn("Assuming equidistant grid")
        dz_1 = np.ones(np.size(f,-3)) / dz

    return der_6th(f, dz_1, axis=-3)

def xder2_6th(f, dx=None, dx_1=None, dx_tilde=None):
    """
    xder2_6th(f, dx=None, dx_1=None)

    Compute the 2nd order derivative, 6th order accurate in x.

    Parameters
    ----------
    f : ndarray
        Array for which to compute the derivative.

    dx : float
        Grid-spacing in x. For nonequidistant grids, leave this as None and specify dx_1 (=grid.dx_1) instead. If this is specified, the dx_1 and dx_tilde arguments are ignored.

    dx_1 : ndarray
        Inverse grid spacing. Specify this and leave dx=None if your grid is nonequidistant
    
    dx_tilde : 1D array
        grid.dx_tilde (or dy_tilde or dz_tilde). Is just zero for equidistant grids.
    """

    if f.ndim != 3 and f.ndim != 4:
        print("{0} dimension arrays not handled.".format(str(f.ndim)))
        raise ValueError

    if dx is not None:
        warnings.warn("Assuming equidistant grid")
        dx_1 = np.ones(np.size(f,-1)) / dx
        dx_tilde = np.zeros(np.size(f,-1))

    return der2_6th(f, dx_1, dx_tilde, axis=-1)


def yder2_6th(f, dy=None, dy_1=None, dy_tilde=None):
    """
    Same as xder2_6th, but for y axis
    """

    if f.ndim != 3 and f.ndim != 4:
        print("{0} dimension arrays not handled.".format(str(f.ndim)))
        raise ValueError

    if dy is not None:
        warnings.warn("Assuming equidistant grid")
        dy_1 = np.ones(np.size(f,-2)) / dy
        dy_tilde = np.zeros(np.size(f,-2))

    return der2_6th(f, dy_1, dy_tilde, axis=-2)


def zder2_6th(f, dz=None, dz_1=None, dz_tilde=None):
    """
    Same as xder2_6th, but for z axis
    """

    if f.ndim != 3 and f.ndim != 4:
        print("{0} dimension arrays not handled.".format(str(f.ndim)))
        raise ValueError

    if dz is not None:
        warnings.warn("Assuming equidistant grid")
        dz_1 = np.ones(np.size(f,-3)) / dz
        dz_tilde = np.zeros(np.size(f,-3))

    return der2_6th(f, dz_1, dz_tilde, axis=-3)


def xder6_6th(f, dx=None, dx_1=None):
    """
    xder6_6th(f, dx=None, dx_1=None)

    Compute the 6th order derivative, 6th order accurate in x.

    Parameters
    ----------
    f : ndarray
        Array for which to compute the derivative.

    dx : float
        Grid-spacing in x. For nonequidistant grids, leave this as None and specify dx_1 (=grid.dx_1) instead. If this is specified, the dx_1 argument is ignored.

    dx_1 : ndarray
        Inverse grid spacing. Specify this and leave dx=None if your grid is nonequidistant
    """

    if f.ndim != 3 and f.ndim != 4:
        print("{0} dimension arrays not handled.".format(str(f.ndim)))
        raise ValueError

    if dx is not None:
        warnings.warn("Assuming equidistant grid")
        dx_1 = np.ones(np.size(f,-1)) / dx

    return der6_6th(f, dx_1, axis=-1)


def yder6_6th(f, dy=None, dy_1=None):
    """
    Same as xder6_6th, but for y axis
    """

    if f.ndim != 3 and f.ndim != 4:
        print("{0} dimension arrays not handled.".format(str(f.ndim)))
        raise ValueError

    if dy is not None:
        warnings.warn("Assuming equidistant grid")
        dy_1 = np.ones(np.size(f,-2)) / dy

    return der6_6th(f, dy_1, axis=-2)


def zder6_6th(f, dz=None, dz_1=None):
    """
    Same as xder6_6th, but for z axis
    """

    if f.ndim != 3 and f.ndim != 4:
        print("{0} dimension arrays not handled.".format(str(f.ndim)))
        raise ValueError

    if dz is not None:
        warnings.warn("Assuming equidistant grid")
        dz_1 = np.ones(np.size(f,-3)) / dz

    return der6_6th(f, dz_1, axis=-3)
