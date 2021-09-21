# integration.py
#
# Contains methods for integrating quantities in different coordinate systems.
#
# Authors: S. Candelaresi (iomsn1@gmail.com).
"""
Contains methods for integrating quantities in different coordinate systems.
"""


def integrate(
    quantity,
    dx=1.0,
    dy=1.0,
    dz=1.0,
    x=None,
    y=None,
    z=None,
    coordinate_system="cartesian",
    axis=(0, 1, 2),
):
    """
    Integrate a field along axis 'axis' using the trapezoidal rule.
    Works for Cartesian, cylindrical and spherical coordinates.
    Works with non-uniform grids.

    call signature:

    integrate(quantity, dx=1.0, dy=1.0, dz=1.0, x=None, y=None, z=None,
              coordinate_system='cartesian', axis=[0, 1, 2]):

    Keyword arguments:

    *quantity*:
      Quantity to be integrated over of shape [nz, ny, nx].

    *dx, dy, dz*:
      Grid spacing in the three dimensions. Not needed when x, y, z
      are specified.

    *x, y, z*:
      Radial (x), polar (y) and vertical (z) coordinates, 1d arrays.
      These override dx, dy and dz.
      Can be non-uniform grids.

    *coordinate_system*:
      Coordinate system under which to take the divergence.
      Takes 'cartesian', 'cylindrical' and 'spherical'.

    *axis*:
      Axis along which to integrate.
    """

    import numpy as np
    import scipy as sc
    import scipy.integrate

    # Perform basic checks on the inputs.
    if quantity.ndim != 3:
        print("ERROR: quantity must have shape [nz, ny, nx].")
        raise ValueError
    if not (x is None):
        if x.shape[0] != quantity.shape[2]:
            print("ERROR: shape of vector x does not match last index of quantity.")
            raise ValueError
    if not (y is None):
        if y.shape[0] != quantity.shape[1]:
            print("ERROR: shape of vector y does not match second index of quantity.")
            raise ValueError
    if len(axis) > 3:
        print("ERROR: axis {0} is not valid.".format(axis))
        raise ValueError
    axis = np.array(axis)
    if not set(axis).issubset({0, 1, 2}):
        print("ERROR: axis {0} is not valid.".format(axis))
        raise ValueError

    # Check the grid.
    if (coordinate_system == "cylindrical") and (x is None):
        print("ERROR: need to specify x (radius) for cylindrical coordinates.")
        raise ValueError
    if (coordinate_system == "spherical") and ((x is None) or (y is None)):
        print(
            "ERROR: need to specify x (radius) and y (polar angle) for spherical coordinates."
        )
        raise ValueError

    # Prepare the grid in case None is specified.
    # The origin is not relevant for the integration for trivial metrics.
    if x is None:
        x = np.linspace(0, (quantity.shape[2] - 1) * dx, quantity.shape[2])
    if y is None:
        y = np.linspace(0, (quantity.shape[1] - 1) * dy, quantity.shape[1])
    if z is None:
        z = np.linspace(0, (quantity.shape[0] - 1) * dz, quantity.shape[0])

    # Multiply the Jacobian to our quantity.
    integral = quantity.copy()
    if (coordinate_system == "cylindrical") and (any(axis == 2)):
        integral *= x[np.newaxis, np.newaxis, :]
    if (coordinate_system == "spherical") and (any(axis == 2)):
        integral *= x[np.newaxis, np.newaxis, :] ** 2
    if (coordinate_system == "spherical") and (any(axis == 1)):
        integral *= np.sin(y[np.newaxis, :, np.newaxis])

    # Perform the integration along the spcified axis.
    # Attention: The order of these integrals is important!!!
    if any(axis == 2):
        integral = sc.integrate.trapz(integral, x=x, axis=2)
    if any(axis == 1):
        integral = sc.integrate.trapz(integral, x=y, axis=1)
    if any(axis == 0):
        integral = sc.integrate.trapz(integral, x=z, axis=0)

    return integral
