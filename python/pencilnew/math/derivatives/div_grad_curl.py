# div_grad_curl.py
#
# Contains the vector calculus derivatives.
"""
Compute the divergence, gradient and curl.
"""


def div(f, dx, dy, dz, x=None, y=None, coordinate_system='cartesian'):
    """
    Take divervenge of pencil code vector array f in various coordinate systems.

    Keyword arguments:

    *f*:
      Pencil code vector array f.

    *dx, dy, dz*:
      Grid spacing in the three dimensions.

    *x, y*:
      Radial (x) and polar (y) coordinates, 1d arrays.

    *coordinate_system*:
      Coordinate system under which to take the divergence.
      Takes 'cartesian', 'cylindrical' and 'spherical'.
    """

    import numpy as np
    from .der import xder, yder, zder

    if f.ndim != 4:
        print("div: must have vector 4-D array f[mvar, mz, my, mx] for divergence.")
        raise ValueError

    if coordinate_system == 'cartesian':
        return xder(f[0], dx) + yder(f[1], dy) + zder(f[2], dz)

    if coordinate_system == 'cylindrical':
        if x is None:
            print('ERROR: need to specify x (radius) for cylindrical coordinates.')
            raise ValueError
        # Make sure x has compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        return xder(x*f[0], dx)/x + yder(f[1], dy)/x + zder(f[2], dz)

    if coordinate_system == 'spherical':
        if (x is None) or (y is None):
            print('ERROR: need to specify x (radius) and y (polar angle) for spherical coordinates.')
            raise ValueError
        # Make sure x and y have compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        y = y[np.newaxis, :, np.newaxis]
        return xder(x**2*f[0], dx)/x**2 + yder(np.sin(y)*f[1], dy)/(x*np.sin(y)) + zder(f[2], dz)/(x*np.sin(y))

    print('ERROR: could not recognize coordinate system {0}'.format(coordinate_system))
    raise ValueError


def grad(f, dx, dy, dz, x=None, y=None, coordinate_system='cartesian'):
    """
    Take the curl of a pencil code scalar array f in various coordinate systems.

    Keyword arguments:

    *f*:
      Pencil code scalar array f.

    *dx, dy, dz*:
      Grid spacing in the three dimensions.

    *x, y*:
      Radial (x) and polar (y) coordinates, 1d arrays.

    *coordinate_system*:
      Coordinate system under which to take the divergence.
      Takes 'cartesian', 'cylindrical' and 'spherical'.
    """

    import numpy as np
    from .der import xder, yder, zder

    if f.ndim != 3:
        print("grad: must have scalar 3-D array f[mz, my, mx] for gradient.")
        raise ValueError

    grad_value = np.empty((3,) + f.shape)

    if coordinate_system == 'cartesian':
        grad_value[0] = xder(f, dx)
        grad_value[1] = yder(f, dy)
        grad_value[2] = zder(f, dz)

    if coordinate_system == 'cylindrical':
        if x is None:
            print('ERROR: need to specify x (radius) for cylindrical coordinates.')
            raise ValueError
        # Make sure x has compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        grad_value[0] = xder(f, dx)
        grad_value[1] = yder(f, dy)/x
        grad_value[2] = zder(f, dz)

    if coordinate_system == 'spherical':
        if (x is None) or (y is None):
            print('ERROR: need to specify x (radius) and y (polar angle) for spherical coordinates.')
            raise ValueError
        # Make sure x and y have compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        y = y[np.newaxis, :, np.newaxis]
        grad_value[0] = xder(f, dx)
        grad_value[1] = yder(f, dy)/x
        grad_value[2] = zder(f, dz)/(x*np.sin(y))

    return grad_value


def curl(f, dx, dy, dz, x=None, y=None, run2D=False, coordinate_system='cartesian'):
    """
    Take the curl of a pencil code vector array f in various coordinate systems.

    Keyword arguments:

    *f*:
      Pencil code scalar array f.

    *dx, dy, dz*:
      Grid spacing in the three dimensions.

    *x, y*:
      Radial (x) and polar (y) coordinates, 1d arrays.

    *run2D*:
      Deals with pure 2-D snapshots (solved the (x,z)-plane pb).
      !Only for Cartesian grids at the moment!

    *coordinate_system*:
      Coordinate system under which to take the divergence.
      Takes 'cartesian', 'cylindrical' and 'spherical'.
      !Does not work for 2d runs yet!
    """

    import numpy as np
    from .der import xder, yder, zder

    if f.shape[0] != 3:
        print("curl: must have vector 4-D array f[3, mz, my, mx] for curl.")
        raise ValueError

    curl_value = np.empty_like(f)
    if (dy != 0. and dz != 0.):
        # 3-D case
        if coordinate_system == 'cartesian':
            curl_value[0] = yder(f[2], dy) - zder(f[1], dz)
            curl_value[1] = zder(f[0], dz) - xder(f[2], dx)
            curl_value[2] = xder(f[1], dx) - yder(f[0], dy)
        if coordinate_system == 'cylindrical':
            if x is None:
                print('ERROR: need to specify x (radius) for cylindrical coordinates.')
                raise ValueError
            # Make sure x has compatible dimensions.
            x = x[np.newaxis, np.newaxis, :]
            curl_value[0] = (1/x)*yder(f[2], dy) - zder(f[1], dz)
            curl_value[1] = zder(f[0], dz) - xder(f[2], dx)
            curl_value[2] = (1/x)*xder(x*f[1], dx) - (1/x)*yder(f[0], dy)
        if coordinate_system == 'spherical':
            if (x is None) or (y is None):
                print('ERROR: need to specify x (radius) and y (polar angle) for spherical coordinates.')
                raise ValueError
            # Make sure x and y have compatible dimensions.
            x = x[np.newaxis, np.newaxis, :]
            y = y[np.newaxis, :, np.newaxis]
            curl_value[0] = (yder(np.sin(y)*f[2], dy) - zder(f[1], dz))/(x*np.sin(y))
            curl_value[1] = (zder(f[0], dz)/np.sin(y) - xder(x*f[2], dx))/x
            curl_value[2] = (xder(x*f[1], dx) - yder(f[0], dy))/x
    elif dy == 0.:
        # 2-D case in the xz-plane
        curl_value[0] = zder(f, dz, run2D)[0] - xder(f, dx)[2]
    else:
        # 2-D case in the xy-plane
        curl_value[0] = xder(f, dx)[1] - yder(f, dy)[0]

    return curl_value


def curl2(f, dx, dy, dz, x=None, y=None, coordinate_system='cartesian'):
    """
    Take the double curl of a pencil code vector array f.

    *x, y*:
      Radial (x) and polar (y) coordinates, 1d arrays.

    *run2D*:
      Deals with pure 2-D snapshots (solved the (x,z)-plane pb).
      !Only for Cartesian grids at the moment!

    *coordinate_system*:
      Coordinate system under which to take the divergence.
      Takes 'cartesian' and 'cylindrical'.
    """

    import numpy as np
    from .der import xder, yder, zder, xder2, yder2, zder2

    if (f.ndim != 4 or f.shape[0] != 3):
        print("curl2: must have vector 4-D array f[3, mz, my, mx] for curl2.")
        raise ValueError

    if coordinate_system == 'spherical':
        print('ERROR: curl2 currently not implemented for spherical coordinates.')

    curl2_value = np.empty(f.shape)

    if coordinate_system == 'cartesian':
        curl2_value[0] = xder(yder(f[1], dy) + zder(f[2], dz), dx) \
                              -yder2(f[0], dy) - zder2(f[0], dz)
        curl2_value[1] = yder(xder(f[0], dx) + zder(f[2], dz), dy) \
                              -xder2(f[1], dx) - zder2(f[1], dz)
        curl2_value[2] = zder(xder(f[0], dx) + yder(f[1], dy), dz) \
                              -xder2(f[2], dx) - yder2(f[2], dy)
    if coordinate_system == 'cylindrical':
        if x is None:
            print('ERROR: need to specify x (radius) for cylindrical coordinates.')
            raise ValueError
        # Make sure x has compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        curl2_value[0] = yder(f[1], dy)/x**2 + xder(yder(f[1], dy), dx)/x \
                         - yder2(f[0], dy)/x**2 - zder2(f[0], dz) + xder(zder(f[2], dz), dx)
        curl2_value[1] = yder(zder(f[2], dz), dy)/x - zder2(f[1], dz) + f[1]/x**2 \
                         - xder(f[1], dx)/x - xder2(f[1], dx) + xder(yder(f[0], dy), dx)/x \
                         - yder(f[0], dy)/x**2
        curl2_value[2] = zder(f[0], dz)/x + xder(zder(f[0], dz), dx) - zder(f[2], dx)/x \
                         - xder2(f[2], dx) - yder2(f[2], dy)/x**2 + yder(zder(f[1], dz), dy)/x

    return curl2_value


def del2(f, dx, dy, dz):
    """
    Calculate del2.
    """

    from .der import xder2, yder2, zder2

    return xder2(f, dx) + yder2(f, dy) + zder2(f, dz)


def del6(f, dx, dy, dz):
    """
    Calculate del6 (defined here as d^6/dx^6 + d^6/dy^6 + d^6/dz^6, rather
    than del2^3) of a scalar f for hyperdiffusion.
    """

    from .der import xder6, yder6, zder6

    return xder6(f, dx) + yder6(f, dy) + zder6(f, dz)
