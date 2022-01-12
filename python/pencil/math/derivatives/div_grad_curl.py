# div_grad_curl.py
#
# Contains the vector calculus derivatives.
#
# Authors: Wladimir Lyra (wlyra@caltech.edu)
#          Simon Candelaresi (iomsn1@gmail.com)
#          Callum Reid (apollocreid@gmail.com)
#
# TODO: Include non-equidistant coordinates.
"""
Compute the divergence, gradient and curl.
"""


def div(f, dx, dy, dz, x=None, y=None, coordinate_system="cartesian"):
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
    from pencil.math.derivatives.der import xder, yder, zder

    if f.ndim != 4:
        print("div: must have vector 4-D array f[mvar, mz, my, mx] for divergence.")
        raise ValueError

    if coordinate_system == "cartesian":
        return xder(f[0], dx) + yder(f[1], dy) + zder(f[2], dz)

    if coordinate_system == "cylindrical":
        if x is None:
            print("ERROR: need to specify x (radius) for cylindrical coordinates.")
            raise ValueError
        # Make sure x has compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        return xder(x * f[0], dx) / x + yder(f[1], dy) / x + zder(f[2], dz)

    if coordinate_system == "spherical":
        if (x is None) or (y is None):
            print(
                "ERROR: need to specify x (radius) and y (polar angle) for spherical coordinates."
            )
            raise ValueError
        # Make sure x and y have compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        y = y[np.newaxis, :, np.newaxis]
        return (
            xder(x ** 2 * f[0], dx) / x ** 2
            + yder(np.sin(y) * f[1], dy) / (x * np.sin(y))
            + zder(f[2], dz) / (x * np.sin(y))
        )

    print("ERROR: could not recognize coordinate system {0}".format(coordinate_system))
    raise ValueError


def grad(f, dx, dy, dz, x=None, y=None, coordinate_system="cartesian"):
    """
    Take the gradient of a pencil code scalar array f in various coordinate systems.

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
    from pencil.math.derivatives.der import xder, yder, zder

    if f.ndim != 3:
        print("grad: must have scalar 3-D array f[mz, my, mx] for gradient.")
        raise ValueError

    grad_value = np.zeros((3,) + f.shape)

    if coordinate_system == "cartesian":
        grad_value[0] = xder(f, dx)
        grad_value[1] = yder(f, dy)
        grad_value[2] = zder(f, dz)

    if coordinate_system == "cylindrical":
        if x is None:
            print("ERROR: need to specify x (radius) for cylindrical coordinates.")
            raise ValueError
        # Make sure x has compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        grad_value[0] = xder(f, dx)
        grad_value[1] = yder(f, dy) / x
        grad_value[2] = zder(f, dz)

    if coordinate_system == "spherical":
        if (x is None) or (y is None):
            print(
                "ERROR: need to specify x (radius) and y (polar angle) for spherical coordinates."
            )
            raise ValueError
        # Make sure x and y have compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        y = y[np.newaxis, :, np.newaxis]
        grad_value[0] = xder(f, dx)
        grad_value[1] = yder(f, dy) / x
        grad_value[2] = zder(f, dz) / (x * np.sin(y))

    return grad_value


def curl(f, dx, dy, dz, x=None, y=None, run2D=False, coordinate_system="cartesian"):
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
    from pencil.math.derivatives.der import xder, yder, zder

    if f.shape[0] != 3:
        print("curl: must have vector 4-D array f[3, mz, my, mx] for curl.")
        raise ValueError

    curl_value = np.zeros_like(f)

    if dy != 0.0 and dz != 0.0:
        # 3-D case
        if coordinate_system == "cartesian":
            curl_value[0] = yder(f[2], dy) - zder(f[1], dz)
            curl_value[1] = zder(f[0], dz) - xder(f[2], dx)
            curl_value[2] = xder(f[1], dx) - yder(f[0], dy)
        if coordinate_system == "cylindrical":
            if x is None:
                print("ERROR: need to specify x (radius) for cylindrical coordinates.")
                raise ValueError
            # Make sure x has compatible dimensions.
            x = x[np.newaxis, np.newaxis, :]
            curl_value[0] = (1 / x) * yder(f[2], dy) - zder(f[1], dz)
            curl_value[1] = zder(f[0], dz) - xder(f[2], dx)
            curl_value[2] = (1 / x) * xder(x * f[1], dx) - (1 / x) * yder(f[0], dy)
        if coordinate_system == "spherical":
            if (x is None) or (y is None):
                print(
                    "ERROR: need to specify x (radius) and y (polar angle) for spherical coordinates."
                )
                raise ValueError
            # Make sure x and y have compatible dimensions.
            x = x[np.newaxis, np.newaxis, :]
            y = y[np.newaxis, :, np.newaxis]
            curl_value[0] = (yder(np.sin(y) * f[2], dy) - zder(f[1], dz)) / (
                x * np.sin(y)
            )
            curl_value[1] = (zder(f[0], dz) / np.sin(y) - xder(x * f[2], dx)) / x
            curl_value[2] = (xder(x * f[1], dx) - yder(f[0], dy)) / x
    elif dy == 0.0:
        # 2-D case in the xz-plane
        curl_value[0] = zder(f, dz, run2D)[0] - xder(f, dx)[2]
    else:
        # 2-D case in the xy-plane
        curl_value[0] = xder(f, dx)[1] - yder(f, dy)[0]

    return curl_value


def curl2(f, dx, dy, dz, x=None, y=None, coordinate_system="cartesian"):
    """
    Take the double curl of a pencil code vector array f.

    *x, y*:
      Radial (x) and polar (y) coordinates, 1d arrays.

    *run2D*:
      Deals with pure 2-D snapshots (solved the (x,z)-plane pb).
      !Only for Cartesian grids at the moment!

    *coordinate_system*:
      Coordinate system under which to take the divergence.
      Takes 'cartesian', 'cylindrical' and 'spherical'.
    """

    import numpy as np
    from pencil.math.derivatives.der import xder, yder, zder, xder2, yder2, zder2

    if f.ndim != 4 or f.shape[0] != 3:
        print("curl2: must have vector 4-D array f[3, mz, my, mx] for curl2.")
        raise ValueError

    curl2_value = np.zeros(f.shape)

    if coordinate_system == "cartesian":
        curl2_value[0] = (
            xder(yder(f[1], dy) + zder(f[2], dz), dx)
            - yder2(f[0], dy)
            - zder2(f[0], dz)
        )
        curl2_value[1] = (
            yder(xder(f[0], dx) + zder(f[2], dz), dy)
            - xder2(f[1], dx)
            - zder2(f[1], dz)
        )
        curl2_value[2] = (
            zder(xder(f[0], dx) + yder(f[1], dy), dz)
            - xder2(f[2], dx)
            - yder2(f[2], dy)
        )
    if coordinate_system == "cylindrical":
        if x is None:
            print("ERROR: need to specify x (radius) for cylindrical coordinates.")
            raise ValueError
        # Make sure x has compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        curl2_value[0] = (
            yder(f[1], dy) / x ** 2
            + xder(yder(f[1], dy), dx) / x
            - yder2(f[0], dy) / x ** 2
            - zder2(f[0], dz)
            + xder(zder(f[2], dz), dx)
        )
        curl2_value[1] = (
            yder(zder(f[2], dz), dy) / x
            - zder2(f[1], dz)
            + f[1] / x ** 2
            - xder(f[1], dx) / x
            - xder2(f[1], dx)
            + xder(yder(f[0], dy), dx) / x
            - yder(f[0], dy) / x ** 2
        )
        curl2_value[2] = (
            zder(f[0], dz) / x
            + xder(zder(f[0], dz), dx)
            - zder(f[2], dx) / x
            - xder2(f[2], dx)
            - yder2(f[2], dy) / x ** 2
            + yder(zder(f[1], dz), dy) / x
        )
    if coordinate_system == "spherical":
        if x is None or y is None:
            print(
                "ERROR: need to specify x (radius) and y (polar angle) for spherical coordinates."
            )
            raise ValueError
        # Make sure x and y have compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        y = y[np.newaxis, :, np.newaxis]
        curl2_value[0] = (
            yder(np.sin(y) * (xder(x * f[1], dx) - yder(f[0], dy)) / x, dy)
            - zder((zder(f[0], dz) / np.sin(y) - xder(x * f[2], dx)) / x, dz)
        ) / (x * np.sin(y))
        curl2_value[1] = (
            zder((yder(np.sin(y) * f[2], dy) - zder(f[1], dz)) / (x * np.sin(y)), dz)
            / np.sin(y)
            - xder((xder(x * f[1], dx) - yder(f[0], dy)), dx)
        ) / x
        curl2_value[2] = (
            xder((zder(f[0], dz) / np.sin(y) - xder(x * f[2], dx)), dx)
            - yder((yder(np.sin(y) * f[2], dy) - zder(f[1], dz)) / (x * np.sin(y)), dy)
        ) / x
    return curl2_value


def del2(f, dx, dy, dz, x=None, y=None, coordinate_system="cartesian"):
    """
    Calculate del2, the Laplacian of a scalar field f.

    *x, y*:
      Radial (x) and polar (y) coordinates, 1d arrays.

    *run2D*:
      Deals with pure 2-D snapshots (solved the (x,z)-plane pb).
      !Only for Cartesian grids at the moment!

    *coordinate_system*:
      Coordinate system under which to take the divergence.
      Takes 'cartesian', 'cylindrical' and 'spherical'.
    """

    import numpy as np
    from pencil.math.derivatives.der import xder2, yder2, zder2, xder, yder

    if coordinate_system == "cartesian":
        del2_value = xder2(f, dx) + yder2(f, dy) + zder2(f, dz)
    if coordinate_system == "cylindrical":
        if x is None:
            print("ERROR: need to specify x (radius)")
            raise ValueError
        # Make sure x has compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        del2_value = (
            xder(f, dx) / x + xder2(f, dx) + yder2(f, dy) / (x ** 2) + zder2(f, dz)
        )
    if coordinate_system == "spherical":
        if x is None or y is None:
            print("ERROR: need to specify x (radius) and y (polar angle)")
            raise ValueError
        # Make sure x and y have compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        y = y[np.newaxis, :, np.newaxis]
        del2_value = (
            2 * xder(f, dx) / x
            + xder2(f, dx)
            + np.cos(y) * yder(f, dy) / ((x ** 2) * np.sin(y))
            + yder2(f, dy) / (x ** 2)
            + zder2(f, dz) / ((x * np.sin(y)) ** 2)
        )
    return del2_value


def del2v(f, dx, dy, dz, x=None, y=None, coordinate_system="cartesian"):
    """
    Calculate del2, the Laplacian of a vector field f.

    *x, y*:
      Radial (x) and polar (y) coordinates, 1d arrays.

    *run2D*:
      Deals with pure 2-D snapshots (solved the (x,z)-plane pb).
      !Only for Cartesian grids at the moment!

    *coordinate_system*:
      Coordinate system under which to take the divergence.
      Takes 'cartesian', 'cylindrical' and 'spherical'.
    """

    if f.shape[0] != 3:
        print(
            "Vector Laplacian: must have a vector 4D array f(3, mz, my, mx) for Vector Laplacian"
        )
        raise ValueError

    import numpy as np
    from pencil.math.derivatives.der import xder2, yder2, zder2, yder, zder

    del2v_value = np.zeros(f.shape)

    if coordinate_system == "cartesian":
        del2v_value[0] = xder2(f[0], dx) + yder2(f[0], dy) + zder2(f[0], dz)
        del2v_value[1] = xder2(f[1], dx) + yder2(f[1], dy) + zder2(f[1], dz)
        del2v_value[2] = xder2(f[2], dx) + yder2(f[2], dy) + zder2(f[2], dz)
    if coordinate_system == "cylindrical":
        if x is None:
            print("Error: need to specify x (radius)")
            raise ValueError
        # Make sure x has compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        del2v_value[0] = (
            del2(f[0], dx, dy, dz, x=x, y=y, coordinate_system="cylindrical")
            - f[0] / (x ** 2)
            - 2 * yder(f[1], dy) / x ** 2
        )
        del2v_value[1] = (
            del2(f[1], dx, dy, dz, x=x, y=y, coordinate_system="cylindrical")
            - f[1] / (x ** 2)
            + 2 * yder(f[0], dy) / x ** 2
        )
        del2v_value[2] = del2(
            f[2], dx, dy, dz, x=x, y=y, coordinate_system="cylindrical"
        )
    if coordinate_system == "spherical":
        if x is None or y is None:
            print("ERROR: need to specify x (radius) and y (polar angle)")
        # Make sure x and y have compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        y = y[np.newaxis, :, np.newaxis]
        del2v_value[0] = (
            del2(f[0], dx, dy, dz, x=x, y=y, coordinate_system="spherical")
            - 2 * f[0] / x ** 2
            - 2 * yder(f[1], dy) / x ** 2
            - 2 * np.cos(y) * f[1] / (x ** 2 * np.sin(y))
            - 2 * zder(f[2], dz) / (x ** 2 * np.sin(y))
        )
        del2v_value[1] = (
            del2(f[1], dx, dy, dz, x=x, y=y, coordinate_system="spherical")
            - f[1] / (x * np.sin(y)) ** 2
            + 2 * yder(f[0], dy) / (x ** 2)
            - (2 * np.cos(y)) * zder(f[2], dz) / (x * np.sin(y)) ** 2
        )
        del2v_value[2] = (
            del2(f[2], dx, dy, dz, x=x, y=y, coordinate_system="spherical")
            - f[2] / (x * np.sin(y)) ** 2
            + 2 * zder(f[0], dz) / (np.sin(y) * x ** 2)
            - (2 * np.cos(y)) * zder(f[1], dz) / (x * np.sin(y)) ** 2
        )
    return del2v_value


def curl3(f, dx, dy, dz, x=None, y=None, coordinate_system="cartesian"):
    """
    Take the triple curl of a pencil code vector array f.

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
    from pencil.math.derivatives.der import (
        xder,
        yder,
        zder,
        xder2,
        yder2,
        zder2,
        xder3,
        yder3,
        zder3,
    )

    if f.ndim != 4 or f.shape[0] != 3:
        print("curl3: must have vector 4-D array f[3, mz, my, mx] for curl3.")
        raise ValueError

    curl3_value = np.zeros(f.shape)

    if coordinate_system == "cartesian":
        curl3_value[0] = (
            zder(xder2(f[1], dx) + yder2(f[1], dy), dz)
            + zder3(f[1], dz)
            - yder(xder2(f[2], dx) + zder2(f[2], dz), dy)
            - yder3(f[2], dy)
        )
        curl3_value[1] = (
            xder3(f[2], dx)
            + xder(yder2(f[2], dy) + zder2(f[2], dz), dx)
            - zder(xder2(f[0], dx) + yder2(f[0], dy), dz)
            - zder3(f[0], dz)
        )
        curl3_value[2] = (
            yder(xder2(f[0], dx) + zder2(f[0], dz), dy)
            + yder3(f[0], dy)
            - xder3(f[1], dx)
            - xder(yder2(f[1], dy) + zder2(f[1], dz), dx)
        )
    if coordinate_system == "cylindrical":
        if x is None:
            print("ERROR: need to specify x (radius) for cylindrical coordinates.")
            raise ValueError
        # Make sure x has compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        curl3_value[0] = (
            2 * yder(zder(f[0], dz), dy) / x ** 2
            - yder(zder(f[2], dz), dy) / x ** 2
            - xder2(yder(f[2], dy), dx) / x
            - yder3(f[2], dy) / x ** 3
            + yder2(zder(f[1], dz), dy) / x ** 2
            - yder(zder2(f[2], dz), dy) / x
            + zder3(f[1], dz)
            - zder(f[1], dz) / x ** 2
            + xder(zder(f[1], dz), dx) / x
            + xder2(zder(f[1], dz), dx)
        )
        curl3_value[1] = (
            2 * yder(zder(f[1], dz), dy) / x ** 2
            - yder2(zder(f[0], dz), dy) / x ** 2
            - zder3(f[0], dz)
            + xder(zder2(f[2], dz), dx)
            - xder(zder(f[0], dz), dx) / x
            + zder(f[0], dz) / x ** 2
            - xder2(zder(f[0], dz), dx)
            + xder(zder(f[2], dz), dx) / x
            - zder(f[2], dz) / x ** 2
            + xder3(f[2], dx)
            + xder(yder2(f[2], dy), dx) / x ** 2
            - 2 * yder2(f[2], dy) / x ** 3
        )
        curl3_value[2] = (
            -xder(zder2(f[1], dz), dx)
            - zder2(f[1], dz) / x
            + xder(f[1], dx) / x ** 2
            - f[1] / x ** 3
            - xder3(f[1], dx)
            + xder2(yder(f[0], dy), dx) / x
            - xder(yder(f[0], dy), dx) / x ** 2
            + yder(f[0], dy) / x ** 3
            - yder2(f[1], dy) / x ** 3
            - xder(yder2(f[1], dy), dx) / x ** 2
            + yder3(f[0], dy) / x ** 3
            + yder(zder2(f[0], dz), dy) / x
            - 2 * xder2(f[1], dx) / x
        )
    if coordinate_system == "spherical":
        if x is None or y is None:
            print(
                "ERROR: need to specify x (radius) and y (polar angle) for spherical coordinates"
            )
            raise ValueError
        # Make sure x and y have compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        y = y[np.newaxis, :, np.newaxis]
        # TODO: This one needs some fixing.
        curl3_value[0] = (-1 / x ** 3) * (
            np.sin(y)
            * (
                (np.cos(y) ** 2 - 1) * (yder2(zder(f[1], dz), dy))
                + (-np.cos(y) ** 2 * np.sin(y) + np.sin(y)) * (yder3(f[2], dy))
                - x ** 2 * (xder2(zder(f[1], dz), dx))
                - zder3(f[1], dz)
                + np.sin(y) * (xder2(yder(f[2], dy), dx)) * x ** 2
                + (yder(zder2(f[2], dz), dy)) * np.sin(y)
                + (-6 * np.cos(y) ** 3 + 6 * np.cos(y)) * (yder2(f[2], dy))
                - 2 * x * (xder(zder(f[1], dz), dx))
                - 3 * np.sin(y) * np.cos(y) * (yder(zder(f[1], dz), dy))
                + np.cos(y) * (xder2(f[2], dx)) * x ** 2
                + 2 * np.sin(y) * (xder(yder(f[2], dy), dx)) * x
                + (zder2(f[2], dz)) * np.cos(y)
                + (-2 * np.cos(y) ** 2 + 1) * (zder(f[1], dz))
                + (11 * np.cos(y) ** 2 * np.sin(y) - 4 * np.sin(y)) * (yder(f[2], dy))
                + 2 * np.cos(y) * (xder(f[2], dx)) * x
                + (6 * np.cos(y) ** 3 - 5 * np.cos(y)) * f[2]
            )
        )
        curl3_value[1] = (
            1
            / (x ** 3 * np.sin(y))
            * (
                (np.cos(y) ** 2 - 1) * (yder2(zder(f[0], dz), dy))
                + (-x * np.cos(y) ** 2 + x) * np.sin(y) * (xder(yder2(f[2], dy), dx))
                - (xder2(zder(f[0], dz), dx)) * x ** 2
                - zder3(f[0], dz)
                + x ** 3 * (xder3(f[2], dx)) * np.sin(y)
                + x * (xder(zder2(f[2], dz), dx)) * np.sin(y)
                + (-2 * np.cos(y) ** 2 + 2) * (yder(zder(f[1], dz), dy))
                + 3 * x ** 2 * xder2(f[2], dx) * np.sin(y)
                + zder2(f[2], dz) * np.sin(y)
                + (2 * x * np.cos(y) ** 2 - x) * np.sin(y) * (xder(f[2], dx))
                + (3 * np.cos(y) ** 3 - 3 * np.cos(y)) * (yder(f[2], dy))
                + 2 * np.sin(y) * np.cos(y) * (zder(f[1], dz))
                + (-2 * np.cos(y) ** 2 + 1) * np.sin(y) * f[2]
            )
        )
        curl3_value[2] = (
            1
            / x ** 3
            * (
                (yder3(f[0], dy)) * (1 - np.cos(y) ** 2)
                + (x * np.cos(y) ** 2 - x) * (xder(yder2(f[1], dy), dx))
                + (xder2(yder(f[1], dy), dx)) * x ** 2
                + yder(zder2(f[0], dz), dy)
                - (xder3(f[1], dx)) * x ** 3
                - (xder(zder2(f[1], dz), dx)) * x
                + (np.cos(y) ** 2 - 1) * (yder2(f[1], dy))
                + 3 * (yder2(f[1], dy)) * np.sin(y) * np.cos(y)
                - 3 * (xder2(f[1], dx)) * x ** 2
                - 3 * (xder(yder(f[1], dy), dx)) * np.sin(y) * np.cos(y) * x
                + zder2(f[1], dz)
                - 2 * (yder(zder(f[2], dz), dy)) * np.sin(y)
                + (2 * np.cos(y) ** 2 - 1) * (yder(f[1], dy))
                + (-2 * np.cos(y) ** 2 + x) * (xder(f[1], dx))
                - 3 * (yder(f[1], dy)) * np.sin(y) * np.cos(y)
                - 2 * (zder(f[2], dz)) * np.cos(y)
                + (-2 * np.cos(y) ** 2 + 1) * f[1]
            )
        )
    return curl3_value


def del6(f, dx, dy, dz):
    """
    Calculate del6 (defined here as d^6/dx^6 + d^6/dy^6 + d^6/dz^6, rather
    than del2^3) of a scalar f for hyperdiffusion.
    """

    from pencil.math.derivatives.der import xder6, yder6, zder6

    return xder6(f, dx) + yder6(f, dy) + zder6(f, dz)


def gij(f, dx, dy, dz, nder=6):
    """
    Calculate del6 (defined here as d^6/dx^6 + d^6/dy^6 + d^6/dz^6, rather
    than del2^3) of a scalar f for hyperdiffusion.
    """
    import numpy as np

    if f.ndim != 4 or f.shape[0] != 3:
        print("curl3: must have vector 4-D array f[3, mz, my, mx] for curl3.")
        raise ValueError
    gij = np.array(f, f, f)
    for i in range(3):
        if nder == 1:
            from pencil.math.derivatives.der import xder, yder, zder

            gij[i, 0] = xder(f[i], dx, dy, dz)
            gij[i, 1] = yder(f[i], dx, dy, dz)
            gij[i, 2] = zder(f[i], dx, dy, dz)
        elif nder == 2:
            from pencil.math.derivatives.der import xder2, yder2, zder2

            gij[i, 0] = xder2(f[i], dx, dy, dz)
            gij[i, 1] = yder2(f[i], dx, dy, dz)
            gij[i, 2] = zder2(f[i], dx, dy, dz)
        elif nder == 3:
            from pencil.math.derivatives.der import xder3, yder3, zder3

            gij[i, 0] = xder3(f[i], dx, dy, dz)
            gij[i, 1] = yder3(f[i], dx, dy, dz)
            gij[i, 2] = zder3(f[i], dx, dy, dz)
        elif nder == 4:
            from pencil.math.derivatives.der import xder4, yder4, zder4

            gij[i, 0] = xder4(f[i], dx, dy, dz)
            gij[i, 1] = yder4(f[i], dx, dy, dz)
            gij[i, 2] = zder4(f[i], dx, dy, dz)
        elif nder == 5:
            from pencil.math.derivatives.der import xder5, yder5, zder5

            gij[i, 0] = xder5(f[i], dx, dy, dz)
            gij[i, 1] = yder5(f[i], dx, dy, dz)
            gij[i, 2] = zder5(f[i], dx, dy, dz)
        elif nder == 6:
            from pencil.math.derivatives.der import xder6, yder6, zder6

            gij[i, 0] = xder6(f[i], dx, dy, dz)
            gij[i, 1] = yder6(f[i], dx, dy, dz)
            gij[i, 2] = zder6(f[i], dx, dy, dz)

    return gij


def traceless_strain(
    f, dx, dy, dz, x=None, y=None, z=None, coordinate_system="cartesian"
):

    import numpy as np

    if f.ndim != 4 or f.shape[0] != 3:
        print("curl3: must have vector 4-D array f[3, mz, my, mx] for curl3.")
        raise ValueError

    uij = gij(f, dx, dy, dz, nder=1)
    divu = div(f, dx, dy, dz, x=x, y=y, coordinate_system=coordinate_system)
    sij = uij.copy()
    for j in range(3):
        sij[j, j] = uij[j, j] - (1.0 / 3.0) * divu
        for i in range(j + 1, 3):
            sij[i, j] = 0.5 * (uij[i, j] + uij[j, i])
            sij[j, i] = sij[i, j]

    if coordinate_system == "cylindrical":
        if x is None:
            print("ERROR: need to specify x (radius) for cylindrical coordinates.")
            raise ValueError
        # Make sure x has compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        sij[1, 2] = sij[1, 2] - 0.5 * f[2] / x
        sij[2, 2] = sij[2, 2] + 0.5 * f[1] / x
        sij[2, 1] = sij[1, 2]

    if coordinate_system == "spherical":
        if (x is None) or (y is None):
            print(
                "ERROR: need to specify x (radius) and y (polar angle) for spherical coordinates."
            )
            raise ValueError
        # Make sure x and y have compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        y = y[np.newaxis, :, np.newaxis]
        # sij[1,1] remains unchanged in spherical coordinates
        sij[1, 2] = sij[1, 2] - 0.5 / x * f[2]
        sij[1, 3] = sij[1, 3] - 0.5 / x * f[3]
        sij[2, 1] = sij[1, 2]
        sij[2, 2] = sij[2, 2] + f[1] / x
        # sij[2,3]=sij[2,3]-.5/x*cotth(m)*f[3]
        sij[3, 1] = sij[1, 3]
        sij[3, 2] = sij[2, 3]
        # sij[3,3]=sij[3,3]+f[1]/x+cotth(m)/x*f[2]
        # if (lshear_ROS) then
        #  sij(:,1,2)=sij(:,1,2)+Sshear
        #  sij(:,2,1)=sij(:,2,1)+Sshear
    return sij


# def sij2(f, dx, dy, dz, x=None, y=None, z=None,
#                     coordinate_system='cartesian'):
#
#    if (f.ndim != 4 or f.shape[0] != 3):
#        print("curl3: must have vector 4-D array f[3, mz, my, mx] for curl3.")
#        raise ValueError
# ):
#    uij = gij(f, dx, dy, dz, nder=1)
#! divu -> uij2
#      call div_mn(uij,sij2,uu)
#    sij = traceless_strain(f, dx, dy, dz, x=x, y=y, z=z,
#                     coordinate_system=coordinate_system)
#      call traceless_strain(uij,sij2,sij,uu,lshear_rateofstrain)
#! sij^2
#      call multm2_sym_mn(sij,sij2)
#
