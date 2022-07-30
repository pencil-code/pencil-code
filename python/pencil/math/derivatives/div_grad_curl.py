# div_grad_curl.py
#
# Contains the vector calculus derivatives.
#
# Authors: Wladimir Lyra (wlyra@caltech.edu)
#          Simon Candelaresi (iomsn1@gmail.com)
#          Callum Reid (apollocreid@gmail.com)
#          Kishore G. (kishore96@gmail.com)
#
"""
Compute the divergence, gradient and curl.
"""

import warnings


def div(
    f,
    dx=None,
    dy=None,
    dz=None,
    x=None,
    y=None,
    coordinate_system="cartesian",
    grid=None,
):
    """
    div(f, dx=None, dy=None, dz=None, x=None, y=None, coordinate_system="cartesian", grid=None)

    Take divergence of pencil code vector array f in various coordinate systems.

    Parameters
    ----------
    f : ndarray
        Pencil code vector array f.

    grid : pencil.read.grids.Grid
        Pencil grid object. See pc.read.grid().

    coordinate_system : string
        Coordinate system under which to take the divergence.
        Takes 'cartesian', 'cylindrical' and 'spherical'.

    Deprecated parameters (only for backwards compatibility)
    --------------------------------------------------------
    dx, dy, dz : floats
        Grid spacing in the three dimensions. These will not have any effect if grid!=None

    x, y : ndarrays
        Radial (x) and polar (y) coordinates, 1d arrays. These will not have any effect if grid!=None

    """

    import numpy as np
    from pencil.math.derivatives.der import xder, yder, zder

    if f.ndim != 4:
        print("div: must have vector 4-D array f[mvar, mz, my, mx] for divergence.")
        raise ValueError

    if grid is not None:
        x = grid.x
        y = grid.y
        dx_1 = grid.dx_1
        dy_1 = grid.dy_1
        dz_1 = grid.dz_1
    else:
        warnings.warn(
            "Assuming equidistant grid. To silence this warning, please pass a Pencil grid object instead of specifying dx,dy,..."
        )
        dx_1 = np.ones(np.size(f, -1)) / dx
        dy_1 = np.ones(np.size(f, -2)) / dy
        dz_1 = np.ones(np.size(f, -3)) / dz

    if coordinate_system == "cartesian":
        return xder(f[0], dx_1=dx_1) + yder(f[1], dy_1=dy_1) + zder(f[2], dz_1=dz_1)

    if coordinate_system == "cylindrical":
        if x is None:
            print("ERROR: need to specify x (radius) for cylindrical coordinates.")
            raise ValueError
        # Make sure x has compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        return (
            xder(x * f[0], dx_1=dx_1) / x
            + yder(f[1], dy_1=dy_1) / x
            + zder(f[2], dz_1=dz_1)
        )

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
            xder(x**2 * f[0], dx_1=dx_1) / x**2
            + yder(np.sin(y) * f[1], dy_1=dy_1) / (x * np.sin(y))
            + zder(f[2], dz_1=dz_1) / (x * np.sin(y))
        )

    print("ERROR: could not recognize coordinate system {0}".format(coordinate_system))
    raise ValueError


def grad(
    f,
    dx=None,
    dy=None,
    dz=None,
    x=None,
    y=None,
    coordinate_system="cartesian",
    grid=None,
):
    """
    grad(f, dx=None, dy=None, dz=None, x=None, y=None, coordinate_system="cartesian", grid=None)

    Take the gradient of a pencil code scalar array f in various coordinate systems.

    Parameters
    ----------
    f : ndarray
        Pencil code scalar array f.

    grid : pencil.read.grids.Grid
        Pencil grid object. See pc.read.grid().

    coordinate_system : string
        Coordinate system under which to take the divergence.
        Takes 'cartesian', 'cylindrical' and 'spherical'.

    Deprecated parameters (only for backwards compatibility)
    --------------------------------------------------------
    dx, dy, dz : floats
        Grid spacing in the three dimensions. These will not have any effect if grid!=None

    x, y : ndarrays
        Radial (x) and polar (y) coordinates, 1d arrays. These will not have any effect if grid!=None
    """

    import numpy as np
    from pencil.math.derivatives.der import xder, yder, zder

    if f.ndim != 3:
        print("grad: must have scalar 3-D array f[mz, my, mx] for gradient.")
        raise ValueError

    if grid is not None:
        x = grid.x
        y = grid.y
        dx_1 = grid.dx_1
        dy_1 = grid.dy_1
        dz_1 = grid.dz_1
    else:
        warnings.warn(
            "Assuming equidistant grid. To silence this warning, please pass a Pencil grid object instead of specifying dx,dy,..."
        )
        dx_1 = np.ones(np.size(f, -1)) / dx
        dy_1 = np.ones(np.size(f, -2)) / dy
        dz_1 = np.ones(np.size(f, -3)) / dz

    grad_value = np.zeros((3,) + f.shape)

    if coordinate_system == "cartesian":
        grad_value[0] = xder(f, dx_1=dx_1)
        grad_value[1] = yder(f, dy_1=dy_1)
        grad_value[2] = zder(f, dz_1=dz_1)

    if coordinate_system == "cylindrical":
        if x is None:
            print("ERROR: need to specify x (radius) for cylindrical coordinates.")
            raise ValueError
        # Make sure x has compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        grad_value[0] = xder(f, dx_1=dx_1)
        grad_value[1] = yder(f, dy_1=dy_1) / x
        grad_value[2] = zder(f, dz_1=dz_1)

    if coordinate_system == "spherical":
        if (x is None) or (y is None):
            print(
                "ERROR: need to specify x (radius) and y (polar angle) for spherical coordinates."
            )
            raise ValueError
        # Make sure x and y have compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        y = y[np.newaxis, :, np.newaxis]
        grad_value[0] = xder(f, dx_1=dx_1)
        grad_value[1] = yder(f, dy_1=dy_1) / x
        grad_value[2] = zder(f, dz_1=dz_1) / (x * np.sin(y))

    return grad_value


def curl(
    f,
    dx=None,
    dy=None,
    dz=None,
    x=None,
    y=None,
    run2D=False,
    coordinate_system="cartesian",
    grid=None,
):
    """
    curl(f, dx=None, dy=None, dz=None, x=None, y=None, run2D=False, coordinate_system="cartesian", grid=None)

    Take the curl of a pencil code vector array f in various coordinate systems.

    Parameters
    ----------
    f : ndarray
        Pencil code scalar array f.

    grid : pencil.read.grids.Grid
        Pencil grid object. See pc.read.grid().

    run2D : bool
        Deals with pure 2-D snapshots.
        !Only for Cartesian grids at the moment!
        Requires grid!=None.

    coordinate_system : string
        Coordinate system under which to take the divergence.
        Takes 'cartesian', 'cylindrical' and 'spherical'.
        !Does not work for 2d runs yet!

    Deprecated parameters (only for backwards compatibility)
    --------------------------------------------------------
    dx, dy, dz : floats
        Grid spacing in the three dimensions. These will not have any effect if grid!=None

    x, y : ndarrays
        Radial (x) and polar (y) coordinates, 1d arrays. These will not have any effect if grid!=None
    """

    import numpy as np
    from pencil.math.derivatives.der import xder, yder, zder

    if f.shape[0] != 3:
        print("curl: must have vector 4-D array f[3, mz, my, mx] for curl.")
        raise ValueError

    if grid is not None:
        x = grid.x
        y = grid.y
        dx_1 = grid.dx_1
        dy_1 = grid.dy_1
        dz_1 = grid.dz_1
    else:
        warnings.warn(
            "Assuming equidistant grid. To silence this warning, please pass a Pencil grid object instead of specifying dx,dy,..."
        )
        dx_1 = np.ones(np.size(f, -1)) / dx
        dy_1 = np.ones(np.size(f, -2)) / dy
        dz_1 = np.ones(np.size(f, -3)) / dz

    if run2D:
        if grid is None:
            raise NotImplementedError(
                "run2D=True supported only when a grid object is passed"
            )

        if coordinate_system != "cartesian":
            """
            KG: I haven't checked if it works correctly for non-Cartesian coordinates, but it looks like it should work.
            TODO
            """
            warnings.warn(
                "run2D=True may not work correctly for non-Cartesian coordinates.",
                RuntimeWarning,
            )

        # We insert a dummy axis so that we can use the same code as in the 3D case for calculating the curl.
        if grid.dx == grid.Lx:
            # 2-D case in the yz-plane
            newax_size = np.size(grid.dx_1)
            f = np.stack((f,) * newax_size, -1)
        elif grid.dy == grid.Ly:
            # 2-D case in the xz-plane
            newax_size = np.size(grid.dy_1)
            f = np.stack((f,) * newax_size, -2)
        elif grid.dz == grid.Lz:
            # 2-D case in the xy-plane
            newax_size = np.size(grid.dz_1)
            f = np.stack((f,) * newax_size, -3)
        else:
            raise RuntimeError("Unable to detect which axis was omitted in the 2D run.")

    curl_value = np.zeros_like(f)

    if coordinate_system == "cartesian":
        curl_value[0] = yder(f[2], dy_1=dy_1) - zder(f[1], dz_1=dz_1)
        curl_value[1] = zder(f[0], dz_1=dz_1) - xder(f[2], dx_1=dx_1)
        curl_value[2] = xder(f[1], dx_1=dx_1) - yder(f[0], dy_1=dy_1)
    if coordinate_system == "cylindrical":
        if x is None:
            print("ERROR: need to specify x (radius) for cylindrical coordinates.")
            raise ValueError
        # Make sure x has compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        curl_value[0] = (1 / x) * yder(f[2], dy_1=dy_1) - zder(f[1], dz_1=dz_1)
        curl_value[1] = zder(f[0], dz_1=dz_1) - xder(f[2], dx_1=dx_1)
        curl_value[2] = (1 / x) * xder(x * f[1], dx_1=dx_1) - (1 / x) * yder(
            f[0], dy_1=dy_1
        )
    if coordinate_system == "spherical":
        if (x is None) or (y is None):
            print(
                "ERROR: need to specify x (radius) and y (polar angle) for spherical coordinates."
            )
            raise ValueError
        # Make sure x and y have compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        y = y[np.newaxis, :, np.newaxis]
        curl_value[0] = (yder(np.sin(y) * f[2], dy_1=dy_1) - zder(f[1], dz_1=dz_1)) / (
            x * np.sin(y)
        )
        curl_value[1] = (
            zder(f[0], dz_1=dz_1) / np.sin(y) - xder(x * f[2], dx_1=dx_1)
        ) / x
        curl_value[2] = (xder(x * f[1], dx_1=dx_1) - yder(f[0], dy_1=dy_1)) / x

    if run2D:
        # Remove the dummy axis we inserted earlier
        if grid.dx == grid.Lx:
            curl_value = curl_value[..., 0]
        elif grid.dy == grid.Ly:
            curl_value = curl_value[..., 0, :]
        elif grid.dz == grid.Lz:
            curl_value = curl_value[..., 0, :, :]

    return curl_value


def curl2(
    f,
    dx=None,
    dy=None,
    dz=None,
    x=None,
    y=None,
    coordinate_system="cartesian",
    grid=None,
):
    """
    curl2(f, dx=None, dy=None, dz=None, x=None, y=None, coordinate_system="cartesian", grid=None)

    Take the double curl of a pencil code vector array f.

    Parameters
    ----------
    f : ndarray
        Pencil code vector array f.

    grid : pencil.read.grids.Grid
        Pencil grid object. See pc.read.grid().

    coordinate_system : string
        Coordinate system under which to take the divergence.
        Takes 'cartesian', 'cylindrical' and 'spherical'.

    Deprecated parameters (only for backwards compatibility)
    --------------------------------------------------------
    dx, dy, dz : floats
        Grid spacing in the three dimensions. These will not have any effect if grid!=None

    x, y : ndarrays
        Radial (x) and polar (y) coordinates, 1d arrays. These will not have any effect if grid!=None
    """

    import numpy as np
    from pencil.math.derivatives.der import xder, yder, zder, xder2, yder2, zder2

    if f.ndim != 4 or f.shape[0] != 3:
        print("curl2: must have vector 4-D array f[3, mz, my, mx] for curl2.")
        raise ValueError

    if grid is not None:
        x = grid.x
        y = grid.y
        dx_1 = grid.dx_1
        dy_1 = grid.dy_1
        dz_1 = grid.dz_1
        dx_tilde = grid.dx_tilde
        dy_tilde = grid.dy_tilde
        dz_tilde = grid.dz_tilde
    else:
        warnings.warn(
            "Assuming equidistant grid. To silence this warning, please pass a Pencil grid object instead of specifying dx,dy,..."
        )
        dx_1 = np.ones(np.size(f, -1)) / dx
        dy_1 = np.ones(np.size(f, -2)) / dy
        dz_1 = np.ones(np.size(f, -3)) / dz
        dx_tilde = np.zeros(np.size(f, -1))
        dy_tilde = np.zeros(np.size(f, -2))
        dz_tilde = np.zeros(np.size(f, -3))

    curl2_value = np.zeros(f.shape)

    if coordinate_system == "cartesian":
        curl2_value[0] = (
            xder(yder(f[1], dy_1=dy_1) + zder(f[2], dz_1=dz_1), dx_1=dx_1)
            - yder2(f[0], dy_1=dy_1, dy_tilde=dy_tilde)
            - zder2(f[0], dz_1=dz_1, dz_tilde=dz_tilde)
        )
        curl2_value[1] = (
            yder(xder(f[0], dx_1=dx_1) + zder(f[2], dz_1=dz_1), dy_1=dy_1)
            - xder2(f[1], dx_1=dx_1, dx_tilde=dx_tilde)
            - zder2(f[1], dz_1=dz_1, dz_tilde=dz_tilde)
        )
        curl2_value[2] = (
            zder(xder(f[0], dx_1=dx_1) + yder(f[1], dy_1=dy_1), dz_1=dz_1)
            - xder2(f[2], dx_1=dx_1, dx_tilde=dx_tilde)
            - yder2(f[2], dy_1=dy_1, dy_tilde=dy_tilde)
        )
    if coordinate_system == "cylindrical":
        if x is None:
            print("ERROR: need to specify x (radius) for cylindrical coordinates.")
            raise ValueError
        # Make sure x has compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        curl2_value[0] = (
            yder(f[1], dy_1=dy_1) / x**2
            + xder(yder(f[1], dy_1=dy_1), dx_1=dx_1) / x
            - yder2(f[0], dy_1=dy_1, dy_tilde=dy_tilde) / x**2
            - zder2(f[0], dz_1=dz_1, dz_tilde=dz_tilde)
            + xder(zder(f[2], dz_1=dz_1), dx_1=dx_1)
        )
        curl2_value[1] = (
            yder(zder(f[2], dz_1=dz_1), dy_1=dy_1) / x
            - zder2(f[1], dz_1=dz_1, dz_tilde=dz_tilde)
            + f[1] / x**2
            - xder(f[1], dx_1=dx_1) / x
            - xder2(f[1], dx_1=dx_1, dx_tilde=dx_tilde)
            + xder(yder(f[0], dy_1=dy_1), dx_1=dx_1) / x
            - yder(f[0], dy_1=dy_1) / x**2
        )
        curl2_value[2] = (
            zder(f[0], dz_1=dz_1) / x
            + xder(zder(f[0], dz_1=dz_1), dx_1=dx_1)
            - zder(f[2], dx_1=dx_1) / x
            - xder2(f[2], dx_1=dx_1, dx_tilde=dx_tilde)
            - yder2(f[2], dy_1=dy_1, dy_tilde=dy_tilde) / x**2
            + yder(zder(f[1], dz_1=dz_1), dy_1=dy_1) / x
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
            yder(
                np.sin(y) * (xder(x * f[1], dx_1=dx_1) - yder(f[0], dy_1=dy_1)) / x,
                dy_1=dy_1,
            )
            - zder(
                (zder(f[0], dz_1=dz_1) / np.sin(y) - xder(x * f[2], dx_1=dx_1)) / x,
                dz_1=dz_1,
            )
        ) / (x * np.sin(y))
        curl2_value[1] = (
            zder(
                (yder(np.sin(y) * f[2], dy_1=dy_1) - zder(f[1], dz_1=dz_1))
                / (x * np.sin(y)),
                dz_1=dz_1,
            )
            / np.sin(y)
            - xder((xder(x * f[1], dx_1=dx_1) - yder(f[0], dy_1=dy_1)), dx_1=dx_1)
        ) / x
        curl2_value[2] = (
            xder(
                (zder(f[0], dz_1=dz_1) / np.sin(y) - xder(x * f[2], dx_1=dx_1)),
                dx_1=dx_1,
            )
            - yder(
                (yder(np.sin(y) * f[2], dy_1=dy_1) - zder(f[1], dz_1=dz_1))
                / (x * np.sin(y)),
                dy_1=dy_1,
            )
        ) / x
    return curl2_value


def del2(
    f,
    dx=None,
    dy=None,
    dz=None,
    x=None,
    y=None,
    coordinate_system="cartesian",
    grid=None,
):
    """
    del2(f, dx=None, dy=None, dz=None, x=None, y=None, coordinate_system="cartesian", grid=None)

    Calculate del2, the Laplacian of a scalar field f.

    Parameters
    ----------
    f : ndarray
        Pencil code vector array f.

    grid : pencil.read.grids.Grid
        Pencil grid object. See pc.read.grid().

    coordinate_system : string
        Coordinate system under which to take the divergence.
        Takes 'cartesian', 'cylindrical' and 'spherical'.

    Deprecated parameters (only for backwards compatibility)
    --------------------------------------------------------
    dx, dy, dz : floats
        Grid spacing in the three dimensions. These will not have any effect if grid!=None

    x, y : ndarrays
        Radial (x) and polar (y) coordinates, 1d arrays. These will not have any effect if grid!=None
    """

    import numpy as np
    from pencil.math.derivatives.der import xder2, yder2, zder2, xder, yder

    if grid is not None:
        x = grid.x
        y = grid.y
        dx_1 = grid.dx_1
        dy_1 = grid.dy_1
        dz_1 = grid.dz_1
        dx_tilde = grid.dx_tilde
        dy_tilde = grid.dy_tilde
        dz_tilde = grid.dz_tilde
    else:
        warnings.warn(
            "Assuming equidistant grid. To silence this warning, please pass a Pencil grid object instead of specifying dx,dy,..."
        )
        dx_1 = np.ones(np.size(f, -1)) / dx
        dy_1 = np.ones(np.size(f, -2)) / dy
        dz_1 = np.ones(np.size(f, -3)) / dz
        dx_tilde = np.zeros(np.size(f, -1))
        dy_tilde = np.zeros(np.size(f, -2))
        dz_tilde = np.zeros(np.size(f, -3))

    if coordinate_system == "cartesian":
        del2_value = (
            xder2(f, dx_1=dx_1, dx_tilde=dx_tilde)
            + yder2(f, dy_1=dy_1, dy_tilde=dy_tilde)
            + zder2(f, dz_1=dz_1, dz_tilde=dz_tilde)
        )
    if coordinate_system == "cylindrical":
        if x is None:
            print("ERROR: need to specify x (radius)")
            raise ValueError
        # Make sure x has compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        del2_value = (
            xder(f, dx_1=dx_1) / x
            + xder2(f, dx_1=dx_1, dx_tilde=dx_tilde)
            + yder2(f, dy_1=dy_1, dy_tilde=dy_tilde) / (x**2)
            + zder2(f, dz_1=dz_1, dz_tilde=dz_tilde)
        )
    if coordinate_system == "spherical":
        if x is None or y is None:
            print("ERROR: need to specify x (radius) and y (polar angle)")
            raise ValueError
        # Make sure x and y have compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        y = y[np.newaxis, :, np.newaxis]
        del2_value = (
            2 * xder(f, dx_1=dx_1) / x
            + xder2(f, dx_1=dx_1, dx_tilde=dx_tilde)
            + np.cos(y) * yder(f, dy_1=dy_1) / ((x**2) * np.sin(y))
            + yder2(f, dy_1=dy_1, dy_tilde=dy_tilde) / (x**2)
            + zder2(f, dz_1=dz_1, dz_tilde=dz_tilde) / ((x * np.sin(y)) ** 2)
        )
    return del2_value


def del2v(
    f,
    dx=None,
    dy=None,
    dz=None,
    x=None,
    y=None,
    coordinate_system="cartesian",
    grid=None,
):
    """
    del2v(f, dx=None, dy=None, dz=None, x=None, y=None, coordinate_system="cartesian", grid=None)

    Calculate del2, the Laplacian of a vector field f.

    Parameters
    ----------
    f : ndarray
        Pencil code vector array f.

    grid : pencil.read.grids.Grid
        Pencil grid object. See pc.read.grid().

    coordinate_system : string
      Coordinate system under which to take the divergence.
      Takes 'cartesian', 'cylindrical' and 'spherical'.

    Deprecated parameters (only for backwards compatibility)
    --------------------------------------------------------
    dx, dy, dz : floats
        Grid spacing in the three dimensions. These will not have any effect if grid!=None

    x, y : ndarrays
        Radial (x) and polar (y) coordinates, 1d arrays. These will not have any effect if grid!=None
    """

    import numpy as np
    from pencil.math.derivatives.der import xder2, yder2, zder2, yder, zder

    if f.shape[0] != 3:
        print(
            "Vector Laplacian: must have a vector 4D array f(3, mz, my, mx) for Vector Laplacian"
        )
        raise ValueError

    if grid is not None:
        x = grid.x
        y = grid.y
        dx_1 = grid.dx_1
        dy_1 = grid.dy_1
        dz_1 = grid.dz_1
        dx_tilde = grid.dx_tilde
        dy_tilde = grid.dy_tilde
        dz_tilde = grid.dz_tilde
        pass_to_del2 = {"grid": grid, "coordinate_system": coordinate_system}
    else:
        warnings.warn(
            "Assuming equidistant grid. To silence this warning, please pass a Pencil grid object instead of specifying dx,dy,..."
        )
        dx_1 = np.ones(np.size(f, -1)) / dx
        dy_1 = np.ones(np.size(f, -2)) / dy
        dz_1 = np.ones(np.size(f, -3)) / dz
        dx_tilde = np.zeros(np.size(f, -1))
        dy_tilde = np.zeros(np.size(f, -2))
        dz_tilde = np.zeros(np.size(f, -3))
        pass_to_del2 = {
            "dx": dx,
            "dy": dy,
            "dz": dz,
            "x": x,
            "y": y,
            "coordinate_system": coordinate_system,
        }

    del2v_value = np.zeros(f.shape)

    if coordinate_system == "cartesian":
        del2v_value[0] = (
            xder2(f[0], dx_1=dx_1, dx_tilde=dx_tilde)
            + yder2(f[0], dy_1=dy_1, dy_tilde=dy_tilde)
            + zder2(f[0], dz_1=dz_1, dz_tilde=dz_tilde)
        )
        del2v_value[1] = (
            xder2(f[1], dx_1=dx_1, dx_tilde=dx_tilde)
            + yder2(f[1], dy_1=dy_1, dy_tilde=dy_tilde)
            + zder2(f[1], dz_1=dz_1, dz_tilde=dz_tilde)
        )
        del2v_value[2] = (
            xder2(f[2], dx_1=dx_1, dx_tilde=dx_tilde)
            + yder2(f[2], dy_1=dy_1, dy_tilde=dy_tilde)
            + zder2(f[2], dz_1=dz_1, dz_tilde=dz_tilde)
        )
    if coordinate_system == "cylindrical":
        if x is None:
            print("Error: need to specify x (radius)")
            raise ValueError
        # Make sure x has compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        del2v_value[0] = (
            del2(f[0], **pass_to_del2)
            - f[0] / (x**2)
            - 2 * yder(f[1], dy_1=dy_1) / x**2
        )
        del2v_value[1] = (
            del2(f[1], **pass_to_del2)
            - f[1] / (x**2)
            + 2 * yder(f[0], dy_1=dy_1) / x**2
        )
        del2v_value[2] = del2(f[2], **pass_to_del2)
    if coordinate_system == "spherical":
        if x is None or y is None:
            print("ERROR: need to specify x (radius) and y (polar angle)")
        # Make sure x and y have compatible dimensions.
        x = x[np.newaxis, np.newaxis, :]
        y = y[np.newaxis, :, np.newaxis]
        del2v_value[0] = (
            del2(f[0], **pass_to_del2)
            - 2 * f[0] / x**2
            - 2 * yder(f[1], dy_1=dy_1) / x**2
            - 2 * np.cos(y) * f[1] / (x**2 * np.sin(y))
            - 2 * zder(f[2], dz_1=dz_1) / (x**2 * np.sin(y))
        )
        del2v_value[1] = (
            del2(f[1], **pass_to_del2)
            - f[1] / (x * np.sin(y)) ** 2
            + 2 * yder(f[0], dy_1=dy_1) / (x**2)
            - (2 * np.cos(y)) * zder(f[2], dz_1=dz_1) / (x * np.sin(y)) ** 2
        )
        del2v_value[2] = (
            del2(f[2], **pass_to_del2)
            - f[2] / (x * np.sin(y)) ** 2
            + 2 * zder(f[0], dz_1=dz_1) / (np.sin(y) * x**2)
            + (2 * np.cos(y)) * zder(f[1], dz_1=dz_1) / (x * np.sin(y)) ** 2
        )
    return del2v_value


def curl3(f, dx, dy, dz, x=None, y=None, coordinate_system="cartesian"):
    """
    curl3(f, dx, dy, dz, x=None, y=None, coordinate_system="cartesian")

    Take the triple curl of a pencil code vector array f. Supports only equidistant grids.

    Parameters
    ----------
    f : ndarray
        Pencil code vector array f.

    dx, dy, dz : floats
        Grid spacing in the three dimensions.

    x, y : ndarrays
        Radial (x) and polar (y) coordinates, 1d arrays.

    coordinate_system : string
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


def del6(f, dx=None, dy=None, dz=None, grid=None):
    """
    del6(f, dx=None, dy=None, dz=None, grid=None)

    Calculate del6 (defined here as d^6/dx^6 + d^6/dy^6 + d^6/dz^6, rather
    than del2^3) of a scalar f for hyperdiffusion.

    Parameters
    ----------
    f : ndarray
        Pencil code scalar array f.

    grid : pencil.read.grids.Grid
        Pencil grid object. See pc.read.grid().

    Deprecated parameters (only for backwards compatibility)
    --------------------------------------------------------
    dx, dy, dz : floats
        Grid spacing in the three dimensions. These will not have any effect if grid!=None
    """

    from pencil.math.derivatives.der import xder6, yder6, zder6
    import numpy as np

    if grid is not None:
        dx_1 = grid.dx_1
        dy_1 = grid.dy_1
        dz_1 = grid.dz_1
    else:
        warnings.warn(
            "Assuming equidistant grid. To silence this warning, please pass a Pencil grid object instead of specifying dx,dy,..."
        )
        dx_1 = np.ones(np.size(f, -1)) / dx
        dy_1 = np.ones(np.size(f, -2)) / dy
        dz_1 = np.ones(np.size(f, -3)) / dz

    return xder6(f, dx_1=dx_1) + yder6(f, dy_1=dy_1) + zder6(f, dz_1=dz_1)


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
