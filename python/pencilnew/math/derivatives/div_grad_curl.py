# div_grad_curl.py
#
# Contains the vector calculus derivatives.
"""
Compute the divergence, gradient and curl.
"""


def div(f, dx, dy, dz):
    """
    Take divervenge of pencil code vector array f.
    """

    from pencilnew.math.derivatives.der import xder, yder, zder

    if (f.ndim != 4):
        print("div: must have vector 4-D array f[mvar ,mz, my, mx] for divergence.")
        raise ValueError

    return xder(f[0, ...], dx) + yder(f[1, ...], dy) + zder(f[2, ...], dz)


def grad(f, dx, dy, dz):
    """
    Take the curl of a pencil code scalar array f.
    """

    import numpy as np
    from pencilnew.math.derivatives.der import xder, yder, zder

    if (f.ndim != 3):
        print("grad: must have scalar 3-D array f[mz, my, mx] for gradient.")
        raise ValueError

    grad_value = np.empty((3,) + f.shape)
    grad_value[0, ...] = xder(f, dx)
    grad_value[1, ...] = yder(f, dy)
    grad_value[2, ...] = zder(f, dz)

    return grad_value


def curl(f, dx, dy, dz, run2D=False):
    """
    Take the curl of a pencil code vector array.
    The run2D parameter deals with pure 2-D snapshots (solved the (x,z)-plane pb).
    """

    import numpy as np
    from pencilnew.math.derivatives.der import xder, yder, zder

    if (f.shape[0] != 3):
        print("curl: must have vector 4-D array f[3,mz,my,mx] for curl.")
        raise ValueError

    curl_value = np.empty_like(f)
    if (dy != 0. and dz != 0.):
        # 3-D case
        curl_value[0, ...] = yder(f[2, ...], dy) - zder(f[1, ...], dz)
        curl_value[1, ...] = zder(f[0, ...], dz) - xder(f[2, ...], dx)
        curl_value[2, ...] = xder(f[1, ...], dx) - yder(f[0, ...], dy)
    elif (dy == 0.):
        # 2-D case in the xz-plane
        curl_value[0, ...] = zder(f, dz, run2D)[0, ...] - xder(f, dx)[2, ...]
    else:
        # 2-D case in the xy-plane
        curl_value[0, ...] = xder(f, dx)[1, ...] - yder(f, dy)[0, ...]

    return curl_value


def curl2(f, dx, dy, dz):
    """
    Take the double curl of a pencil code vector array f.
    CARTESIAN COORDINATES ONLY!!
    """

    import numpy as np
    from pencilnew.math.derivatives.der import xder, yder, zder, xder2, yder2, zder2

    if (f.ndim != 4 or f.shape[0] != 3):
        print("curl2: must have vector 4-D array f[3,mz,my,mx] for curl2.")
        raise ValueError

    curl2_value = np.empty(f.shape)
    curl2_value[0, ...] = xder(yder(f[1, ...], dy) + zder(f[2, ...], dz), dx) \
                          -yder2(f[0, ...], dy) - zder2(f[0, ...], dz)
    curl2_value[1, ...] = yder(xder(f[0, ...], dx) + zder(f[2, ...], dz), dy) \
                          -xder2(f[1, ...], dx) - zder2(f[1, ...], dz)
    curl2_value[2, ...] = zder(xder(f[0, ...], dx) + yder(f[1, ...], dy), dz) \
                          -xder2(f[2, ...], dx) - yder2(f[2, ...], dy)

    return curl2


def del2(f, dx, dy, dz):
    """
    Calculate del6 (defined here as d^6/dx^6 + d^6/dy^6 + d^6/dz^6, rather
    than del2^3) of a scalar f for hyperdiffusion.
    """

    from pencilnew.math.derivatives.der import xder2, yder2, zder2

    return xder2(f, dx) + yder2(f, dy) + zder2(f, dz)


def del6(f, dx, dy, dz):
    """
    Calculate del6 (defined here as d^6/dx^6 + d^6/dy^6 + d^6/dz^6, rather
    than del2^3) of a scalar f for hyperdiffusion.
    """

    from pencilnew.math.derivatives.der import xder6, yder6, zder6

    return xder6(f, dx) + yder6(f, dy) + zder6(f, dz)
