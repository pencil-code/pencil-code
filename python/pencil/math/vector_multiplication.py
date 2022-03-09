# vector_multiplication.py
#
# Various vector multiplication routines.
"""
Contains the methods for the dot products and cross product.
"""


def dot(a, b):
    """
    dot(a, b)

    Take dot product of two pencil-code vectors a and b.

    Parameters
    ----------
    a, b : ndarrays
        Pencil-code vectors with shape [3, mz, my, mx].
    """

    import numpy as np

    if a.ndim != 4 or a.shape[0] != 3 or b.ndim != 4 or b.shape[0] != 3:
        print("error: both vectors must be 4-D array f[3, mz, my, mx] for dot.")
        raise ValueError

    return np.sum(a * b, axis=0)


def dot2(a):
    """
    dot2(a)

    Take dot product of a pencil-code vector with itself.

    Parameters
    ----------
    a : ndarray
        Pencil-code vector with shape [3, mz, my, mx].
    """

    return dot(a, a)


def cross(a, b):
    """
    cross(a, b)

    Take cross of two pencil-code vectors a and b.

    Parameters
    ----------
    a, b : ndarrays
        Pencil-code vectors with shape [3, mz, my, mx].
    """

    import numpy as np

    if a.ndim != 4 or a.shape[0] != 3 or a.shape != b.shape:
        print("error: both vectors must be 4-D array f[3, mz, my, mx] for cross.")
        raise ValueError

    # (a x b)_i = eps_ijk a_j b_k
    return np.cross(a, b, axis=0)
