
__version__ = "$Id$"
import numpy as N

def dot(a,b):
    """take dot product of two pencil-code vectors a & b with shabe

    a.shape = (3,mz,my,mx)

    """

    if (a.ndim != 4 or a.shape[0] != 3 or
        b.ndim != 4 or b.shape[0] != 3):
        print "dot: both vectors must be 4-D array f[3,mz,my,mx] for dot"
        raise ValueError

    return a[0,...]*b[0,...]+ a[1,...]*b[1,...] + a[2,...]*b[2,...]

def dot2(a):
    """take dot product of a pencil-code vector with itself.

    a.shape = (3,mz,my,mx)

    """
    dot2 = dot(a,a)

    return dot2

def cross(a,b):
    """take cross of two pencil-code vectors a & b with shape

    a.shape = (2,mz,my,mx)

    """
    if (a.ndim != 4 or a.shape[0] != 3 or
        a.shape != b.shape):
        print "cross: both vectors must be 4-D array f[3,mz,my,mx] for dot"
        raise ValueError
    
    cross = N.empty(a.shape)

    # a x b = eps_ijk a_j b_k
    cross[0,...] = a[1,...]*b[2,...] - a[2,...]*b[1,...]
    cross[1,...] = a[2,...]*b[0,...] - a[0,...]*b[2,...]
    cross[2,...] = a[0,...]*b[1,...] - a[1,...]*b[0,...]
    
    return cross

    
