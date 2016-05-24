
__version__ = "$Id$"
import numpy as N

def dot(a,b):
    """take dot product of two pencil-code vectors a & b with shabe

    a.shape = (3,mz,my,mx) or avers = (3,:,:), (3,:)

    22-apr-16/fred: generalised to accept 1D and 2D vector arrays

    """

    if (a.shape[0] != 3 or a.ndim<2 or 
        b.ndim != a.ndim or b.shape != a.shape):
        print("dot: both must be vector arrays of same dimension f[3,mz,my,mx], 2D or 3D")
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

    a.shape = (2,mz,my,mx) or avers = (3,:,:), (3,:)

    22-apr-16/fred: generalised to accept 1D and 2D vector arrays

    """
    if ( a.shape[0] != 3 or a.shape != b.shape or a.ndim<2):
        print("cross: both must be vector arrays of same dimension f[3,mz,my,mx], 2D or 1D for dot")
        raise ValueError
    
    cross = N.empty(a.shape)

    # a x b = eps_ijk a_j b_k
    cross[0,...] = a[1,...]*b[2,...] - a[2,...]*b[1,...]
    cross[1,...] = a[2,...]*b[0,...] - a[0,...]*b[2,...]
    cross[2,...] = a[0,...]*b[1,...] - a[1,...]*b[0,...]
    
    return cross

    
