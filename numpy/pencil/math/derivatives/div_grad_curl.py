#
# $Id: div_grad_curl.py,v 1.1 2007-11-19 14:29:08 joishi Exp $
#
"""module to do div, grad, curl (but not 'all that') for pencil-code data.

"""
import numpy as N
from der import *

def div(f,dx,dy,dz):
    """
    take divervenge of pencil code vector array
    """
    if (f.ndim != 4):
        print "div: must have vector 4-D array f[mvar,mz,my,mx] for divergence"
        raise ValueError
   
    return xder(f[0,...],dx) + yder(f[1,...],dy) + zder(f[2,...],dz)

def grad(f,dx,dy,dz):
    """
    take the curl of a pencil code scalar array.
    """
    if (f.ndim != 3):
        print "grad: must have scalar 3-D array f[mz,my,mx] for gradient"
        raise ValueError
    
    grad = N.empty((3,)+f.shape)
    grad[0,...] = xder(f,dx)
    grad[1,...] = yder(f,dy)
    grad[2,...] = zder(f,dz)
    
    return grad


def curl(f,dx,dy,dz):
    """
    take the curl of a pencil code vector array.
    """
    if (f.ndim != 4 or f.shape[0] != 3):
        print "curl: must have vector 4-D array f[3,mz,my,mx] for curl"
        raise ValueError

    curl = N.empty(f.shape)
    curl[0,...] = yder(f[2,...],dy) - zder(f[1,...],dz)
    curl[1,...] = zder(f[0,...],dz) - xder(f[2,...],dx)
    curl[2,...] = xder(f[1,...],dx) - yder(f[0,...],dy)
    
    return curl

def curl2(f,dx,dy,dz):
    """
    take the double curl of a pencil code vector array.

    CARTESIAN COORDINATES ONLY!!
    """
    if (f.ndim != 4 or f.shape[0] != 3):
        print "curl2: must have vector 4-D array f[3,mz,my,mx] for curl2"
        raise ValueError

    curl2 = N.empty(f.shape)
    curl2[0,...] = xder(yder(f[1,...],dy) + zder(f[2,...],dz),dx) \
                   - yder2(f[0,...],dy) - zder2(f[0,...],dz)
    curl2[1,...] = yder(xder(f[0,...],dx) + zder(f[2,...],dz),dy) \
                   - xder2(f[1,...],dx) - zder2(f[1,...],dz)
    curl2[2,...] = zder(xder(f[0,...],dx) + yder(f[1,...],dy),dz) \
                   - xder2(f[2,...],dx) - yder2(f[2,...],dy)
    
    return curl2

