#
# $Id$
#
"""module to do div, grad, curl (but not 'all that') for pencil-code data.

"""
import numpy as N
from der import *
from sys import exit

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

def curl(f,dx,dy,dz,run2D=False):
    """
    take the curl of a pencil code vector array.
    23-fev-2009/dintrans+morin: introduced the run2D parameter to deal
    with pure 2-D snapshots (solved the (x,z)-plane pb)
    """
    if (f.shape[0] != 3):
        print "curl: must have vector 4-D array f[3,mz,my,mx] for curl"
        raise ValueError

    curl = N.empty_like(f)
    if (dy != 0. and dz != 0.):
# 3-D case
      curl[0,...] = yder(f[2,...],dy) - zder(f[1,...],dz)
      curl[1,...] = zder(f[0,...],dz) - xder(f[2,...],dx)
      curl[2,...] = xder(f[1,...],dx) - yder(f[0,...],dy)
    elif (dy == 0.):
# 2-D case in the (x,z)-plane
# f[...,nz,1,nx] if run2D=False or f[...,nz,nx] if run2D=True
      curl[0,...] = zder(f,dz,run2D)[0,...] - xder(f,dx)[2,...]
    else:
# 2-D case in the (x,y)-plane
# f[...,1,ny,nx] if run2D=False or f[...,ny,nx] if run2D=True
      curl[0,...] = xder(f,dx)[1,...] - yder(f,dy)[0,...]
    
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

def del2(f,dx,dy,dz):
    """taken from pencil code's sub.f90 
    !  calculate del6 (defined here as d^6/dx^6 + d^6/dy^6 + d^6/dz^6, rather
    !  than del2^3) of a scalar for hyperdiffusion
    """
    del2 = xder2(f,dx)
    del2 = del2 + yder2(f,dy)
    del2 = del2 + zder2(f,dz)

    return del2

def del6(f,dx,dy,dz):
    """taken from pencil code's sub.f90 
    !  calculate del6 (defined here as d^6/dx^6 + d^6/dy^6 + d^6/dz^6, rather
    !  than del2^3) of a scalar for hyperdiffusion
    """
    del6 = xder6(f,dx)
    del6 = del6 + yder6(f,dy)
    del6 = del6 + zder6(f,dz)

    return del6
