#
# $Id$
#
"""module to do div, grad, curl (but not 'all that') for pencil-code data.

"""
import numpy as N
from pencil.math.derivatives.der import *
from sys import exit
from pencil.files.param import read_param
from pencil.files.grid import read_grid
from pencil.files.dim import read_dim


def grad(f,dx,dy,dz,x=[],y=[],z=[]):
    """
    take the gradient of a pencil code scalar array.
    """
    if (f.ndim != 3):
        print("grad: must have scalar 3-D array f[mz,my,mx] for gradient")
        raise ValueError
    gd  = read_grid(quiet=True)
    if len(x) < 1:
        x = gd.x
    if len(y) < 1:
        y = gd.y
    if len(z) < 1:
        z = gd.z

    grad = N.empty((3,)+f.shape)
    grad[0,...] = xder(f,dx,x=x,y=y,z=z)
    grad[1,...] = yder(f,dy,x=x,y=y,z=z)
    grad[2,...] = zder(f,dz,x=x,y=y,z=z)

    return grad

def div(f,dx,dy,dz,x=[],y=[],z=[]):
    """
    take divergence of pencil code vector array
    """
    if (f.ndim != 4):
        print("div: must have vector 4-D array f[mvar,mz,my,mx] for divergence")
        raise ValueError
    param = read_param(quiet=True)
    gd  = read_grid(quiet=True)
    if len(x) < 1:
        x = gd.x
    if len(y) < 1:
        y = gd.y
    if len(z) < 1:
        z = gd.z

    div = xder(f[0,...],dx,x=x,y=y,z=z) +\
          yder(f[1,...],dy,x=x,y=y,z=z) +\
          zder(f[2,...],dz,x=x,y=y,z=z)

    if param.coord_system == 'cylindric':
        div += f[0,...]/x
    if param.coord_system == 'spherical':
        sin_y = N.sin(y)
        cos_y = N.cos(y)
        i_sin = N.where(N.abs(sin_y) < 1e-5)[0]
        if i_sin.size > 0:
            cos_y[i_sin] = 0.; sin_y[i_sin] = 1
        x_1, cotth = N.meshgrid(1./gd.x, cos_y/sin_y)
        div += 2*f[0,...]*x_1 + f[1,...]*x_1*cotth
   
    return div

def laplacian(f,dx,dy,dz,x=[],y=[],z=[]):
    """
    take the laplacian of a pencil code scalar array
    """
    param = read_param(quiet=True)
    gd  = read_grid(quiet=True)
    if len(x) < 1:
        x = gd.x
    if len(y) < 1:
        y = gd.y
    if len(z) < 1:
        z = gd.z

    laplacian = N.empty(f.shape)
    laplacian = xder2(f,dx,x=x,y=y,z=z)+\
                yder2(f,dy,x=x,y=y,z=z)+\
                zder2(f,dz,x=x,y=y,z=z)

    if param.coord_system == 'cylindric':
        laplacian += xder(f,dx,x=x,y=y,z=z)/x
    if param.coord_system == 'spherical':
        sin_y = N.sin(y)
        cos_y = N.cos(y)
        i_sin = N.where(N.abs(sin_y) < 1e-5)[0]
        if i_sin.size > 0:
            cos_y[i_sin] = 0.; sin_y[i_sin] = 1
        x_2, cotth = N.meshgrid(1./x**2, cos_y/sin_y)
        laplacian += 2*xder(f,dx,x=x,y=y,z=z)/x +\
                       yder(f,dy,x=x,y=y,z=z)*x_2*cotth

    return laplacian

def curl(f,dx,dy,dz,x=[],y=[],z=[],run2D=False):
    """
    take the curl of a pencil code vector array.
    23-fev-2009/dintrans+morin: introduced the run2D parameter to deal
    with pure 2-D snapshots (solved the (x,z)-plane pb)
    """
    if (f.shape[0] != 3):
        print("curl: must have vector 4-D array f[3,mz,my,mx] for curl")
        raise ValueError
    param = read_param(quiet=True)
    gd  = read_grid(quiet=True)
    if len(x) < 1:
        x = gd.x
    if len(y) < 1:
        y = gd.y
    if len(z) < 1:
        z = gd.z

    curl = N.empty_like(f)
    if (dy != 0. and dz != 0.):
    # 3-D case
        curl[0,...] = yder(f[2,...],dy,x=x,y=y,z=z) -\
                      zder(f[1,...],dz,x=x,y=y,z=z)
        curl[1,...] = zder(f[0,...],dz,x=x,y=y,z=z) -\
                      xder(f[2,...],dx,x=x,y=y,z=z)
        curl[2,...] = xder(f[1,...],dx,x=x,y=y,z=z) -\
                      yder(f[0,...],dy,x=x,y=y,z=z)
    elif (dy == 0.):
    # 2-D case in the (x,z)-plane
    # f[...,nz,1,nx] if run2D=False or f[...,nz,nx] if run2D=True
        curl[0,...] = zder(f,dz,x=x,y=y,z=z,run2D=run2D)[0,...] -\
                      xder(f,dx,x=x,y=y,z=z)[2,...]
    else:
    # 2-D case in the (x,y)-plane
    # f[...,1,ny,nx] if run2D=False or f[...,ny,nx] if run2D=True
        curl[0,...] = xder(f,dx,x=x,y=y,z=z)[1,...] -\
                      yder(f,dy,x=x,y=y,z=z)[0,...]

    if param.coord_system == 'cylindric':
        curl[2,...] += f[1,...]/x
    if param.coord_system == 'spherical':
        sin_y = N.sin(y)
        cos_y = N.cos(y)
        i_sin = N.where(N.abs(sin_y) < 1e-5)[0]
        if i_sin.size > 0:
            cos_y[i_sin] = 0.; sin_y[i_sin] = 1
        x_1, cotth = N.meshgrid(1./x, cos_y/sin_y)
        curl[0,...] += f[2,...]*x_1*cotth
        curl[1,...] -= f[2,...]/x
        curl[2,...] += f[1,...]/x

    return curl

def curl2(f,dx,dy,dz,x=[],y=[],z=[]):
    """
    take the double curl of a pencil code vector array.
    """
    if (f.ndim != 4 or f.shape[0] != 3):
        print("curl2: must have vector 4-D array f[3,mz,my,mx] for curl2")
        raise ValueError
    param = read_param(quiet=True)
    gd  = read_grid(quiet=True)
    if len(x) < 1:
        x = gd.x
    if len(y) < 1:
        y = gd.y
    if len(z) < 1:
        z = gd.z

    curl2 = N.empty(f.shape)
    curl2[0,...] = xder(yder(f[1,...],dy,x=x,y=y,z=z) +
                        zder(f[2,...],dz,x=x,y=y,z=z),dx,x=x,y=y,z=z) -\
                   yder2(f[0,...],dy,x=x,y=y,z=z) -\
                   zder2(f[0,...],dz,x=x,y=y,z=z)
    curl2[1,...] = yder(xder(f[0,...],dx,x=x,y=y,z=z) +
                        zder(f[2,...],dz,x=x,y=y,z=z),dy,x=x,y=y,z=z) -\
                   xder2(f[1,...],dx,x=x,y=y,z=z) -\
                   zder2(f[1,...],dz,x=x,y=y,z=z)
    curl2[2,...] = zder(xder(f[0,...],dx,x=x,y=y,z=z) +
                        yder(f[1,...],dy,x=x,y=y,z=z),dz,x=x,y=y,z=z) -\
                   xder2(f[2,...],dx,x=x,y=y,z=z) -\
                   yder2(f[2,...],dy,x=x,y=y,z=z)

    if param.coord_system == 'cylindric':
        curl2[0,...] +=                    yder(f[1,...],dy,x=x,y=y,z=z)/x**2
        curl2[1,...] += f[1,...]/gd.x**2 - xder(f[1,...],dx,x=x,y=y,z=z)/x
        curl2[2,...] +=                   (zder(f[0,...],dz,x=x,y=y,z=z) -
                                           xder(f[2,...],dx,x=x,y=y,z=z))/x
    if param.coord_system == 'spherical':
        sin_y = N.sin(y)
        cos_y = N.cos(y)
        i_sin = N.where(N.abs(sin_y) < 1e-5)[0]
        if i_sin.size > 0:
            cos_y[i_sin] = 0.; sin_y[i_sin] = 1
        x_1 ,cotth  = N.meshgrid(   1./x, cos_y/sin_y)
        sin2th, x_2 = N.meshgrid(1./x**2, 1/sin_y**2 )
        curl2[0,...] += (yder(f[1,...],dy,x=x,y=y,z=z) +
                         zder(f[2,...],dz,x=x,y=y,z=z))/x +\
              x_1*cotth*(xder(f[1,...],dx,x=x,y=y,z=z) -
                         yder(f[0,...],dy,x=x,y=y,z=z) + f[1,...]/x )
        curl2[1,...] +=  zder(f[2,...],dz,x=x,y=y,z=z)*x_1*cotth -\
                       2*xder(f[1,...],dx,x=x,y=y,z=z)/x
        curl2[2,...] += x_2*sin2th*f[2,...] - \
                       2*xder(f[2,...],dx,x=x,y=y,z=z)/x - (
                         yder(f[2,...],dy,x=x,y=y,z=z) +
                         zder(f[1,...],dz,x=x,y=y,z=z))*x_1*cotth

    return curl2

def del2(f,dx,dy,dz,x=[],y=[],z=[]):
    """taken from pencil code's sub.f90 
    !  calculate del6 (defined here as d^6/dx^6 + d^6/dy^6 + d^6/dz^6, rather
    !  than del2^3) of a scalar for hyperdiffusion
    Duplcation of laplacian why? Fred - added curvelinear
    """
    param = read_param(quiet=True)
    gd  = read_grid(quiet=True)
    if len(x) < 1:
        x = gd.x
    if len(y) < 1:
        y = gd.y
    if len(z) < 1:
        z = gd.z

    del2 =        xder2(f,dx,x=x,y=y,z=z)
    del2 = del2 + yder2(f,dy,x=x,y=y,z=z)
    del2 = del2 + zder2(f,dz,x=x,y=y,z=z)

    if param.coord_system == 'cylindric':
        del2 += xder(f,dx,x=x,y=y,z=z)/x
    if param.coord_system == 'spherical':
        sin_y = N.sin(y)
        cos_y = N.cos(y)
        i_sin = N.where(N.abs(sin_y) < 1e-5)[0]
        if i_sin.size > 0:
            cos_y[i_sin] = 0.; sin_y[i_sin] = 1
        x_2, cotth = N.meshgrid(1./x**2, cos_y/sin_y)
        del2 += 2*xder(f,dx,x=x,y=y,z=z)/x +\
                  yder(f,dy,x=x,y=y,z=z)*x_2*cotth

    return del2

def del6(f,dx,dy,dz,x=[],y=[],z=[]):
    """taken from pencil code's sub.f90 
    !  calculate del6 (defined here as d^6/dx^6 + d^6/dy^6 + d^6/dz^6, rather
    !  than del2^3) of a scalar for hyperdiffusion
    """
    gd  = read_grid(quiet=True)
    if len(x) < 1:
        x = gd.x
    if len(y) < 1:
        y = gd.y
    if len(z) < 1:
        z = gd.z
    del6 =        xder6(f,dx,x=x,y=y,z=z)
    del6 = del6 + yder6(f,dy,x=x,y=y,z=z)
    del6 = del6 + zder6(f,dz,x=x,y=y,z=z)

    return del6
