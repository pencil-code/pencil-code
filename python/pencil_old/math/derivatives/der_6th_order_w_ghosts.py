# $Id$
#
#
# calculate 6th order spatial derivatives with ghost zones assumed to be in the arrays
#
# Author: J. Oishi (joishi@amnh.org). based on A. Brandenburg's IDL routines
# 
#

#import numpy as N

"""
6th Order derivatives. currently only equidistant grids are supported.
"""
import numpy as N
from pencil_old.files.param import read_param
from pencil_old.files.grid import read_grid
from pencil_old.files.dim import read_dim

def xder_6th(f,dx,x=[],y=[],z=[],param=[],dim=[]):

    if (f.ndim != 3 and f.ndim != 4):
        print("%s dimension arrays not handled." % (str(f.ndim)))
        raise ValueError

    if not param:
        param=read_param(quiet=True)
    if not dim:
        dim=read_dim()
    if len(x) < 1:
        gd  = read_grid(quiet=True)
        x = gd.x
    dx=N.gradient(x)
    if (dim.nx!=1):
        dx2 = 1./(60.*dx)
    dfdx = N.zeros_like(f)
    l1 = 3
    l2 = f.shape[-1]-3
    if (l2 > l1 and dim.nx!=1):
        for l in range(l1,l2):
            dfdx[...,l] = dx2[l]*( +45.*(f[...,l+1]-f[...,l-1])
                                    -9.*(f[...,l+2]-f[...,l-2])
                                    +   (f[...,l+3]-f[...,l-3]) )
    else:
        dfdx = 0.
    return dfdx

def yder_6th(f,dy,x=[],y=[],z=[],param=[],dim=[]):

    if (f.ndim != 3 and f.ndim != 4):
        print("%s dimension arrays not handled." % (str(f.ndim)))
        raise ValueError

    if not param:
        param=read_param(quiet=True)
    if not dim:
        dim=read_dim()

    if len(y) < 1:
        gd  = read_grid(quiet=True)
        y = gd.y
    dy=N.gradient(y)
    if (dim.ny!=1):
        dy2 = 1./(60.*dy)
    dfdy = N.zeros_like(f)
    m1 = 3
    m2 = f.shape[-2]-3

    if (m2 > m1 and dim.ny != 1):
        for m in range(m1,m2):
            dfdy[...,m,:] = dy2[m]*( +45.*(f[...,m+1,:]-f[...,m-1,:]) 
                                      -9.*(f[...,m+2,:]-f[...,m-2,:]) 
                                      +   (f[...,m+3,:]-f[...,m-3,:]) )
    else:
        dfdy = 0.
    if param.coord_system == ('cylindric' or 'spherical'):
        if len(x) < 1:
            gd=read_grid(quiet=True)
            x=gd.x
        dfdy /= x
        
    return dfdy

def zder_6th(f,dz,x=[],y=[],z=[],run2D=False,param=[],dim=[]):
    
    if (f.ndim != 3 and f.ndim != 4):
        print("%s dimension arrays not handled." % (str(f.ndim)))
        raise ValueError

    if not param:
        param=read_param(quiet=True)
    if not dim:
        dim=read_dim()

    if len(z) < 1:
        gd  = read_grid(quiet=True)
        z = gd.z
    dz=N.gradient(z)
    if (dim.nz!=1):
        dz2 = 1./(60.*dz)
    dfdz = N.zeros_like(f)

    n1 = 3
    if run2D:
        n2 = f.shape[1]-3
    else:
        n2 = f.shape[-3]-3

    if (n2 > n1 and dim.nz!=1):
        if (run2D):
            # f[...,z,x] or f[...,z,y]
            for n in range(n1,n2):
                dfdz[...,n,:] = dz2[n]*(+45.*(f[...,n+1,:]-f[...,n-1,:])
                                         -9.*(f[...,n+2,:]-f[...,n-2,:])
                                            +(f[...,n+3,:]-f[...,n-3,:]) )

        else:
            # f[...,z,y,x]
            for n in range(n1,n2):
                dfdz[...,n,:,:] = dz2[n]*(+45.*(f[...,n+1,:,:]-f[...,n-1,:,:])
                                           -9.*(f[...,n+2,:,:]-f[...,n-2,:,:])
                                              +(f[...,n+3,:,:]-f[...,n-3,:,:]) )
    else:
        dfdz=0
    if param.coord_system == 'spherical':
        if (len(x) or len(y)) < 1:
            gd=read_grid(quiet=True)
            x=gd.x; y=gd.y
        sin_y = N.sin(y)
        siny1 = 1./sin_y
        i_sin = N.where(N.abs(sin_y) < 1e-5)[0]
        if i_sin.size > 0:
            siny1[i_sin] = 0.
        x_1, sin1th = N.meshgrid(1./x, siny1)
        dfdz *= x_1*sin1th

    return dfdz

def xder2_6th(f,dx,x=[],y=[],z=[],param=[],dim=[]):

    if (f.ndim != 3 and f.ndim != 4):
        print("%s dimension arrays not handled." % (str(f.ndim)))
        raise ValueError

    if not param:
        param=read_param(quiet=True)
    if not dim: 
        dim=read_dim()

    dx = N.gradient(x)
    if (dim.nx!=1):
        dx2 = 1./(180.*dx**2.)
    dfdx = N.zeros_like(f)
    l1 = 3
    l2 = f.shape[-1]-3

    if (l2 > l1 and dim.nx!=1):
        for l in range(l1,l2): 
            dfdx[...,l] = dx2[l]*(-490.* f[...,l]
                                  +270.*(f[...,l-1]+f[...,l+1])
                                  - 27.*(f[...,l-2]+f[...,l+2])
                                  +  2.*(f[...,l-3]+f[...,l+3]) )
    else:
        dfdx = 0.

    return dfdx

def yder2_6th(f,dy,x=[],y=[],z=[],param=[],dim=[]):

    if (f.ndim != 3 and f.ndim != 4):
        print("%s dimension arrays not handled." % (str(f.ndim)))
        raise ValueError
    
    if not param:
        param=read_param(quiet=True)
    if not dim:
        dim=read_dim()

    dy = N.gradient(y)
    if (dim.ny!=1):
        dy2 = 1./(180.*dy**2.)
    dfdy = N.zeros_like(f)
    m1 = 3
    m2 = f.shape[-2]-3

    if (m2 > m1 and dim.ny!=1):
        for m in range(m1,m2):
            dfdy[...,m,:] = dy2[m]*(-490.* f[...,m,:]
                                    +270.*(f[...,m-1,:]+f[...,m+1,:])
                                    - 27.*(f[...,m-2,:]+f[...,m+2,:])
                                    +  2.*(f[...,m-3,:]+f[...,m+3,:]) )
    else:
        dfdy = 0.
    if param.coord_system == ('cylindric' or 'spherical'):
        if (len(x) or len(y)) < 1:
            gd=read_grid(quiet=True)
            x=gd.x; y=gd.y
        dfdy /= x**2

    return dfdy

def zder2_6th(f,dz,x=[],y=[],z=[],param=[],dim=[]):

    if (f.ndim != 3 and f.ndim != 4):
        print("%s dimension arrays not handled." % (str(f.ndim)))
        raise ValueError

    if (len(z) < 1):
        gd=read_grid(quiet=True)
        z=gd.z
    if not param:
        param=read_param(quiet=True)
    if not dim:
        dim=read_dim()

    dz = N.gradient(z)
    if (dim.nz!=1):
        dz2 = 1./(180.*dz**2.)
    dfdz = N.zeros_like(f)
    n1 = 3
    n2 = f.shape[-3]-3
    
    if (n2 > n1 and dim.nz!=1):
        for n in range(n1,n2): 
            dfdz[...,n,:,:] = dz2[n]*(-490.* f[...,n,:,:]
                                      +270.*(f[...,n-1,:,:]+f[...,n+1,:,:])
                                      - 27.*(f[...,n-2,:,:]+f[...,n+2,:,:])
                                      +  2.*(f[...,n-3,:,:]+f[...,n+3,:,:]) )
    else:
        dfdz = 0.
    if param.coord_system == 'spherical':
        if (len(x) or len(y)) < 1:
            gd=read_grid(quiet=True)
            x=gd.x; y=gd.y
        sin_y = N.sin(y)
        siny1 = 1./sin_y
        i_sin = N.where(N.abs(sin_y) < 1e-5)[0]
        if i_sin.size > 0:
            siny1[i_sin] = 0.
        x_2, sin2th = N.meshgrid(1./x**2, siny1**2)
        dfdz *= x_2*sin2th

    return dfdz

def xder6_6th(f,dx,x=[],y=[],z=[]):
    if (f.ndim != 3 and f.ndim != 4):
        print("%s dimension arrays not handled." % (str(f.ndim)))
        raise ValueError

    l1 = 3
    l2 = f.shape[-1] - 3
    dx = N.gradient(x)
    if (dim.nx!=1):
        fac=1/dx**6
    d6fdx = N.zeros_like(f)

    dim=read_dim()

    if (l2 > l1 and dim.nx!=1):
        for l in range(l1,l2): 
            d6fdx[...,l] = fac[l]*( - 20.0* f[...,l] 
                                    + 15.0*(f[...,l+1]+f[...,l-1]) 
                                    -  6.0*(f[...,l+2]+f[...,l-2]) 
                                    +      (f[...,l+3]+f[...,l-3]))

    return d6fdx


def yder6_6th(f,dy,x=[],y=[],z=[]):
    if (f.ndim != 3 and f.ndim != 4):
        print("%s dimension arrays not handled." % (str(f.ndim)))
        raise ValueError

    m1 = 3
    m2 = f.shape[-2] - 3
    dy = N.gradient(y)
    if (dim.ny!=1):
        fac=1/dy**6
    d6fdy = N.zeros_like(f)
    
    dim=read_dim()

    if (m2 > m1 and dim.ny!=1):
        for m in range(m1,m2):
            d6fdy[...,m1:m2,:] = fac[m]*(- 20.0* f[...,m,:] 
                                         + 15.0*(f[...,m+1,:]+f[...,m-1,:]) 
                                         -  6.0*(f[...,m+2,:]+f[...,m-2,:]) 
                                         +      (f[...,m+3,:]+f[...,m-3,:]))
    
    return d6fdy

def zder6_6th(f,dz,x=[],y=[],z=[]):
    if (f.ndim != 3 and f.ndim != 4):
        print("%s dimension arrays not handled." % (str(f.ndim)))
        raise ValueError

    n1 = 3
    n2 = f.shape[-3] - 3
    dz = N.gradient(z)
    if (dim.nz!=1):
        fac=1/dz**6
    d6fdz = N.zeros_like(f)

    dim = read_dim()

    if (n2 > n1 and dim.nz!=1):
        for n in range(n1,n2):
            d6fdz[...,n,:,:] = fac[n]*(- 20.0* f[...,n,:,:] 
                                       + 15.0*(f[...,n+1,:,:]+f[...,n-1,:,:]) 
                                       -  6.0*(f[...,n+2,:,:]+f[...,n-2,:,:]) 
                                       +      (f[...,n+3,:,:]+f[...,n-3,:,:]))
    
    return d6fdz
    
