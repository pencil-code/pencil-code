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
from pencil.files.param import read_param
from pencil.files.grid import read_grid

def xder_6th(f,dx,x=[],y=[],z=[]):

    if (f.ndim != 3 and f.ndim != 4):
        print("%s dimension arrays not handled." % (str(f.ndim)))
        raise ValueError

    dx=N.gradient(x)
    dx2 = 1./(60.*dx)
    dfdx = N.zeros_like(f)
    l1 = 3
    l2 = f.shape[-1]-3
    if (l2 > l1):
        for l in range(l1,l2):
            dfdx[...,l] = dx2[l]*( +45.*(f[...,l+1]-f[...,l-1])
                                     -9.*(f[...,l+2]-f[...,l-2])
                                     +(f[...,l+3]-f[...,l-3]) )
    else:
        dfdx = 0.
    return dfdx

def yder_6th(f,dy,x=[],y=[],z=[]):

    if (f.ndim != 3 and f.ndim != 4):
        print("%s dimension arrays not handled." % (str(f.ndim)))
        raise ValueError

    param=read_param(quiet=True)

    dy=N.gradient(y)
    dy2 = 1./(60.*dy)
    dfdy = N.zeros_like(f)
    m1 = 3
    m2 = f.shape[-2]-3

    if (m2 > m1):
        for m in range(m1,m2):
            dfdy[...,m,:] = dy2[m]*( +45.*(f[...,m+1,:]-f[...,m-1,:]) 
                                      -9.*(f[...,m+2,:]-f[...,m-2,:]) 
                                      +(f[...,m+3,:]-f[...,m-3,:]) )
    else:
        dfdy = 0.
    if param.coord_system == ('cylindric' or 'spherical'):
        if len(x) < 1:
            gd=read_grid(quiet=True)
            x=gd.x
        dfdy /= x
        
    return dfdy

def zder_6th(f,dz,x=[],y=[],z=[],run2D=False):
    
    if (f.ndim != 3 and f.ndim != 4):
        print("%s dimension arrays not handled." % (str(f.ndim)))
        raise ValueError

    param=read_param(quiet=True)

    dz=N.gradient(z)
    dz2 = 1./(60.*dz)
    dfdz = N.zeros_like(f)
    n1 = 3
    if run2D:
        n2 = f.shape[1]-3
    else:
        n2 = f.shape[-3]-3


    if (n2 > n1):
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

def xder2_6th(f,dx,x=[],y=[],z=[]):

    if (f.ndim != 3 and f.ndim != 4):
        print("%s dimension arrays not handled." % (str(f.ndim)))
        raise ValueError

    l1 = 3
    l2 = f.shape[-1]-3
    dx = N.gradient(x[l1:l2])
    dx2 = 1./(180.*dx**2.)
    dfdx = N.zeros_like(f)

    if (l2 > l1):
        dfdx[...,l1:l2] = dx2*(-490.*f[...,l1:l2]
                              +270.*(f[...,l1-1:l2-1]+f[...,l1+1:l2+1])
                              - 27.*(f[...,l1-2:l2-2]+f[...,l1+2:l2+2])
                              +  2.*(f[...,l1-3:l2-3]+f[...,l1+3:l2+3]) )
    else:
        dfdx = 0.

    return dfdx

def yder2_6th(f,dy,x=[],y=[],z=[]):

    if (f.ndim != 3 and f.ndim != 4):
        print("%s dimension arrays not handled." % (str(f.ndim)))
        raise ValueError
    
    param=read_param(quiet=True)

    m1 = 3
    m2 = f.shape[-2]-3
    dy = N.gradient(y) #[m1:m2])
    dy2 = 1./(180.*dy**2.)
    dfdy = N.zeros_like(f)

    if (m2 > m1):
        for m in range(m1,m2+1):
            dfdy[...,m,:] = dy2[m]*(-490.*f[...,m,:]
                                     +270.*(f[...,m-1,:]+f[...,m+1,:])
                                     - 27.*(f[...,m-2,:]+f[...,m+2,:])
                                     +  2.*(f[...,m-3,:]+f[...,m+3,:]) )
        #dfdy[...,m1:m2,:] = dy2*(-490.*f[...,m1:m2,:]
        #                      +270.*(f[...,m1-1:m2-1,:]+f[...,m1+1:m2+1,:])
        #                      - 27.*(f[...,m1-2:m2-2,:]+f[...,m1+2:m2+2,:])
        #                      +  2.*(f[...,m1-3:m2-3,:]+f[...,m1+3:m2+3,:]) )
    else:
        dfdy = 0.
    if param.coord_system == ('cylindric' or 'spherical'):
        if (len(x) or len(y)) < 1:
            gd=read_grid(quiet=True)
            x=gd.x; y=gd.y
        dfdy /= x**2

    return dfdy

def zder2_6th(f,dz,x=[],y=[],z=[]):

    if (f.ndim != 3 and f.ndim != 4):
        print("%s dimension arrays not handled." % (str(f.ndim)))
        raise ValueError

    param=read_param(quiet=True)

    n1 = 3
    n2 = f.shape[-3]-3
    dz = N.gradient(z[n1:n2])
    dz2 = 1./(180.*dz**2.)
    dfdz = N.zeros_like(f)

    if (n2 > n1):
        dfdz[...,n1:n2,:,:] = dz2*(-490.*f[...,n1:n2,:,:]
                              +270.*(f[...,n1-1:n2-1,:,:]+f[...,n1+1:n2+1,:,:])
                              - 27.*(f[...,n1-2:n2-2,:,:]+f[...,n1+2:n2+2,:,:])
                              +  2.*(f[...,n1-3:n2-3,:,:]+f[...,n1+3:n2+3,:,:]) )
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
    dx = N.gradient(x[l1:l2])
    fac=1/dx**6
    d6fdx = N.zeros_like(f)

    if (l2 > l1):
        d6fdx[...,l1:l2] = fac*(- 20.0* f[...,l1:l2] 
                     + 15.0*(f[...,l1+1:l2+1]+f[...,l1-1:l2-1]) 
                     -  6.0*(f[...,l1+2:l2+2]+f[...,l1-2:l2-2]) 
                     +      (f[...,l1+3:l2+3]+f[...,l1-3:l2-3]))

    return d6fdx


def yder6_6th(f,dy,x=[],y=[],z=[]):
    if (f.ndim != 3 and f.ndim != 4):
        print("%s dimension arrays not handled." % (str(f.ndim)))
        raise ValueError

    m1 = 3
    m2 = f.shape[-2] - 3

    dy = N.gradient(y[m1:m2])
    fac=1/dy**6
    d6fdy = N.zeros_like(f)

    if (m2 > m1):
        d6fdy[...,m1:m2,:] = fac*(- 20.0* f[...,m1:m2,:] 
                                 + 15.0*(f[...,m1+1:m2+1,:]+f[...,m1-1:m2-1,:]) 
                                 -  6.0*(f[...,m1+2:m2+2,:]+f[...,m1-2:m2-2,:]) 
                                 +      (f[...,m1+3:m2+3,:]+f[...,m1-3:m2-3,:]))
    
    return d6fdy

def zder6_6th(f,dz,x=[],y=[],z=[]):
    if (f.ndim != 3 and f.ndim != 4):
        print("%s dimension arrays not handled." % (str(f.ndim)))
        raise ValueError

    n1 = 3
    n2 = f.shape[-3] - 3
    dz = N.gradient(z[n1:n2])
    fac=1/dz**6

    d6fdy = N.zeros_like(f)

    if (n2 > n1):
        d6fdy[...,n1:n2,:,:] = fac*(- 20.0* f[...,n1:n2,:,:] 
                       + 15.0*(f[...,n1+1:n2+1,:,:]+f[...,n1-1:n2-1,:,:]) 
                       -  6.0*(f[...,n1+2:n2+2,:,:]+f[...,n1-2:n2-2,:,:]) 
                       +      (f[...,n1+3:n2+3,:,:]+f[...,n1-3:n2-3,:,:]))
    
    return d6fdy
    
