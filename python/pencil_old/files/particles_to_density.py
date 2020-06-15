##
##  $Id$
##
##  Convert positions of particles to a grid density field.
##
##  Author: Wladimir Lyra (adapted from Anders Johansen's IDL script) 
##
import numpy as np
import os.path
from .dim import read_dim
from .param import read_param
import sys

def particles_to_density(xxp,yyp,zzp,x,y,z,density=True):
    dim = read_dim()

    if dim.precision == 'D':
        precision = 'd'
    else:
        precision = 'f'
    nnp = np.zeros([dim.mz,dim.my,dim.mx])

    npar=len(xxp)

    par = read_param(quiet=True)

    dx = np.gradient(x)
    dy = np.gradient(y)
    dz = np.gradient(z)
#
    dx1=1.0/dx
    dy1=1.0/dy
    dz1=1.0/dz
#
    dx2=1.0/dx**2
    dy2=1.0/dy**2
    dz2=1.0/dz**2
#
    if (par.grid_func[0]=='linear'): dx1_pt = dx1[0]
    if (par.grid_func[1]=='linear'): dy1_pt = dy1[0]
    if (par.grid_func[2]=='linear'): dz1_pt = dz1[0]
#
    for k in range(npar):

        xp=xxp[k]
        yp=yyp[k]
        zp=zzp[k]

        if (par.grid_func[0]=='linear'):
            ix0 = int(round((xp-x[0])*dx1_pt))
        else:
            ix0=find_index_bisect(xp,x)
        if (ix0 == dim.l2+1): ix0=ix0-1
        if (ix0 == dim.l1-1): ix0=ix0+1
        dx_1=dx1[ix0]
        dx_2=dx2[ix0]
#
        if (par.grid_func[1]=='linear'):
            iy0 = int(round((yp-y[0])*dy1_pt))
        else:
            iy0=find_index_bisect(yp,y)
        if (iy0 == dim.m2+1): iy0=iy0-1
        if (iy0 == dim.m1-1): iy0=iy0+1
        dy_1=dy1[iy0]
        dy_2=dy2[iy0]
#
        if (par.grid_func[2]=='linear'):
            iz0 = int(round((zp-z[0])*dz1_pt))
        else:
            iz0=find_index_bisect(zp,z)
        if (iz0 == dim.n2+1): iz0=iz0-1
        if (iz0 == dim.n1-1): iz0=iz0+1
        dz_1=dz1[iz0]
        dz_2=dz2[iz0]
#
        ixx0 = ix0-1
        ixx1 = ix0+1
#
        iyy0 = iy0-1
        iyy1 = iy0+1
#
        izz0 = iz0-1
        izz1 = iz0+1
#
        for ixx in np.arange(ixx0,ixx1+1):
            for iyy in np.arange(iyy0,iyy1+1):
                for izz in np.arange(izz0,izz1+1):

                    if ( ((ixx-ix0) == -1) or ((ixx-ix0) == +1) ):
                        weight_x = 1.125 - 1.5*abs(xp-x[ixx])  *dx_1 + 0.5*abs(xp-x[ixx])**2*dx_2
                    else:
                        if (dim.nx != 1): weight_x = 0.75 - (xp-x[ixx])**2*dx_2

                    if ( ((iyy-iy0) == -1) or ((iyy-iy0) == +1) ):
                        weight_y = 1.125 - 1.5*abs(yp-y[iyy])  *dy_1 + 0.5*abs(yp-y[iyy])**2*dy_2
                    else:
                        if (dim.ny != 1): weight_y = 0.75 - (yp-y[iyy])**2*dy_2

                    if ( ((izz-iz0) == -1) or ((izz-iz0) == +1) ):
                        weight_z = 1.125 - 1.5*abs(zp-z[izz])  *dz_1 + 0.5*abs(zp-z[izz])**2*dz_2
                    else:
                        if (dim.nz != 1): weight_z = 0.75 - (zp-z[izz])**2*dz_2

                    if density then:
                        if type(par.rhop_swarm) is float:
                            weight=par.rhop_swarm
                        else:
                            weight=par.rhop_swarm[k]
                    else:
                        weight=1.0


                    if (dim.nx != 1): weight=weight*weight_x
                    if (dim.ny != 1): weight=weight*weight_y
                    if (dim.nz != 1): weight=weight*weight_z
#
                    nnp[izz,iyy,ixx]=nnp[izz,iyy,ixx] + weight

    return nnp
#​
# ---------------------------------------------------------------
#
def find_index_bisect(qpar, q):
#
    jl=0
    ju=len(q)-1
#
    while ((ju-jl)>1):
        jm=(ju+jl)//2
        if (qpar > q[jm]):
            jl=jm
        else:
            ju=jm

    if (qpar-q[jl] <= q[ju]-qpar):
        iq0=jl
    else:
        iq0=ju

    return iq0
