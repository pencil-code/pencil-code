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
import sys

def particles_to_density(xxp,yyp,zzp,x,y,z):
    dim = read_dim()
    if dim.precision == 'D':
        precision = 'd'
    else:
        precision = 'f'
    nnp = np.zeros([dim.mz,dim.my,dim.mx])
                
    npar=len(xxp)
    
    for k in range(npar):

        xp=xxp[k]
        yp=yyp[k]
        zp=zzp[k]

        ix0=find_index_bisect(xp,x)
        iy0=find_index_bisect(yp,y)
        iz0=find_index_bisect(zp,z)

        ixx0=ix0-1
        ixx1=ix0+1

        iyy0=iy0-1
        iyy1=iy0+1

        izz0=iz0-1
        izz1=iz0+1

        if (dim.nx > 1):
            dx=x[ix0]-x[ixx0]
        else:
            dx=1.0

        if (dim.ny > 1):
            dy=y[iy0]-y[iyy0]
        else:
            dy=1.0

        if (dim.nz > 1):
            dz=z[iz0]-z[izz0]
        else:
            dz=1.0

        dx_1=1.0/dx
        dy_1=1.0/dy
        dz_1=1.0/dz
#
        dx_2=1.0/dx**2
        dy_2=1.0/dy**2
        dz_2=1.0/dz**2
    
#if (not ghost) then begin
#  mx=nx+6 & l1=3 & l2=l1+nx-1
#  x2=fltarr(mx)*one
#  x2[l1:l2]=x
#  for l=l1-1,   0,-1 do x2[l]=x2[l+1]-dx
#  for l=l2+1,mx-1,+1 do x2[l]=x2[l-1]+dx
#  x=x2
                                                                                         
#  my=ny+6 & m1=3 & m2=m1+ny-1
#  y2=fltarr(my)*one
#  y2[m1:m2]=y
#  for m=m1-1,   0,-1 do y2[m]=y2[m+1]-dy
#  for m=m2+1,my-1,+1 do y2[m]=y2[m-1]+dy
#  y=y2
#;                                                                                         
#    if (dim.nz == 1): 
#        for n in np.arange(3):
#            nn=n+1
#            z[dim.n1-nn]=z[dim.n1]-nn*dz
#            z[dim.n2+nn]=z[dim.n2]+nn*dz

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

                    weight=1.0
                    if (dim.nx != 1): weight=weight*weight_x
                    if (dim.ny != 1): weight=weight*weight_y
                    if (dim.nz != 1): weight=weight*weight_z
#                                                                                                                                                                                                         
                    nnp[izz,iyy,ixx]=nnp[izz,iyy,ixx] + weight

    return nnp
#---------------------------------------------------------------
def find_index_bisect(qpar,q):

    jl=0
    ju=len(q)-1

    while ((ju-jl)>1):
        jm=(ju+jl)/2
        if (qpar > q[jm]):
            jl=jm
        else:
            ju=jm

    if (qpar-q[jl] <= q[ju]-qpar):
        iq0=jl
    else:
        iq0=ju

    return iq0
