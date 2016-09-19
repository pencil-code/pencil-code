##
##  $Id$
##
##  Convert positions of particles to a grid density field.
##
##  Author: Wladimir Lyra (adapted from Anders Johansen's IDL script) 
##
import numpy as np
import os.path
#Pencil routines
from pencil.files.dim import read_dim
from pencil.files.grid import read_grid
import sys

def particles_to_density(xp,yp,zp,x,y,z):
    dim = read_dim()
    if dim.precision == 'D':
        precision = 'd'
    else:
        precision = 'f'

    npar=len(xp)
    #g=read_grid()

    if (dim.nx > 1): 
        dx=x[1]-x[0]
    else:
        dx=1.0
        
    if (dim.ny > 1):
        dy=y[1]-y[0]
    else:
        dy=1.0
        
    if (dim.nz > 1):
        dz=z[1]-z[0]
    else:
        dz=1.0

    dx_1=1.0/dx
    dy_1=1.0/dy
    dz_1=1.0/dz
    
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
    if (dim.nz == 1): 
        for n in np.arange(3):
            nn=n+1
            z[dim.n1-nn]=z[dim.n1]-nn*dz
            z[dim.n2+nn]=z[dim.n2]+nn*dz


    nnp=np.zeros((dim.mx,dim.my,dim.mz)) #(dim.mx,dim.my,dim.mz))

    for k in np.arange(npar): 

        ix = np.round((xp[k]-x[0])*dx_1)
        iy = np.round((yp[k]-y[0])*dy_1)
        iz = np.round((zp[k]-z[0])*dz_1)

        if (ix == dim.l2+1): ix=ix-1
        if (iy == dim.m2+1): iy=iy-1
        if (iz == dim.n2+1): iz=iz-1
        if (ix == dim.l1-1): ix=ix+1
        if (iy == dim.m1-1): iy=iy+1
        if (iz == dim.n1-1): iz=iz+1

##  Particles are assigned to the nearest grid point.

        nnp[ix,iy,iz]=nnp[ix,iy,iz]+1.0

    #print 'nnp minmax=',nnp.min(), nnp.max()
    nnp=nnp[dim.l1:dim.l2+1,dim.m1:dim.m2+1,dim.n1:dim.n2+1]

    if (dim.nz == 1): nnp=nnp[:,:,0]
##  Convert from count to density.
#
    return nnp
