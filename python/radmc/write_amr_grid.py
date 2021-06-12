import pencil as pc
import numpy as np
import math
from astropy.constants import au
#
datadir = '../data'
unit_length=au.to('cm').value
par=pc.read.param(datadir=datadir)
#
# According to the RADMC3D manual, the way the variable coordsystem works is
#
# If coordsystem<100 the coordinate system is cartesian.
# If 100<=coordsystem<200 thecoordinate system is spherical (polar).
# If 200<=coordsystem<300 the coordinate system is cylindrical.
#
xdim=unit_length
print(par.coord_system)
if (par.coord_system == 'cartesian'):
    coordsystem = 99
    ydim=unit_length
    zdim=unit_length
elif (par.coord_system == 'cylindric'):
    coordsystem = 200
    ydim=1.
    zdim=unit_length
elif (par.coord_system == 'spherical'):
    coordsystem = 100
    ydim=1.
    zdim=1.
else:
    print("the world is flat and we never got here")
    #break
            
grid=pc.read.grid(trim=True,datadir=datadir)
dim=pc.read.dim(datadir=datadir)

iformat=1
grid_style=0
gridinfo=0
if (dim.nx > 1):
    incl_x=1
    dx=np.gradient(grid.x)
else:
    incl_x=0
    dx=np.repeat(grid.dx,dim.nx)
if (dim.ny > 1):
    incl_y=1
    dy=np.gradient(grid.y)
else:
    incl_y=0
    dy=np.repeat(grid.dy,dim.ny)
if (dim.nz > 1):
    incl_z=1
    dz=np.gradient(grid.z)
else:
    incl_z=0
    dz=np.repeat(grid.dz,dim.nz)
    
f = open('amr_grid.inp', 'w')

f.write(str(iformat)+'\n')
f.write(str(grid_style)+'\n')
f.write(str(coordsystem)+'\n')
f.write(str(gridinfo)+'\n')
f.write(str(incl_x)+' '+str(incl_y)+' '+str(incl_z)+'\n')
f.write(str(dim.nx)+' '+str(dim.ny)+' '+str(dim.nz)+'\n')
for i in range(dim.nx):
    f.write(str('%.9e' % ((grid.x[i]-dx[i]/2)*xdim))+'\n')
f.write(str('%.9e' % ((grid.x[dim.nx-1]+dx[dim.nx-1]/2)*xdim))+'\n')    
#            
for m in range(dim.ny):
    f.write(str('%.9e' % ((grid.y[m]-dy[m]/2)*ydim))+'\n')
f.write(str('%.9e' % ((grid.y[dim.ny-1]+dy[dim.ny-1]/2)*ydim))+'\n')    
#    
for n in range(dim.nz):
    f.write(str('%.9e' % ((grid.z[n]-dz[n]/2)*zdim))+'\n')
f.write(str('%.9e' % ((grid.z[dim.nz-1]+dz[dim.nz-1]/2)*zdim))+'\n')

f.close()


