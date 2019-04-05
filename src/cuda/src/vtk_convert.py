#
# This python script can convert the datafiles into a vtk file to be used with Paraview. TODO
#

#Python modules
import numpy as np
import struct

#Own modules
import read_snapshot as rvar


def convert_arrays_to_vtk(xx, yy, zz, lnrho, uu_x, uu_y, uu_z, uu_tot, destination):
	# This subroutine will convert the array into a VTK file. 
	#
	# Set grid dimensions
	dimx = len(xx)
	dimy = len(yy)
	dimz = len(zz)
	dim = dimx * dimy * dimz
	dx = (np.max(xx) - np.min(xx))/(dimx-1)
	dy = (np.max(yy) - np.min(yy))/(dimy-1)
	dz = (np.max(zz) - np.min(zz))/(dimz-1)

	#Start writing the vtk file
	fd = open(destination + '.vtk', 'wb')
	fd.write('# vtk DataFile Version 2.0\n')
	fd.write('VAR files\n')
	fd.write('BINARY\n')
	fd.write('DATASET STRUCTURED_POINTS\n')
	fd.write('DIMENSIONS {0:9} {1:9} {2:9}\n'.format(dimx, dimy, dimz))
	fd.write('ORIGIN {0:8.12} {1:8.12} {2:8.12}\n'.format(xx[0], yy[0], zz[0]))
	fd.write('SPACING {0:8.12} {1:8.12} {2:8.12}\n'.format(dx, dy, dz))
	fd.write('POINT_DATA {0:9}\n'.format(dim))

	#Write the grid binaries
	print 'writing lnrho...'
	fd.write('SCALARS lnrho float\n')
	fd.write('LOOKUP_TABLE default\n')  
	#TODO: Check the ordering  
	for k in range(dimz):
		for j in range(dimy):
			for i in range(dimx):
				fd.write(struct.pack(">f", lnrho[i,j,k]))

        print 'writing uu_x, uu_y and uu_z...'
        fd.write('VECTORS vfield float\n')
        for k in range(dimz):
		for j in range(dimy):
			for i in range(dimx):
				fd.write(struct.pack(">f", uu_x[i,j,k]))
				fd.write(struct.pack(">f", uu_y[i,j,k]))
				fd.write(struct.pack(">f", uu_z[i,j,k]))

        print 'writing uu_tot...'
        fd.write('SCALARS uu_tot float\n')
        fd.write('LOOKUP_TABLE default\n')  
        for k in range(dimz):
                for j in range(dimy):
                        for i in range(dimx):
                                fd.write(struct.pack(">f", uu_tot[i,j,k]))

	fd.close()

########Parse arguments#########
import sys

#Set defaults
include_ghost=0

#Check user parametres

for idx in range (1, len(sys.argv)):
        if '--include_ghost' in sys.argv[idx]:
                include_ghost = int(sys.argv[idx].split("=")[1])

print "include_ghost = " + str(include_ghost)

#Read datafiles

xx0, yy0, zz0, lnrho0, uu_x0, uu_y0, uu_z0, uu_tot0 = rvar.read_whole_grid(snapshot=0, include_ghost=include_ghost)
xx, yy, zz, lnrho, uu_x, uu_y, uu_z, uu_tot = rvar.read_whole_grid(include_ghost=include_ghost)

#Write vtk file

convert_arrays_to_vtk(xx0, yy0, zz0, lnrho0, uu_x0, uu_y0, uu_z0, uu_tot0, "asth_grid0")

convert_arrays_to_vtk(xx, yy, zz, lnrho, uu_x, uu_y, uu_z, uu_tot0, "asth_grid")

