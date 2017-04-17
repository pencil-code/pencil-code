import numpy as np

import derivs as der 

def len_tsfile(datafile):
	#Calculates the length of the time series file
	tsinf = open(datafile, "r")

	line = tsinf.readline() # Skip the title line

	nstep = 0
	for line in tsinf:
		nstep = nstep + 1 

	tsinf.close()

	nstep = nstep - 1 #remove the extra step

	print nstep

	return nstep

def draw_isosurface(data, contour_number):
	mlab.contour3d(data, opacity = 0.5, contours = contour_number) 

def read_nn_inf(datafile):
	#READ nx ny nz
	nninf = open(datafile, "r")

	line = nninf.readline()
	line = line.split()

	nx = int(line[0])
	ny = int(line[1])
	nz = int(line[2])

	nninf.close()

	print " nx =", nx, " ny =", ny, " nz =", nz

	return nx, ny, nz

def read_ts_inf(datafile):
	#READ nstep
	#READ step t dt urms umax umin uxmin uxmax uymin uymax uzmin uzmax Marms Mamax
	nstep = len_tsfile(datafile)

	tsinf = open(datafile, "r")

	line = tsinf.readline()
	line = line.split()

	headers = line

	steps = np.zeros((nstep+1), dtype=int)
	var_list = np.zeros((nstep+1, len(headers)-1))

	kk = 0
	for line in tsinf:
		line = line.split()
		steps[kk] = int(line[0])
		#print steps[kk]
		for ii in np.arange(len(headers)-1):
			var_list[kk, ii] = float(line[ii+1])
		kk = kk+1 

	tsinf.close()

	return var_list, steps, headers 

def read_grid_dat(filename, nx, ny, nz):
	#Read and reconstruct the binary data

	datafile = open(filename, "rb")

	dtype = np.dtype(np.single) #C-type floating point number 
	all_numbers = np.fromfile(datafile, dtype) #Read the binary file

	datafile.close()

	#Copy data into cube form
	
	data_cube = np.zeros((nx, ny, nz))

	for i in range(0,nx): 
		for j in range(0,ny): 
			for k in range(0,nz): 
				ind = i + j*nx + k*nx*ny
				data_cube[i,j,k] = all_numbers[ind]
				#print all_numbers[ind]
	
	#Remove ghost zones

	grid = np.zeros((nx-6, ny-6, nz-6))

	xtop = nx-3
	ytop = ny-3
	ztop = nz-3
	grid = data_cube[3:xtop,3:ytop,3:ztop]

	#Clear useless arrays to free memory
	all_numbers = []
	data_cube = []


	return grid


#Read data

nx, ny, nz = read_nn_inf("data/nn.inf")

xx = np.arange((nx-6))
yy = np.arange((ny-6))
zz = np.arange((nz-6))

unit_lenght = 1.0
unit_density = 1.0
unit_velocity = 1.0

Lbox = 2.0
dx = Lbox / (nx - 6)
dy = dx
dz = dx
print "dx = ", dx

lnrho0 = read_grid_dat("data/density0.dat", nx, ny, nz)
lnrho  = read_grid_dat("data/density.dat", nx, ny, nz)
velx0    = read_grid_dat("data/velx0.dat", nx, ny, nz)
velx     = read_grid_dat("data/velx.dat", nx, ny, nz)
vely0    = read_grid_dat("data/vely0.dat", nx, ny, nz)
vely     = read_grid_dat("data/vely.dat", nx, ny, nz)
velz0    = read_grid_dat("data/velz0.dat", nx, ny, nz)
velz     = read_grid_dat("data/velz.dat", nx, ny, nz)

vel_tot = np.sqrt( np.power(velx,2.0) + np.power(vely,2.0) + np.power(velz,2.0) )
vel_tot0 = np.sqrt( np.power(velx0,2.0) + np.power(vely0,2.0) + np.power(velz0,2.0) )

rho0 = np.exp(lnrho0)
rho = np.exp(lnrho)

#Import arrays in yt

from yt.mods import *
from yt.frontends.stream.api import load_uniform_grid

data = dict(Density = rho, Velocity_x = velx, Velocity_y = vely, Velocity_z = velz, Velocity = vel_tot)
bbox = np.array([[-Lbox/2.0, Lbox/2.0], [-Lbox/2.0, Lbox/2.0], [-Lbox/2.0, Lbox/2.0]])
pf = load_uniform_grid(data, rho.shape, unit_lenght, bbox=bbox)

#pc = PlotCollection(pf, [0.0, 0.0, 0.0])
'''
#Save slice plots 

SlicePlot(pf, 'x', 'Density').save()
SlicePlot(pf, 'y', 'Density').save()
SlicePlot(pf, 'z', 'Density').save()
SlicePlot(pf, 'x', 'Velocity').save()
SlicePlot(pf, 'y', 'Velocity').save()
SlicePlot(pf, 'z', 'Velocity').save()
SlicePlot(pf, 'x', 'Velocity_x').save()
SlicePlot(pf, 'y', 'Velocity_x').save()
SlicePlot(pf, 'z', 'Velocity_x').save()
SlicePlot(pf, 'x', 'Velocity_y').save()
SlicePlot(pf, 'y', 'Velocity_y').save()
SlicePlot(pf, 'z', 'Velocity_y').save()
SlicePlot(pf, 'x', 'Velocity_z').save()
SlicePlot(pf, 'y', 'Velocity_z').save()
SlicePlot(pf, 'z', 'Velocity_z').save()


#Save projection plots

ProjectionPlot(pf, 'x', 'Density').save()
ProjectionPlot(pf, 'y', 'Density').save()
ProjectionPlot(pf, 'z', 'Density').save()
ProjectionPlot(pf, 'x', 'Velocity').save()
ProjectionPlot(pf, 'y', 'Velocity').save()
ProjectionPlot(pf, 'z', 'Velocity').save()
ProjectionPlot(pf, 'x', 'Velocity_x').save()
ProjectionPlot(pf, 'y', 'Velocity_x').save()
ProjectionPlot(pf, 'z', 'Velocity_x').save()
ProjectionPlot(pf, 'x', 'Velocity_y').save()
ProjectionPlot(pf, 'y', 'Velocity_y').save()
ProjectionPlot(pf, 'z', 'Velocity_y').save()
ProjectionPlot(pf, 'x', 'Velocity_z').save()
ProjectionPlot(pf, 'y', 'Velocity_z').save()
ProjectionPlot(pf, 'z', 'Velocity_z').save()
'''
#
#Save isosurface images
#

field = 'Velocity'

# Find the bounds for your field
dd = pf.h.all_data()
mi, ma = dd.quantities["Extrema"](field)[0]

# Instantiate the ColorTransferfunction.
tf = ColorTransferFunction((mi, ma))

# Set up the camera parameters: center, looking direction, width, resolution
c = (pf.domain_right_edge + pf.domain_left_edge)/2.0
L = np.array([1.0, 1.0, 1.0])
W = 1.0
N = 512

# Create a camera object
cam = pf.h.camera(c, L, W, N, tf, fields = [field])

# Now let's add some isocontours, and take a snapshot, saving the image
# to a file.
tf.add_layers(10, colormap = 'RdBu_r')
im = cam.snapshot('test_rendering.png')

# To add the domain box to the image:
nim = cam.draw_domain(im)
nim.write_png('test_rendering_with_domain.png')

#
#Null grid data
#

density0 = []
density  = []
velx0    = []
velx     = []
vely0    = []
vely     = []
velz0    = []
velz     = []

diff_density  = [] 
diff_velx     = [] 
diff_vely     = [] 
diff_velz     = [] 



