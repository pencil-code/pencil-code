import numpy as np
import pylab as plt
import scipy


def read_slices_dat(filebase, nx, ny, nz, nslice):
	#Read and reconstruct the binary data

	slice_cube = np.zeros((nx, ny, nslice))

	for sli in range(0,nslice):
		try:

			filename = filebase + str(sli) + ".ani"

			print filename

			datafile = open(filename, "rb")

			dtype = np.dtype(np.single) #C-type floating point number 
			all_numbers = np.fromfile(datafile, dtype) #Read the binary file

			datafile.close()

			#Copy data slice cube form
	
			for i in range(0,nx): 
				for j in range(0,nz): 
					ind = i + j*nx
					slice_cube[i,j,sli] = all_numbers[ind]
		except IOError:
			print "No " + filename + " found!"

	return slice_cube

def save_figure(X, Z, aslice, levels, xx, yy, zz, title, filename):
        if np.size(levels) > 0:	
		fig = plt.figure()
		##ax = fig.add_subplot(111, aspect="equal")
		cont = plt.contourf(X, Z, aslice, levels, cmap=plt.get_cmap("hot"))
		cont.cmap.set_under('white')
		cont.cmap.set_over('black')
	        plt.xlim( np.min(xx), np.max(zz) )
	        plt.ylim( np.min(zz), np.max(zz) )
	        #plt.imshow(aslice)
		plt.colorbar()
	        plt.title(title)
	        plt.savefig(filename)
		#plt.show()
	        plt.close()
	else:
		print filename + " NOT SAVED -> NO CONTOUR LEVELS"

def save_figure_onedim(axis, aslice, title, filename, slicemin, slicemax):
        plt.figure()
        plt.plot(axis, aslice, "r-")
        maxloc = np.argmax(aslice)
        minloc = np.argmin(aslice)
        plt.plot(axis[maxloc], aslice[maxloc],"b^")
        plt.plot(axis[minloc], aslice[minloc],"k*")
        plt.title(title)
        plt.ylim( slicemin, slicemax )
        plt.savefig(filename)
        plt.close()

def save_as_images_pylab(slice_cube, steps, xx, yy, zz, filebase, title, onedim=0, compare=0):
	
	slice_max = np.nanmax(slice_cube[:,:,steps/2])
	slice_min = np.nanmin(slice_cube[:,:,steps/2])
        print slice_max, slice_min
	
        #Determine contour levels
        nlevs = 100.0
        spacing = (slice_max- slice_min)/nlevs
        print str(slice_min) + " " + str(slice_max) + " " + str(spacing)
	#NOTE; spacing == 0.0 for lnrho with initial values (2.0 for all points)
	if (spacing == 0.0): 
		spacing = 0.01
        levels = np.arange(slice_min, slice_max, spacing)

        X, Z = np.meshgrid(xx, zz)

        if onedim:
                onedimtxt = "onedim_"
        else:
                onedimtxt = ""

        if onedim:
		if compare:
                        slicemin = np.nanmin(slice_cube[:, zz.size/2, :]) 
                        slicemax = np.nanmax(slice_cube[:, zz.size/2, :])
		else:
			slicemin = np.nanmin(slice_cube[:, zz.size/2, :])
                	slicemax = np.nanmax(slice_cube[:, zz.size/2, :])
	
	for sli in range(0,steps):
		number = "0000"
		number = number[:len(number) - len(str(sli))] + str(sli)
        	filename = filebase + onedimtxt + number + ".png"

		print filename

                if onedim:
                        zz = np.array(zz)
                        if compare:
				aslice = (slice_cube[:, zz.size/2, sli] 
					  - slice_cube[:, zz.size/2, 0])
			else:
                		aslice = slice_cube[:, zz.size/2, sli]

                	save_figure_onedim(xx, aslice, title, filename, slicemin, slicemax)
                else:
               		aslice = slice_cube[:,:,sli]
                	save_figure(X, Z, aslice, levels, xx, yy, zz, title, filename)


	



########Parse arguments#########
import sys

nslice = 390 #Default values for arguments
onedim = 0
compare = 0

print "BEFORE nslice = " + str(nslice)

print sys.argv

for idx in range (1, len(sys.argv)):
	if '--nslice' in sys.argv[idx]:
		nslice = int(sys.argv[idx].split("=")[1])
		print "LOOP nslice = " + str(nslice)
        elif '--onedim' in sys.argv[idx]:
                onedim = int(sys.argv[idx].split("=")[1])
        elif '--compare' in sys.argv[idx]:
                compare = int(sys.argv[idx].split("=")[1])

print "nslice = " + str(nslice)

print "onedim=", onedim

print "compare = ", compare 

#Rad grid info
grid_file = open("data/grid_info.ac", "r")

line = grid_file.readline()
line = line.split()

nx = int(line[0])
ny = int(line[1])
nz = int(line[2])
grid_size = int(line[3])
comp_dom_size_x = int(line[4])
comp_dom_size_y = int(line[5])
comp_dom_size_z = int(line[6])
pad_size = int(line[7])
bound_size = int(line[8])
dx = float(line[9])
dy = float(line[10])
dz = float(line[11])
forcing = float(line[12])
grid_file.close()

nx = comp_dom_size_x +6#include bounds, purkka
ny = comp_dom_size_y +6
nz = comp_dom_size_z +6
#Commented away the minus sixes. Messes with the data reading of slice binaries 
#and produces the effect we see. 
xx = np.arange(nx)
yy = np.arange(ny)
zz = np.arange(nz)

#XY-plane

animation = read_slices_dat("data/animation/lnrho_z_", nx, ny, nz, nslice)
save_as_images_pylab(animation, nslice, xx, yy, zz, "animation_images/lnrho_z_", 
                     'LNRHO, sliced z axis', onedim=onedim, compare=compare)

animation = read_slices_dat("data/animation/uu_x_z_", nx, ny, nz, nslice)
save_as_images_pylab(animation, nslice, xx, yy, zz, "animation_images/uu_x_z_", 
                     'UU_X, sliced z axis', onedim=onedim, compare=compare)

animation = read_slices_dat("data/animation/uu_y_z_", nx, ny, nz, nslice)
save_as_images_pylab(animation, nslice, xx, yy, zz, "animation_images/uu_y_z_", 
                     'UU_Y, sliced z axis', onedim=onedim, compare=compare)

animation = read_slices_dat("data/animation/uu_z_z_", nx, ny, nz, nslice)
save_as_images_pylab(animation, nslice, xx, yy, zz, "animation_images/uu_z_z_", 
                     'UU_Z, sliced z axis', onedim=onedim, compare=compare)

animation = read_slices_dat("data/animation/uu_z_", nx, ny, nz, nslice)
save_as_images_pylab(animation, nslice, xx, yy, zz, "animation_images/uu_z_", 
                     'UU total, sliced z axis', onedim=onedim, compare=compare)


#compare_images(animation2, 20)














































