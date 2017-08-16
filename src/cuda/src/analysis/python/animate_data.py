import numpy as np
import pylab as plt
import scipy
#import Image

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

def save_as_images(slice_cube, steps, filebase):
	
	slice_max = np.nanmax(slice_cube[:,:,steps/2])
	slice_min = np.nanmin(slice_cube[:,:,steps/2])

	for sli in range(0,steps):
		if sli >= 0:
			number = "00000" + str(sli)
		if sli >= 10:
  			number = "0000"  + str(sli)
		if sli >= 100:
  			number = "000"   + str(sli)
		if sli >= 1000:
  			number = "00"    + str(sli)
		if sli >= 10000:
  			number = "0"     + str(sli)
		if sli >= 100000:
  			number = ""      + str(sli)

        	filename = filebase + number + ".png"

		print filename
		
		aslice = slice_cube[:,:,sli]

		# Scale array into image range 
		aslice_pix = ((aslice+abs(slice_min))/(slice_max+abs(slice_min))) * 255.0
		aslice_pix = np.uint8(aslice_pix) #TARKSITA SKAALAUS

		im = scipy.misc.toimage(aslice_pix)
		im.save(filename)
	#im.show()

        print slice_max, slice_min

def save_figure(X, Z, aslice, levels, xx, yy, zz, title, filename):
        if np.size(levels) > 0:
		plt.figure()
		plt.contourf(X, Z, aslice, levels)
	        plt.xlim( np.min(xx), np.max(zz) )
	        plt.ylim( np.min(zz), np.max(zz) )
	        #plt.imshow(aslice)
		plt.colorbar()
	        plt.title(title)
	        plt.savefig(filename)
	        plt.close()
	else:
		print filename + " NOT SAVED -> NO CONTOUR LEVELS"

def save_figure_onedim(axis, aslice, title, filename):
        plt.figure()
        plt.plot(axis, aslice, "r-")
        plt.title(title)
        plt.savefig(filename)
        plt.close()

def save_as_images_pylab(slice_cube, steps, xx, yy, zz, filebase, title, onedim=0):
	
	slice_max = np.nanmax(slice_cube[:,:,steps/2])
	slice_min = np.nanmin(slice_cube[:,:,steps/2])
        print slice_max, slice_min

        #Determine contour levels
        nlevs = 100.0
        spacing = (slice_max- slice_min)/nlevs
        print str(slice_min) + " " + str(slice_max) + " " + str(spacing)
        levels = np.arange(slice_min, slice_max, spacing)

        X, Z = np.meshgrid(xx, zz)

        if onedim:
                onedimtxt = "_onedim_"
        else:
		onedimtxt = ""

	for sli in range(0,steps):
		if sli >= 0:
			number = "00000" + str(sli)
		if sli >= 10:
  			number = "0000"  + str(sli)
		if sli >= 100:
  			number = "000"   + str(sli)
		if sli >= 1000:
  			number = "00"    + str(sli)
		if sli >= 10000:
  			number = "0"     + str(sli)
		if sli >= 100000:
  			number = ""      + str(sli)

        	filename = filebase + onedimtxt + number + ".png"

		print filename
		
                if onedim:
                   aslice = slice_cube[:, zz.size[0]/2, sli]

                   save_figure_onedim(xx, aslice, title, filename) 
                else:
  		   aslice = slice_cube[:,:,sli]

		   save_figure(X, Z, aslice, levels, xx, yy, zz, title, filename)




def compare_images(slice_cube, step):
	
	slice_max = np.nanmax(slice_cube)
		
	aslice = slice_cube[:,:,step]

	# Scale array into image range 

	aslice_pix = (aslice/slice_max) * 255.0

	aslice_pix_int = np.uint8(aslice_pix)				
	im = scipy.misc.toimage(aslice_pix_int)
	im.show()

	plt.figure()
	plt.imshow(aslice)
	plt.figure()
	plt.imshow(aslice_pix)
	plt.figure()
	plt.imshow(aslice_pix_int)
	plt.show()

	print aslice_pix_int[50,:]
	print aslice_pix_int[:,50]


########Parse arguments#########
import sys

nslice = 10 #Default values for arguments
onedim = 0

for idx in range (1, len(sys.argv)):
	if '--nslice' in sys.argv[idx]:
		nslice = int(sys.argv[idx].split("=")[1])
        elif '--onedim' in sys.argv[idx]:
                onedim = int(sys.argv[idx].split("=")[1])
print "nslice = " + str(nslice)

print "onedim=", onedim

nx, ny, nz = read_nn_inf("data/nn.inf")

#Commented away the minus sixes. Messes with the data reading of slice binaries 
#and produces the effect we see. 

nx = nx # -6 
ny = ny # -6 
nz = nz # -6 

xx = np.arange(nx)
yy = np.arange(ny)
zz = np.arange(nz)

#XY-plane

animation1 = read_slices_dat("data/animation/density_xy_", nx, ny, nz, nslice)
animation2 = read_slices_dat("data/animation/uu_xy_", nx, ny, nz, nslice)
animation3 = read_slices_dat("data/animation/velx_xy_", nx, ny, nz, nslice)
animation4 = read_slices_dat("data/animation/vely_xy_", nx, ny, nz, nslice)
animation5 = read_slices_dat("data/animation/velz_xy_", nx, ny, nz, nslice)

save_as_images_pylab(animation1, nslice, xx, yy, zz, "animation_images/density_xy_", 'Density', onedim=onedim)
save_as_images_pylab(animation2, nslice, xx, yy, zz, "animation_images/uu_xy_", 'Velocity', onedim=onedim)
save_as_images_pylab(animation3, nslice, xx, yy, zz, "animation_images/velx_xy_", 'Velocity (x)', onedim=onedim)
save_as_images_pylab(animation4, nslice, xx, yy, zz, "animation_images/vely_xy_", 'Velocity (y)', onedim=onedim)
save_as_images_pylab(animation5, nslice, xx, yy, zz, "animation_images/velz_xy_", 'Velocity (z)', onedim=onedim)

#XZ-plane

animation1 = read_slices_dat("data/animation/density_xz_", nx, ny, nz, nslice)
animation2 = read_slices_dat("data/animation/uu_xz_", nx, ny, nz, nslice)
animation3 = read_slices_dat("data/animation/velx_xz_", nx, ny, nz, nslice)
animation4 = read_slices_dat("data/animation/vely_xz_", nx, ny, nz, nslice)
animation5 = read_slices_dat("data/animation/velz_xz_", nx, ny, nz, nslice)

#save_as_images(animation1, nslice, "animation_images/density_im_xz_")
#save_as_images(animation2, nslice, "animation_images/uu_im_xz_")
#save_as_images(animation3, nslice, "animation_images/velx_im_xz_")
#save_as_images(animation4, nslice, "animation_images/vely_im_xz_")
#save_as_images(animation5, nslice, "animation_images/velz_im_xz_")

save_as_images_pylab(animation1, nslice, xx, yy, zz, "animation_images/density_xz_", 'Density', onedim=onedim)
save_as_images_pylab(animation2, nslice, xx, yy, zz, "animation_images/uu_xz_", 'Velocity', onedim=onedim)
save_as_images_pylab(animation3, nslice, xx, yy, zz, "animation_images/velx_xz_", 'Velocity (x)', onedim=onedim)
save_as_images_pylab(animation4, nslice, xx, yy, zz, "animation_images/vely_xz_", 'Velocity (y)', onedim=onedim)
save_as_images_pylab(animation5, nslice, xx, yy, zz, "animation_images/velz_xz_", 'Velocity (z)', onedim=onedim)

#YZ-plane

animation1 = read_slices_dat("data/animation/density_yz_", nx, ny, nz, nslice)
animation2 = read_slices_dat("data/animation/uu_yz_", nx, ny, nz, nslice)
animation3 = read_slices_dat("data/animation/velx_yz_", nx, ny, nz, nslice)
animation4 = read_slices_dat("data/animation/vely_yz_", nx, ny, nz, nslice)
animation5 = read_slices_dat("data/animation/velz_yz_", nx, ny, nz, nslice)

save_as_images_pylab(animation1, nslice, xx, yy, zz, "animation_images/density_yz_", 'Density', onedim=onedim)
save_as_images_pylab(animation2, nslice, xx, yy, zz, "animation_images/uu_yz_", 'Velocity', onedim=onedim)
save_as_images_pylab(animation3, nslice, xx, yy, zz, "animation_images/velx_yz_", 'Velocity (x)', onedim=onedim)
save_as_images_pylab(animation4, nslice, xx, yy, zz, "animation_images/vely_yz_", 'Velocity (y)', onedim=onedim)
save_as_images_pylab(animation5, nslice, xx, yy, zz, "animation_images/velz_yz_", 'Velocity (z)', onedim=onedim)

#compare_images(animation2, 20)

