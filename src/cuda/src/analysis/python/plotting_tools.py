import numpy as np
import pylab as plt

def plot_simple(var0, var, slice_plane, slice_pos, curve_axis, title, xx, yy, zz):
	#Simple plotting tool to look to compare var files
	print "Plotting... " + title
	title0 = title+"0"

	#plt.rc('text', usetex=True)
	#plt.rc('font', family='serif')

	if slice_plane == "xy":
		plt.figure()
		plt.subplot(121)
		plt.imshow(var0[:,:,slice_pos])
		plt.xlabel('y')
		plt.ylabel('x')
		plt.title(title0)
		plt.subplot(122)
		plt.imshow(var[:,:,slice_pos])
		plt.xlabel('y')
		plt.ylabel('x')
		plt.title(title)
	elif slice_plane == "xz":
		plt.figure()
		plt.subplot(121)
		plt.imshow(var0[:,slice_pos,:])
		plt.xlabel('z')
		plt.ylabel('x')
		plt.title(title0)
		plt.subplot(122)
		plt.imshow(var[:,slice_pos,:])
		plt.xlabel('z')
		plt.ylabel('x')
		plt.title(title)
	elif slice_plane == "yz":
		plt.figure()
		plt.subplot(121)
		plt.imshow(var0[slice_pos,:,:])
		plt.xlabel('z')
		plt.ylabel('y')
		plt.title(title0)
		plt.subplot(122)
		plt.imshow(var[slice_pos,:,:])
		plt.xlabel('z')
		plt.ylabel('y')
		plt.title(title)
	elif slice_plane == "all":
		plt.figure()
                plt.subplot(131)
                plt.imshow(var[:,:,slice_pos])
                plt.xlabel('y')
                plt.ylabel('x')
                plt.title(title)

                plt.subplot(132)
                plt.imshow(var[:,slice_pos,:])
                plt.xlabel('z')
                plt.ylabel('x')
                plt.title(title)

		plt.subplot(133)
		plt.imshow(var[slice_pos,:,:])
		plt.xlabel('z')
		plt.ylabel('y')
                plt.title(title)
		

	#
	# Transform this so that all directions can be plotted here for the comparison. 
	#
	if curve_axis == "x":
		plt.figure()
		plt.plot(xx, var0[:,slice_pos,slice_pos], xx, var[:,slice_pos,slice_pos])
		plt.xlabel('x')
		plt.ylabel(title)
		plt.title(title)
#		plt.figure()
#		plt.plot(xx, diff_var[:,slice_pos,slice_pos])
        elif curve_axis == "y":
		plt.figure()
		plt.plot(yy, var0[slice_pos,:,slice_pos], yy, var[slice_pos,:,slice_pos])
		plt.xlabel('y')
		plt.ylabel(title)
		plt.title(title)
	elif curve_axis == "z":
		plt.figure()
		plt.plot(zz, var0[slice_pos,slice_pos,:], zz, var[slice_pos,slice_pos,:])
		plt.xlabel('z')
		plt.ylabel(title)
		plt.title(title)
	elif curve_axis == "xyz":
		plt.figure()
		varx0 = var0[:,slice_pos,slice_pos]
		vary0 = var0[slice_pos,:,slice_pos]
		varz0 = var0[slice_pos,slice_pos,:]
		varx = var[:,slice_pos,slice_pos]
		vary = var[slice_pos,:,slice_pos]
		varz = var[slice_pos,slice_pos,:]
		#plt.plot(xx, varx0, 'r-', yy, vary0, 'r--', zz, varz0, 'r-.', xx, varx, 'b-', yy, vary, 'b--', zz, varz, 'b-.')
		plt.plot(xx, varx0, 'r-', label="x-axis 0")
		plt.plot(yy, vary0, 'r--', label="y-axis 0") 
		plt.plot(zz, varz0, 'r-.', label="z-axis 0") 
		plt.plot(xx, varx, 'b-', label="x-axis")
		plt.plot(yy, vary, 'b--', label="y-axis") 
		plt.plot(zz, varz, 'b-.', label="z-axis")
		plt.xlabel('x/y/z')
		plt.ylabel(title)
		plt.title(title)
		plt.legend()

def plot_once(var, slice_plane, slice_pos, title):
	#Simple plotting tool to look to compare var files
	print "Plotting... " + title

	if slice_plane == "xy":
		plt.figure()
		plt.imshow(var[:,:,slice_pos])
		plt.xlabel('y')
		plt.ylabel('x')
		plt.title(title)
	elif slice_plane == "xz":
		plt.figure()
		plt.imshow(var[:,slice_pos,:])
		plt.xlabel('z')
		plt.ylabel('x')
		plt.title(title)
	elif slice_plane == "yz":
		plt.figure()
		plt.imshow(var[slice_pos,:,:])
		plt.xlabel('z')
		plt.ylabel('y')
		plt.title(title)
	elif slice_plane == "all":
		plt.figure()
                plt.subplot(131)
                plt.imshow(var[:,:,slice_pos])
                plt.xlabel('y')
                plt.ylabel('x')
                plt.title(title)

                plt.subplot(132)
                plt.imshow(var[:,slice_pos,:])
                plt.xlabel('z')
                plt.ylabel('x')
                plt.title(title)

		plt.subplot(133)
		plt.imshow(var[slice_pos,:,:])
		plt.xlabel('z')
		plt.ylabel('y')
                plt.title(title)
		

