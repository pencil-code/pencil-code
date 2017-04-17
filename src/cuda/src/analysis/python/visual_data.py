import numpy as np
import pylab as plt
#from enthought.mayavi import mlab

from mpl_toolkits.mplot3d import Axes3D

import derivs as der 
import powerspectrum as pows
import read_snapshot as rvar
import plotting_tools as acplt

#def draw_isosurface(data, contour_number):
#	mlab.contour3d(data, opacity = 0.5, contours = contour_number) 


########Parse arguments#########
import sys

#Set defaults
view = "uu_z" 
slice_plane = "xy"
slice_pos = 64
curve_axis = 'x'
include_ghost=0

#Check user parametres

for idx in range (1, len(sys.argv)):
	if '--view' in sys.argv[idx]:
		view = sys.argv[idx].split("=")[1]
	elif '--slice_pos' in sys.argv[idx]:
		slice_pos = int(sys.argv[idx].split("=")[1])
	elif '--slice_plane' in sys.argv[idx]:
                slice_plane = sys.argv[idx].split("=")[1]
	elif '--curve_axis' in sys.argv[idx]:
                curve_axis = sys.argv[idx].split("=")[1]
	elif '--include_ghost' in sys.argv[idx]:
		include_ghost = int(sys.argv[idx].split("=")[1])

print "view = " + view
print "slice_plane = " + slice_plane
print "curve_axis = " + curve_axis
print "slice_pos = " + str(slice_pos)
print "include_ghost = " + str(include_ghost)

#Read datafiles

xx0, yy0, zz0, lnrho0, uu_x0, uu_y0, uu_z0, uu_tot0 = rvar.read_whole_grid(snapshot=0, include_ghost=include_ghost)
xx, yy, zz, lnrho, uu_x, uu_y, uu_z, uu_tot = rvar.read_whole_grid(include_ghost=include_ghost)

#Compare

diff_lnrho  = lnrho - lnrho0 
diff_uu_x = uu_x - uu_x0 
diff_uu_y = uu_y - uu_y0 
diff_uu_z = uu_z - uu_z0 
diff_uu_tot = uu_tot - uu_tot0 

#Inspect the arrays

rvar.check_inf_nan(lnrho0, uu_x0, uu_y0, uu_z0, uu_tot0, snapshot=0)
rvar.check_inf_nan(lnrho, uu_x, uu_y, uu_z, uu_tot)

#Plot data

if view == "uu_x":
	acplt.plot_simple(uu_x0, uu_x, slice_plane, slice_pos, curve_axis, view, xx, yy, zz)

if view == "uu_y":
#	draw_isosurface(diff_uu_y, 10)
	acplt.plot_simple(uu_y0, uu_y, slice_plane, slice_pos, curve_axis, view, xx, yy, zz)

if view == "uu_z":
#	draw_isosurface(diff_uu_z, 10)
	acplt.plot_simple(uu_z0, uu_z, slice_plane, slice_pos, curve_axis, view, xx, yy, zz)

if view == "lnrho":
	#draw_isosurface(lnrho, 10)
	acplt.plot_simple(lnrho0, lnrho, slice_plane, slice_pos, curve_axis, view, xx, yy, zz)

if view == "exp_lnrho":
#	draw_isosurface(np.exp(lnrho), 10)
	acplt.plot_simple(np.exp(lnrho0), np.exp(lnrho), slice_plane, slice_pos, curve_axis, view, xx, yy, zz)

if view == "exp_lnrho_isinf":
	print np.nonzero(np.isinf(np.exp(lnrho)))
	fig = plt.figure()
	ax = Axes3D(fig) 
	X, Y, Z = np.nonzero(np.isinf(np.exp(lnrho)))
	ax.scatter(X, Y, Z)
	plt.show()

if view == "uu_tot":
#	draw_isosurface(uu_tot, 3)
	acplt.plot_simple(uu_tot0, uu_tot, slice_plane, slice_pos, curve_axis, view, xx, yy, zz)

if view == "uu_x_lnrho":
        acplt.plot_simple(uu_x, lnrho, slice_plane, slice_pos, curve_axis, view, xx, yy, zz)

if view == "uu_y_lnrho":
        acplt.plot_simple(uu_y, lnrho, slice_plane, slice_pos, curve_axis, view, xx, yy, zz)	

if view == "uu_z_lnrho":
        acplt.plot_simple(uu_z, lnrho, slice_plane, slice_pos, curve_axis, view, xx, yy, zz)

if view == "helicity":
	oox, ooy, ooz = der.curl(uu_x, uu_y, uu_z, dx, dy, dz)
	helical_grid = der.dot(uu_x, uu_y, uu_z, oox, ooy, ooz)
        oo_tot = np.sqrt( np.power(oox,2.0) + np.power(ooy,2.0) + np.power(ooz,2.0) )
	#draw_isosurface(helical_grid, 3)
	plt.figure() 
	plt.imshow(helical_grid[:,slice_pos,:])
	plt.figure() 
	plt.imshow(oox[:,slice_pos,:])
	plt.figure() 
	plt.imshow(ooy[:,slice_pos,:])
	plt.figure() 
	plt.imshow(ooz[:,slice_pos,:])
	plt.figure() 
	plt.imshow(oo_tot[:,slice_pos,:])
'''
if view == "powersp_uutot":
	vkx, vky, vkz, kx, ky, kz = pows.powerspectrum(uu_tot, NX, NY, NZ, Lbox)
	plt.figure() 
	plt.plot(kx, vkx)
	plt.figure() 
	plt.plot(ky, vky)
	plt.figure() 
	plt.plot(kz, vkz)

if view == "powersp_uu_x":
	vkx, vky, vkz, kx, ky, kz = pows.powerspectrum(uu_x, NX, NY, NZ, Lbox)
	plt.figure() 
	plt.plot(kx, vkx)
	plt.figure() 
	plt.plot(ky, vky)
	plt.figure() 
	plt.plot(kz, vkz)

if view == "powersp_uu_y":
	vkx, vky, vkz, kx, ky, kz = pows.powerspectrum(uu_y, NX, NY, NZ, Lbox)
	plt.figure() 
	plt.plot(kx, vkx)
	plt.figure() 
	plt.plot(ky, vky)
	plt.figure() 
	plt.plot(kz, vkz)

if view == "powersp_uu_z":
	vkx, vky, vkz, kx, ky, kz = pows.powerspectrum(uu_z, NX, NY, NZ, Lbox)
	plt.figure() 
	plt.plot(kx, vkx)
	plt.figure() 
	plt.plot(ky, vky)
	plt.figure() 
	plt.plot(kz, vkz)
'''
if view == "powersp_all":
	#plt.rc('text', usetex=True)
	#plt.rc('font', family='serif')
	#print "uu_tot..."
	#vkx_tot, vky_tot, vkz_tot, kx_tot, ky_tot, kz_tot = pows.powerspectrum(uu_tot, NX, NY, NZ, Lbox)
	Lbox = COMP_DOMAIN_SIZE_X*DX #Size of the computational box edge.  
	print "uu_x..."
	vkx_x, vky_x, vkz_x, kx_x, ky_x, kz_x = pows.powerspectrum(uu_x, COMP_DOMAIN_SIZE_X, COMP_DOMAIN_SIZE_Y, COMP_DOMAIN_SIZE_Z, Lbox)
	print "uu_y..."
	vkx_y, vky_y, vkz_y, kx_y, ky_y, kz_y = pows.powerspectrum(uu_y, COMP_DOMAIN_SIZE_X, COMP_DOMAIN_SIZE_Y, COMP_DOMAIN_SIZE_Z, Lbox)
	print "uu_z..."
	vkx_z, vky_z, vkz_z, kx_z, ky_z, kz_z = pows.powerspectrum(uu_z, COMP_DOMAIN_SIZE_X, COMP_DOMAIN_SIZE_Y, COMP_DOMAIN_SIZE_Z, Lbox)
	print "All powerspectra taken"

	plt.figure() 
	plt.title("$k_x$")
	plt.xlabel("$k$")
	#plt.plot(kx_tot, vkx_tot, "k-", kx_x, vkx_x, "r-", kx_y, vkx_y, "g-", kx_z, vkx_z, "b-")
	#plt.legend( (r'$| \mathbf{u} |$', r'$u_x$', r'$u_y$', r'$u_z$') )
	plt.plot(kx_x, vkx_x, "r-", kx_y, vkx_y, "g-", kx_z, vkx_z, "b-")
	plt.legend( ('$u_x$', '$u_y$', '$u_z$') )
        plt.axis([0, 20, 0, np.max([vkx_x, vkx_y, vkx_z])])
	plt.figure() 
	plt.title("$k_y$")
	plt.xlabel("$k$")
	#plt.plot(ky_tot, vky_tot, "k-", ky_x, vky_x, "r-", ky_y, vky_y, "g-", ky_z, vky_z, "b-")
	#plt.legend( (r'$| \mathbf{u} |$', r'$u_x$', r'$u_y$', r'$u_z$') )
	plt.plot(ky_x, vky_x, "r-", ky_y, vky_y, "g-", ky_z, vky_z, "b-")
	plt.legend( ('$u_x$', '$u_y$', '$u_z$') )
        plt.axis([0, 20, 0, np.max([vky_x, vky_y, vky_z])])
	plt.figure()
	plt.title("$k_z$") 
	plt.xlabel("$k$")
	#plt.plot(kz_tot, vkz_tot, "k-", kz_x, vkz_x, "r-", kz_y, vkz_y, "g-", kz_z, vkz_z, "b-")
	#plt.legend( (r'$| \mathbf{u} |$', r'$u_x$', r'$u_y$', r'$u_z$') )
	plt.plot(kz_x, vkz_x, "r-", kz_y, vkz_y, "g-", kz_z, vkz_z, "b-")
	plt.legend( ('$u_x$', '$u_y$', '$u_z$') )
        plt.axis([0, 20, 0, np.max([vkz_x, vkz_y, vkz_z])])

        plt.figure()
        plt.title("$k_x$ combined")
        plt.xlabel("$k$")
        kxx = np.sqrt(vkx_x**2.0 + vkx_y**2.0 + vkx_z**2.0)
        plt.plot(kx_x, kxx, "r-")
        plt.axis([0, 20, 0, np.max(kxx)])
        plt.figure()

        plt.title("$k_y$ combined")
        plt.xlabel("$k$")
        kyy = np.sqrt(vky_x**2.0 + vky_y**2.0 + vky_z**2.0)
        plt.plot(ky_x, kyy, "r-")
        plt.axis([0, 20, 0, np.max(kyy)])

        plt.figure()
        plt.title("$k_z$ combined")
        plt.xlabel("$k$")
	kzz = np.sqrt(vkz_x**2.0 + vkz_y**2.0 + vkz_z**2.0)
        plt.plot(kz_x, kzz, "r-")
        plt.axis([0, 20, 0, np.max(kzz)])


#Plot all figures

plt.show()

#Null grid data

lnrho0 = []
lnrho  = []
uu_x0    = []
uu_x     = []
uu_y0    = []
uu_y     = []
uu_z0    = []
uu_z     = []

diff_lnrho  = [] 
diff_uu_x     = [] 
diff_uu_y     = [] 
diff_uu_z     = [] 

