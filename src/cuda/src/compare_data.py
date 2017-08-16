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
dir1="src/data/"
dir2="src/data/"
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
	elif '--dir1' in sys.argv[idx]:
		dir1 = sys.argv[idx].split("=")[1]
	elif '--dir2' in sys.argv[idx]:
		dir2 = sys.argv[idx].split("=")[1]
	elif '--include_ghost' in sys.argv[idx]:
		slice_pos = int(sys.argv[idx].split("=")[1])


print "view = " + view
print "slice_plane = " + slice_plane
print "curve_axis = " + curve_axis
print "slice_pos = " + str(slice_pos)
print "dir1 = " + dir1
print "dir2 = " + dir2
print "include_ghost = " + str(include_ghost)


#Read datafiles

xx1, yy1, zz1, lnrho1, uu_x1, uu_y1, uu_z1, uu_tot1 = rvar.read_whole_grid(data_dir=dir1, include_ghost=include_ghost)
xx2, yy2, zz2, lnrho2, uu_x2, uu_y2, uu_z2, uu_tot2 = rvar.read_whole_grid(data_dir=dir2, include_ghost=include_ghost)

#Compare

diff_lnrho  = lnrho2 - lnrho1 
diff_uu_x = uu_x2 - uu_x1 
diff_uu_y = uu_y2 - uu_y1 
diff_uu_z = uu_z2 - uu_z1 
diff_uu_tot = uu_tot2 - uu_tot1 

#Plot data

if view == "uu_x":
	acplt.plot_simple(uu_x1, uu_x2, slice_plane, slice_pos, curve_axis, view, xx2, yy2, zz2)
	acplt.plot_once(diff_uu_x, slice_plane, slice_pos, "uu_x2 - uu_x1")

if view == "uu_y":
	acplt.plot_simple(uu_y1, uu_y2, slice_plane, slice_pos, curve_axis, view, xx2, yy2, zz2)
	acplt.plot_once(diff_uu_y, slice_plane, slice_pos, "uu_y2 - uu_y1")

if view == "uu_z":
	acplt.plot_simple(uu_z1, uu_z2, slice_plane, slice_pos, curve_axis, view, xx2, yy2, zz2)
	acplt.plot_once(diff_uu_z, slice_plane, slice_pos, "uu_z2 - uu_z1")

if view == "lnrho":
	acplt.plot_simple(lnrho1, lnrho2, slice_plane, slice_pos, curve_axis, view, xx2, yy2, zz2)
	acplt.plot_once(diff_lnrho, slice_plane, slice_pos, "lnrho2 - lnrho1")

if view == "uu_tot":
	acplt.plot_simple(uu_tot1, uu_tot2, slice_plane, slice_pos, curve_axis, view, xx2, yy2, zz2)
	acplt.plot_once(diff_uu_tot, slice_plane, slice_pos, "uu_tot2 - uu_tot1")

#Plot all figures

plt.show()

#Null grid data

lnrho1  = []
lnrho2  = []
uu_x1   = []
uu_x2   = []
uu_y1   = []
uu_y2   = []
uu_z1   = []
uu_z2   = []
uu_tot1 = [] 
uu_tot2 = []

diff_lnrho  = [] 
diff_uu_x   = [] 
diff_uu_y   = [] 
diff_uu_z   = [] 
diff_uu_tot = []

