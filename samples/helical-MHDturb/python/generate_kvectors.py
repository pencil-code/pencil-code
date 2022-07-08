#!/usr/bin/env python3
"""
generate_kvectors.py
Author: Kishore G (kishore96@gmail.com)

Generate table of wave numbers in a given range for helical forcing.
Adapted from samples/helical-MHDturb/idl/generate_kvectors.pro

Anisotropy is achieved by stretching the sphere to an ellipsoid,
so k1 < kref < k2, where kref^2=kx^2+ky^2+kz^2/ez^2.
So, for an 1:2 anisotropy (kz=2*kx) we put ez=2.

For tall boxes, it is useful to allow all smaller kz vectors, so we
put, for a 16x16x256 box, for example, ;dkz=1./16. instead of 1.
"""

# BEGIN configuration
dkx = 1
dky = 1
dkz = 1
ex = 1
ey = 1
ez = 1
kmax = 10
k1 = 4.5
k2 = 5.5
kmaxz = kmax
debug = False #Whether to print debug output
# END configuration

import numpy as np

kx_list = []
ky_list = []
kz_list = []

kav = 0 #Used to keep track of the average of the generated k vectors
for kx in np.arange(-kmax, kmax, dkx):
	for ky in np.arange(-kmax, kmax, dky):
		for kz in np.arange(-kmax, kmax, dkz):
			k = np.sqrt(kx**2 + ky**2 + kz**2)
			kref = np.sqrt( (kx/ex)**2 + (ky/ey)**2 + (kz/ez)**2 )
			if kref > k1 and kref < k2:
				kav += kref
				kx_list.append(kx)
				ky_list.append(ky)
				kz_list.append(kz)

kav = kav/len(kx_list)

with open("k.dat", 'w') as f:
	"""
	Documention for list-directed IO (which forcing.f90 uses to read k.dat):
	https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnc5/index.html
	
	Apparently tabs cause problems in some Fortran compilers, so we have to use spaces in the output.
	"""
	f.write("    {}    {}".format(len(kx_list), kav) )
	f.write("\n")
	
	def print_list(f, l):
		nrec = 6 #number of records per line
		for i in range(0,len(l)):
			if i!= 0 and i%nrec == 0:
				f.write('\n')
			f.write("    {}".format(l[i]))
		f.write('\n')
	
	print_list(f, kx_list)
	print_list(f, ky_list)
	print_list(f, kz_list)

if debug:
	import matplotlib.pyplot as plt
	
	kkx = np.array(kx_list)
	kky = np.array(ky_list)
	kkz = np.array(kz_list)
	
	print("<k> = ", np.average(kkx), np.average(kky), np.average(kkz))
	print("<k^2> = ", np.average(kkx**2), np.average(kky**2), np.average(kkz**2))
	
	kkxm = np.where(kky<0.4, kkx, np.nan)
	kkzm = np.where(kky<0.4, kkz, np.nan)
	fig,ax = plt.subplots()
	ax.scatter(kkxm, kkzm)
	plt.show()
