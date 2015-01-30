#!/usr/bin/env python

import pencil as pc
import numpy as np
import pylab as pyl
import mathgl as mgl
import math as math

from npfile import npfile
import os
import sys
from mpl_toolkits.mplot3d import axes3d as p3
import mayavi.mlab as mlab

testarray = np.array([np.arange(27)])

np.set_printoptions(threshold='nan')

class pcpvar:
		"""
		a class to hold the pvar data
		"""
		def __init__(self,header,npar_loc,ipar,npar_data,footer):
			self.ipar=ipar
			self.npar_data=npar_data
			self.header = header
			self.npar_loc = npar_loc
#			self.stuff = stuff
			self.footer = footer
			
			
def read_npar_loc(casedir='.',datadir='/data',pfile='pvar.dat',proc=0):
	array_npar_loc=np.dtype([('header','<i4'),
							('npar_loc','<i4')])
	npars = np.fromfile(casedir+datadir+'/proc'+str(proc)+'/'+pfile,dtype=array_npar_loc)
	return npars['npar_loc'][0]

def read_npvar(casedir='.',datadir='/data',pfile='pvar.dat',proc=0):
	
	dims = pc.read_dim(datadir=casedir+datadir,proc=proc)
	pdims = pc.read_pdim(casedir+datadir)
	npar_loc = read_npar_loc(casedir,datadir,pfile,proc)
#	print npar_loc,' particles on processor: ',proc
	mvars = pdims.mpaux+pdims.mpvar
	ltot = npar_loc*mvars

	array_shape= np.dtype([('header','<i4'),
							('npar_loc','<i4'),
							('footer','<i4'),
							('header2','<i4'),
							('ipar','<i4',npar_loc),
							('footer2','<i4'),
							('header3','<i4'),
							('fp','<f8',ltot),
							('footer3','<i4'),
							('header4','<i4'),
							('t','<f8'),
							('x','<f8',dims.mx),
							('y','<f8',dims.my),
							('z','<f8',dims.mz),
							('dx','<f8'),
							('dy','<f8'),
							('dz','<f8'),
							('footer4','<i4')])
	
	p_data = np.fromfile(casedir+datadir+'/proc'+str(proc)+'/'+pfile,dtype=array_shape)
	partpars = p_data['fp'].reshape(mvars,npar_loc)
	ipar = np.squeeze(p_data['ipar'])
	part_pos = np.array([partpars[0,:],partpars[1,:],partpars[2,:]])
	xdims = np.squeeze(p_data['x'])
	ydims = np.squeeze(p_data['y'])
	zdims = np.squeeze(p_data['z'])
	uboundx = xdims[-dims.nghostx]-0.5*p_data['dx']
	lboundx = xdims[dims.nghostx]-0.5*p_data['dx']
	uboundy = ydims[-dims.nghosty]-0.5*p_data['dy']
	lboundy = ydims[dims.nghosty]-0.5*p_data['dy']
	uboundz = zdims[-dims.nghostz]-0.5*p_data['dz']
	lboundz = zdims[dims.nghostz]-0.5*p_data['dz']
	procdims = np.squeeze([(uboundx, uboundy, uboundz), (lboundx, lboundy, lboundz)])
	header = p_data['header']
	dl = np.squeeze([p_data['dx'],p_data['dy'],p_data['dz']])
	footer = p_data['footer']
	npar_loc = np.squeeze(p_data['npar_loc'])
	return ipar,part_pos,procdims,dl
	
	
def calc_alpha(datadir='./data',pfile='pvar.dat',proc=0):
	
	params =pc.read_param(datadir,quiet=True)
	tdims = pc.read_pdim(datadir)
	cdim = pc.read_dim(datadir)
	ipar,part_pos,procdims,dl = read_npvar(datadir,pfile,proc)
	np_init = tdims.npar / float(cdim.nx*cdim.ny*cdim.nz)
#	print 'number of particles per grid cell: ',np_init
	dims = pc.read_dim(datadir,proc)
	blocklength = 2.0
	lx = procdims[0,0]-procdims[1,0]	
	ly = procdims[0,1]-procdims[1,1]
	lz = procdims[0,2]-procdims[1,2]
	dx = dl[0] * blocklength
	dy = dl[1] * blocklength
	dz = dl[2] * blocklength
	nx = int(lx /dx)
	ny = int(ly /dy)
	nz = int(lz /dz)
#	print 'one block is '+ str(blocklength**3)+' cells big and should contains on average ' +str(np_init * blocklength**3)+' particles'
	sx = procdims[1,0]
	sy = procdims[1,1]
	sz = procdims[1,2]
#	print 'lx,ly,lz: ',lx,ly,lz
#	print 'dx,dy,dz: ',dx,dy,dz
#	print 'start_x,start_y,start_z: ',sx,sy,sz
	cell_array = np.zeros(nx*ny*nz)
	
	for i in np.arange(len(ipar)):
#		print 'particle '+str(ipar[i])+' is at pos: '+str(part_pos[0,i])+' '+str(part_pos[1,i])+' '+str(part_pos[2,i])
		for j in np.arange(nx-1):
			for k in np.arange(ny-1):
				for l in np.arange(nz-1):
					if (sx + dx*j < part_pos[0,i] <= sx + dx*(j+1) and \
						sy + dy*k < part_pos[1,i] <= sy + dy*(k+1) and \
						sz + dz*l < part_pos[2,i] <= sz + dz*(l+1)):
#						print 'located particle '+ str(ipar[i]) + ' in cell ' + str(l*(n_cells**2) + (k*n_cells)+j) 
						cell_array[l*(nz**2) + (k*ny)+j] += np_init/float(blocklength**3)
	return cell_array
	

def plot_parts(casedir='.',datadir='/data',pfile='pvar.dat',proc=0):

	ipar,part_pos,procdims,dl = read_npvar(datadir=datadir)
	dims = pc.read_dim(datadir=casedir+datadir,proc=proc)
	lx = procdims[0,0]-procdims[1,0]	
	ly = procdims[0,1]-procdims[1,1]
	sx = procdims[1,0]
	sy = procdims[1,1]
	sz = procdims[1,2]
	x,y,z=list(part_pos[0,:]),list(part_pos[1,:]),list(part_pos[2,:])
	s=range(len(x))
	mlab.points3d(x,y,z,s,mode='point')
	mlab.show()
#	fig = pyl.figure(1)
#	ax = fig.gca(projection='3d')
#	plot = ax.plot(part_pos[0,:],part_pos[1,:],part_pos[2,:],'.')

def plot_parts_full(datadir='/data'):
	
	dummy,part_pos = collect_pdata(datadir=datadir)
	x,y,z=list(part_pos[0,:]),list(part_pos[1,:]),list(part_pos[2,:])
	s=range(len(x))
	mlab.points3d(x,y,z,s,mode='point')
	mlab.show()
	#fig = pyl.figure(1)
	#ax = fig.gca(projection='3d')
	#plot = ax.plot(part_pos[0,:],part_pos[1,:],part_pos[2,:],'.')
#	mlab.points3d(part_pos[0,:],part_pos[1,:],part_pos[2,:])
	

def plot_hist(cell_array):
	n, bins, patches = pyl.hist(cell_array, bins=max(cell_array), cumulative=False,normed=1, histtype='bar')

def plot_isosurf(datadir='/data',blocklength=2,surfaces=[1.0],smooth=False,smoothing=5):
	cell_array,nx,ny,nz = calc_tot_alpha(datadir=datadir,blocklength=blocklength)
#	for i in np.arange(len(cell_array)):
#		if cell_array[i]<.2:
#			cell_array[i] = 1.0
#		else:
#			cell_array[i] = 0.0
#		cell_array[i] = math.exp(-cell_array[i])
	if smooth:
		print 'smoothing...'
		scalar = trid_smooth(cell_array,nx,ny,nz,smoothing)
	scalar = np.reshape(cell_array,(nz,ny,nx))
	src = mlab.pipeline.scalar_field(scalar)
	mlab.pipeline.iso_surface(src, contours=surfaces, opacity=0.7)
	#mlab.show()
#	mlab.pipeline.volume(src)
	mlab.pipeline.image_plane_widget(src,
                            plane_orientation='z_axes',
                            slice_index=10,
                        )
	mlab.show()


def collect_pdata(casedir='.',datadir='/data'):
	procs = np.arange(32)
	for i in procs:
		dom_ipar,dom_part_pos = read_npvar_red(casedir=casedir,datadir=datadir,proc=i)
		if i == 0:
			ipars = dom_ipar
			part_pos = dom_part_pos
		else:
			ipars = np.hstack((ipars,dom_ipar))
			part_pos = np.hstack((part_pos,dom_part_pos))

	return ipars,part_pos

def calc_tot_alpha(casedir='.',datadir='/data',blocklength=4):
	ipars,part_pos = collect_pdata(casedir,datadir)
	params =pc.read_param(casedir+datadir,quiet=True)
	cdim = pc.read_dim(casedir+datadir)
	pdims = pc.read_pdim(casedir+datadir)
	lx = params.lx0
	ly = params.ly0
	lz = params.lz0
	nx = cdim.nx/blocklength
	ny = cdim.ny/blocklength
	nz = cdim.nz/blocklength
	dx = lx /float(nx)
	dy = ly /float(ny)
	dz = lz /float(nz)
	sx = params.xyz0[0]
	sy = params.xyz0[1]
	sz = params.xyz0[2]
	cell_array=np.zeros(nx*ny*nz)
	np_init=pdims.npar / float(nx*ny*nz)
#	print 'average particle number in '+ str(blocklength**3)+ ' cells is: '+str(np_init)
	for i in np.arange(len(ipars)):
#		if (i*100)%(10*len(ipars))==0:
#			print str(i/float(len(ipars))*100) + '% of particles have been located'
		j = divmod(part_pos[0,i],dx)[0]
		k = divmod(part_pos[1,i],dy)[0]
		l = divmod(part_pos[2,i],dz)[0]

		cell_array[l*(nz**2) + (k*ny)+j] += np_init/float(blocklength**3)
	return cell_array,nx,ny,nz


def read_npvar_red(casedir='.',datadir='/data',pfile='pvar.dat',proc=0):
	
	dims = pc.read_dim(casedir+datadir,proc)
	pdims = pc.read_pdim(casedir+datadir)
	npar_loc = read_npar_loc(casedir=casedir,datadir=datadir,pfile=pfile,proc=proc)
#	print npar_loc,' particles on processor: ',proc
	mvars = pdims.mpaux+pdims.mpvar
	ltot = npar_loc*mvars

	array_shape= np.dtype([('header','<i4'),
							('npar_loc','<i4'),
							('footer','<i4'),
							('header2','<i4'),
							('ipar','<i4',npar_loc),
							('footer2','<i4'),
							('header3','<i4'),
							('fp','<f8',ltot),
							('footer3','<i4'),
							('header4','<i4'),
							('t','<f8'),
							('x','<f8',dims.mx),
							('y','<f8',dims.my),
							('z','<f8',dims.mz),
							('dx','<f8'),
							('dy','<f8'),
							('dz','<f8'),
							('footer4','<i4')])
	
	p_data = np.fromfile(casedir+datadir+'/proc'+str(proc)+'/'+pfile,dtype=array_shape)
	partpars = p_data['fp'].reshape(mvars,npar_loc)
	ipar = np.squeeze(p_data['ipar'])
	part_pos = np.array([partpars[0,:],partpars[1,:],partpars[2,:]])
	return ipar,part_pos
	
def smoothTriangle(data,degree,dropVals=False):
        """performs moving triangle smoothing with a variable degree."""
        """note that if dropVals is False, output length will be identical
        to input length, but with copies of data at the flanking regions"""
        triangle=np.array(range(degree)+[degree]+range(degree)[::-1])+1
        smoothed=[]
        for i in range(degree,len(data)-degree*2):
                point=data[i:i+len(triangle)]*triangle
                smoothed.append(sum(point)/sum(triangle))
        if dropVals: return smoothed
        smoothed=[smoothed[0]]*(degree+degree/2)+smoothed
        while len(smoothed)<len(data):smoothed.append(smoothed[-1])
        return smoothed

def trid_smooth(cell_array,nx,ny,nz,smoothing=5):
#  smoothing along the x-axis
	arrayx = smoothTriangle(cell_array,smoothing,dropVals=False)
#  assembling the scalar field and swapping x with y axis
	arrayy = np.reshape(cell_array,(nx,ny,nz))
	arrayy = np.swapaxes(arrayy,0,1)
	arrayy = np.reshape(arrayy,(nx*ny*nz))
#  smoothing along y axis
	arrayy = smoothTriangle(arrayy,smoothing,dropVals=False)
	arrayy = np.reshape(arrayy,(ny,nx,nz))
	arrayy = np.swapaxes(arrayy,0,1)
	arrayy = np.reshape(arrayy,(nx*ny*nz))
#  reassembling and swapping first axis (y) with z axis
	arrayz = np.reshape(cell_array,(nx,ny,nz))
	arrayz = np.swapaxes(arrayz,0,2)
	arrayz = np.reshape(arrayz,(nx*ny*nz))
#  smoothing along z axis
	arrayz = smoothTriangle(arrayz,smoothing,dropVals=False)
#  reassembling
	arrayz = np.reshape(arrayz,(nz,ny,nx))
	arrayz = np.swapaxes(arrayz,0,2)
	arrayz = np.reshape(arrayz,(nx*ny*nz))
	target = (arrayx[:]+arrayy[:]+arrayz[:])/3.
	return target
	
	

