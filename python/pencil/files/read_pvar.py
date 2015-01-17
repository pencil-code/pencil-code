#!/usr/bin/env python

import pencil as pc
import numpy as np
import pylab as pyl

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
			
			
def read_npar_loc(datadir='/data',pfile='pvar.dat',proc=0):
	array_npar_loc=np.dtype([('header','<i4'),
							('npar_loc','<i4')])
	npars = np.fromfile(datadir+'/proc'+str(proc)+'/'+pfile,dtype=array_npar_loc)
	return npars['npar_loc'][0]

def read_npvar(datadir='./data',proc=0):
	
	dims = pc.read_dim(datadir,proc)
	pdims = pc.read_pdim(datadir)
	npar_loc = read_npar_loc(datadir,pfile,proc)
	print npar_loc,' particles on processor: ',proc
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
							('x','<f8',26),
							('y','<f8',16),
							('z','<f8',26),
							('dx','<f8'),
							('dy','<f8'),
							('dz','<f8'),
							('footer4','<i4')])
	
	p_data = np.fromfile(datadir+'/proc'+str(proc)+'/'+pfile,dtype=array_shape)
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
	footer = p_data['footer']
	npar_loc = np.squeeze(p_data['npar_loc'])
	return ipar,part_pos,procdims
	
	
def calc_alpha(datadir='./data',pfile='pvar.dat',proc=0,n_cells=4):
	
	ipar,part_pos,procdims = read_npvar(datadir,pfile,proc)
	dims = pc.read_dim(datadir,proc)
	lx = procdims[0,0]-procdims[1,0]	
	ly = procdims[0,1]-procdims[1,1]
	lz = procdims[0,2]-procdims[1,2]
	dx = lx/n_cells
	dy = ly/n_cells
	dz = lz/n_cells
	sx = procdims[1,0]
	sy = procdims[1,1]
	sz = procdims[1,2]
	
	cell_array = np.zeros(n_cells**3)
	print cell_array
	
	for i in np.arange(len(ipar)):
		for j in np.arange(n_cells-1):
			for k in np.arange(n_cells-1):
				for l in np.arange(n_cells-1):
					if (sx + dx*j < part_pos[0,i] <= sx + dx*(j+1) and \
						sy + dy*k < part_pos[1,i] <= sy + dy*(k+1) and \
						sz + dz*l < part_pos[2,i] <= sz + dz*(l+1)):
						print 'located particle '+ str(ipar[i]) + ' in cell ' + str(l*(n_cells**2) + (k*n_cells)+j) 
						cell_array[l*(n_cells**2) + (k*n_cells)+j] += 1
								
	n, bins, patches = pyl.hist(cell_array, 50, normed=1, histtype='stepfilled')
