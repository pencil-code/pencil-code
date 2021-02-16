#!/usr/bin/env python

import pencil_old as pc
import numpy as np

from ..files.npfile import npfile
import os
import sys
#from string import maketrans

def read_pvar(*args, **kwargs):
	""" read pvar files from pencil code. if proc is not provided
		read all processors, if proc is provided read pvar.dat of 
		specified processor.
		
		Parameters:
			varfile=''    e.g. pvar.dat, PVAR9 and so on.
			datadir='./data'
			proc=-1           
			verbose=False    puts out some debugging output
			!!!!! DISCLAIMER!!!!! the next one is potentially dangerous and untested for a low number of particles 
			also, the number of particles put into this should be significantly lower than the original one
			reduce_to=-1  enter a positive integer number here, and the code will produce a subset of the particle files
						from the original pvar files. They are saved in the respective procX folder with the same name
						as the original file with the number of particles prepended (e.g. PVAR9 = 1000_PVAR9)
						This enables the use of subsets of particles in other simulations by e.g.
						using pc_copyvar to copy and subsequently rename the reduced particle files for new runs
	"""
	return pcpvar(*args, **kwargs)

class pcpvar(object):
		"""
		a class to hold the pvar data
		"""
		def __init__(self,
			varfile='pvar.dat',
			datadir='./data',
			proc=-1,
			verbose=False,
			reduce_to=-1):
			keys,places=get_pvarnames(datadir=datadir)
			for i in places:
				setattr(self,keys[int(i)-1].replace('(','P').replace(')','P'),int(i)-1)
			if (proc==-1):
				procdirs = list(filter(lambda s:s.startswith('proc'),os.listdir(datadir)))
				nprocs=len(procdirs)
				ipars,pvars = collect_class_pdata(pfile=varfile,datadir=datadir,nprocs=nprocs,
					verbose=verbose,reduce_to=reduce_to)
			else:
				ipars,pvars = read_class_npvar_red(datadir=datadir,proc=proc,
					verbose=verbose,reduce_to=reduce_to)
			
			setattr(self,'ipars',ipars)
			for i in places:
				setattr(self,keys[int(i)-1][1:].replace('(','P').replace(')','P'),pvars[int(i)-1,:])
				
				
				
def read_npar_loc(datadir='./data',pfile='pvar.dat',proc=0):
	array_npar_loc=np.dtype([('header','<i4')])
	npars = np.fromfile(datadir+'/proc'+str(proc)+'/'+pfile,dtype=array_npar_loc)
	return npars['header'][1]

def collect_pdata(datadir='./data'):
	procs = np.arange(32)
	for i in procs:
		dom_ipar,dom_part_pos = read_npvar_red(datadir=datadir,proc=i)
		if i == 0:
			ipars = dom_ipar
			part_pos = dom_part_pos
		else:
			ipars = np.hstack((ipars,dom_ipar))
			part_pos = np.hstack((part_pos,dom_part_pos))

	return ipars,part_pos


def get_pvarnames(datadir='./data'):
	with open(datadir+'/pvarname.dat') as file:
		lines = [line.rstrip('\n') for line in file]
	keys = []
	places = []
	for line in lines:
		place, key = tuple(line.split())
		keys.append(key)
		places.append(place)
	return keys,places

def read_class_npvar_red(datadir='./data',
	pfile='pvar.dat',
	proc=0,
	verbose=False,
	reduce_to=-1,
	set_reduce=-1):
		
	dims = pc.read_dim(datadir,proc)
	pdims = pc.read_pdim(datadir)
	npar_loc = read_npar_loc(datadir=datadir,pfile=pfile,proc=proc)
	#
	# the Next bit calculates how many particles are written for all but
	# the last processor. The last processor is assigned a number of particles
	# to write so that the required number of particles is obtained
	#
	if (reduce_to > 0):
		if (set_reduce <= 0):
			reductionfactor = float(reduce_to)/float(pdims.npar)
			npar_red = int(round(npar_loc*reductionfactor))
		else:
			npar_red = set_reduce
		if (verbose):
			#print 'reducing '+str(npar_loc)+' to '+str(npar_red)+ ' on proc'+str(proc) 
			print('reducing {} to {} on proc {}'.format(
                              npar_loc, npar_red, proc))
		written_parts=npar_red
	else:
		written_parts=set_reduce
	
	if (verbose):
		#print npar_loc,' particles on processor: ',proc # Python 2
		print(str(npar_loc)+' particles on processor: '+str(proc))
	mvars = pdims.mpaux+pdims.mpvar
	ltot = npar_loc*mvars
	if (dims.precision == 'S'):
		REAL = '<f4'
	else:
		REAL = '<f8'

	array_shape= np.dtype([('header','<i4'),
							('npar_loc','<i4'),
							('footer','<i4'),
							('header2','<i4'),
							('ipar','<i4',npar_loc),
							('footer2','<i4'),
							('header3','<i4'),
							('fp',REAL,ltot),
							('footer3','<i4'),
							('header4','<i4'),
							('t',REAL),
							('x',REAL,dims.mx),
							('y',REAL,dims.my),
							('z',REAL,dims.mz),
							('dx',REAL),
							('dy',REAL),
							('dz',REAL),
							('footer4','<i4')])
	
	
	p_data = np.fromfile(datadir+'/proc'+str(proc)+'/'+pfile,dtype=array_shape)
	partpars = np.array(p_data['fp'].reshape(mvars,npar_loc))
	
	if (reduce_to>0):
		particle_list = map(lambda x: int(x),np.linspace(0.0,npar_loc,num=npar_red,endpoint=False))
		red_parts = partpars[:,particle_list]
		red_shape = np.dtype([('header','<i4'),
							('npar_loc','<i4'),
							('footer','<i4'),
							('header2','<i4'),
							('ipar','<i4',(npar_red),),
							('footer2','<i4'),
							('header3','<i4'),
							('fp',REAL,npar_red*mvars),
							('footer3','<i4'),
							('header4','<i4'),
							('t',REAL),
							('x',REAL,(dims.mx,)),
							('y',REAL,(dims.my,)),
							('z',REAL,(dims.mz,)),
							('dx',REAL),
							('dy',REAL),
							('dz',REAL),
							('footer4','<i4')])
							
		p_red =np.array([(4,
			npar_red,
			4,
			(npar_red*4),
			(np.squeeze(p_data['ipar'][0,:npar_red])),
			(npar_red*4),
			(npar_red*mvars*8),
			(np.squeeze(np.ravel(red_parts))),
			(npar_red*mvars*8),
			(p_data['header4'][0]),
			(p_data['t']),
			(p_data['x']),
			(p_data['y']),
			(p_data['z']),
			(p_data['dx']),
			(p_data['dy']),
			(p_data['dz']),
			p_data['footer4'][0])
			],dtype=red_shape)
			
		p_red.tofile(datadir+'/proc'+str(proc)+'/'+str(reduce_to)+'_'+pfile)
		
	ipar = np.squeeze(p_data['ipar'].reshape(p_data['ipar'].size))
	return ipar, partpars, written_parts
	
def collect_class_pdata(pfile='pvar.dat',datadir='/data',nprocs='0',verbose=False,reduce_to=-1):
	if (nprocs==0):
		#print "this should be greater than zero" # Python 2
		print("this should be greater than zero")
	else:
		procs=range(nprocs)
		
	if (reduce_to > 0):
		parts_to_write=reduce_to
		particle_sum=0
	for i in procs:
		if (not i==procs[-1]):
			dom_ipar,dom_pvar, written_parts = read_class_npvar_red(pfile=pfile,datadir=datadir,proc=i,reduce_to=reduce_to,verbose=verbose)
			if (reduce_to>0):
				parts_to_write -= written_parts
				particle_sum+=written_parts
		else:			
			if (reduce_to>0):
				dom_ipar,dom_pvar, written_parts = read_class_npvar_red(pfile=pfile,datadir=datadir,proc=i,reduce_to=reduce_to,verbose=verbose,set_reduce=parts_to_write)
				particle_sum+=written_parts
				if (verbose):
					print('written particles',particle_sum)
			else:
				dom_ipar,dom_pvar, written_parts = read_class_npvar_red(pfile=pfile,datadir=datadir,proc=i,reduce_to=reduce_to,verbose=verbose)
		if i == 0:
			ipars = dom_ipar
			pvars = dom_pvar
		else:
			ipars = np.hstack((ipars,dom_ipar))
			pvars = np.hstack((pvars,dom_pvar))
		if (verbose):
			#print 'Reading processor '+ str(i)+'.' # Python 2
			print('Reading processor '+ str(i)+'.')
	return ipars,pvars
	

	
	

