#!/usr/bin/env python

import pencil as pc
import numpy as np

from pencil.files.npfile import npfile
import os
import sys
from scipy.io import FortranFile

def read_pstalk(*args, **kwargs):
	return pstalk(*args, **kwargs)

class pstalk(object):
	def __init__(self,
	varfile='particles_stalker.dat',
	datadir='./data',
	proc=-1,
	verbose=False,
	reduce_to=-1):
		
		
		pdata = pc.read_pdim()
		setattr(self,'nstalk',pdata.npar_stalk)
		
		self.find_no_of_procs(proc=proc,datadir=datadir)
		self.get_stalkdim(datadir=datadir)
		self.get_stalktime(datadir=datadir)
		
		for proc in self.procdirs:
			self.read_stalker_data(datadir=datadir,procdir=proc)
			
		self.bundle_stalker_data()
	
	def bundle_stalker_data(self,):
		snap_data = getattr(self,'proc0_snaptime')
		full_data = np.zeros((np.shape(snap_data)[0],self.nstalk+1,self.ndims_stalk))
		for proc in self.procdirs:
			data = getattr(self,proc+'_partdata')
			iddata = getattr(self,proc+'_iddata')
			nr_shot_data = getattr(self,proc+'_nrshot')
			for i in np.arange(np.shape(snap_data)[0]):
				if (nr_shot_data[i]> 0):
					for j in np.arange(np.shape(iddata[i])[0]):
						#~ print data[i][j]
						#~ print i, iddata[i][j]
						#~ print '--------------'
						full_data[i,iddata[i][j],:] = data[i][j]
		full_data[full_data==0.0] = float('NaN')
		setattr(self,'bundled_data',full_data)
					
			#~ data = np.asarray(np.reshape(data,(self.n_stalk_shots,-1,self.ndims_stalk)))
			#~ collect_array = np.concatenate((collect_array,data),axis=1)
			#~ collect_id = np.concatenate((collect_id, iddata),axis=1)
			
	def find_no_of_procs(self,proc,datadir):
		if (proc==-1):
			setattr(self,'procdirs',list(filter(lambda s:s.startswith('proc'),os.listdir(datadir))))
			setattr(self,'nprocs',len(self.procdirs))

	def get_stalkdim(self,datadir):
		with open(datadir+'/particles_stalker_header.dat') as file:
			lines = [line.rstrip('\n') for line in file]
			rmline = lines
		for line in rmline:
			var_stalk = line.split(',')[:-1]
		ndims_stalk = np.shape(var_stalk)[0]
		setattr(self,'ndims_stalk',int(ndims_stalk))
		setattr(self,'var_stalk',var_stalk)
		
	
	def get_stalktime(self,datadir):
		with open(datadir+'/tstalk.dat') as file:
			lines = [line.rstrip('\n') for line in file]
			rmline = lines
		for line in rmline:
			t_stalk = line.split()
			t_last_stalk = t_stalk[0]
			n_stalk_shots = t_stalk[1]
		setattr(self,'t_last_stalk',t_last_stalk)
		setattr(self,'n_stalk_shots',int(n_stalk_shots)-2)
			
			
	def read_stalker_data(self,
		datadir='/data',
		procdir='proc0'):
			
		f = FortranFile(datadir+'/'+procdir+'/particles_stalker.dat','r')
		file_not_ended=True
		firstline = f.read_record([('a','f8'),('b','i4')])
		firstline = f.read_record([('a','f8'),('b','i4')])
		snaptime = []
		nrshot = []
		partdata = []
		iddata = []
		while (True):
			try:
				temp = f.read_record([('snaptime','f8'),('nrshot','i4')])
				snaptime.append(temp['snaptime'])
				nrshot.append(temp['nrshot'])
			except TypeError:
				break
			except ValueError:
				break
			if (nrshot[-1] > 0):
				temp = f.read_record([('iddata','i4')])
				iddata.append(temp['iddata'])
				temp = f.read_record([('partdata','f8')])
				temp = np.reshape(temp,(-1,self.ndims_stalk))
				partdata.append(temp['partdata'])
			else:
				iddata.append([])
				partdata.append([])
		
		partdata=np.asarray(partdata)
		iddata = np.asarray(iddata)
		snaptime = np.asarray(snaptime)
		nrshot = np.array(nrshot)
		oldshape = np.shape(partdata)
		
		partdata = np.reshape(partdata,oldshape)
		
		setattr(self,procdir+'_snaptime',snaptime)
		setattr(self,procdir+'_iddata',iddata)
		setattr(self,procdir+'_partdata',partdata)
		setattr(self,procdir+'_nrshot',nrshot)
				
				
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
	

	
	

