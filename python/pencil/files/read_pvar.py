#!/usr/bin/env python

import pencil as pc
import numpy as np
import pylab as pyl

from npfile import npfile
import os
import sys
from string import maketrans

def read_pvar(*args, **kwargs):
	""" read pvar files from pencil code. if proc is not provided
		read all processors, if proc is provided read pvar.dat of 
		specified processor.
		
		Parameters:
			varfile=''
			datadir='data/'
			proc=-1
			casedir='.'
	"""
	return pcpvar(*args, **kwargs)

class pcpvar(object):
		"""
		a class to hold the pvar data
		"""
		def __init__(self,varfile='',casedir='.',datadir='/data',proc=-1):
			keys,places=get_pvarnames(casedir=casedir,datadir=datadir)
			for i in places:
				setattr(self,keys[int(i)-1].replace('(','P').replace(')','P'),int(i)-1)
			if (proc==-1):
				procdirs = filter(lambda s:s.startswith('proc'),os.listdir(casedir+datadir))
				nprocs=len(procdirs)
				ipars,pvars = collect_class_pdata(casedir=casedir,datadir=datadir,nprocs=nprocs)
			else:
				ipars,pvars = read_class_npvar_red(casedir=casedir,datadir=datadir,proc=proc)
			
			setattr(self,'ipars',ipars)
			for i in places:
				setattr(self,keys[int(i)-1][1:].replace('(','P').replace(')','P'),pvars[int(i)-1,:])
				
				
				
def read_npar_loc(casedir='.',datadir='/data',pfile='pvar.dat',proc=0):
	array_npar_loc=np.dtype([('header','<i4'),
							('npar_loc','<i4')])
	npars = np.fromfile(casedir+datadir+'/proc'+str(proc)+'/'+pfile,dtype=array_npar_loc)
	return npars['npar_loc'][0]


	
	


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


def get_pvarnames(casedir='.',datadir='/data'):
	with open(casedir+datadir+'/pvarname.dat') as file:
		lines = [line.rstrip('\n') for line in file]
	keys = []
	places = []
	for line in lines:
		place, key = tuple(line.split())
		keys.append(key)
		places.append(place)
	return keys,places

def read_class_npvar_red(casedir='.',datadir='/data',pfile='pvar.dat',proc=0):
	
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
	partpars = np.array(p_data['fp'].reshape(mvars,npar_loc))
	ipar = np.squeeze(p_data['ipar'].reshape(p_data['ipar'].size))
	return ipar,partpars
	
def collect_class_pdata(casedir='.',datadir='/data',nprocs='0'):
	if (nprocs==0):
		print "this should be greater than zero"
	else:
		procs=range(nprocs)
	
	for i in procs:
		dom_ipar,dom_pvar = read_class_npvar_red(casedir=casedir,datadir=datadir,proc=i)
		if i == 0:
			ipars = dom_ipar
			pvars = dom_pvar
		else:
			ipars = np.hstack((ipars,dom_ipar))
			pvars = np.hstack((pvars,dom_pvar))

	return ipars,pvars
	

	
	

