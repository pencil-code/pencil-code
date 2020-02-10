# qvar.py
#
# Read QVAR files. Based on the pc_read_qvar.pro IDL script, read_qvar.py Python script, and var.py Python script
#
# Author: Betsy Hernandez (bhern@gmail.com)
#
# modified read_class_nqvar function Mark Richardson (MR) (mrichardson@amnh.org) 05 Apr 18

#!/usr/bin/env python

import pencil_old as pc
import numpy as np
import os

def read_qvar(*args, **kwargs):
    """ Read QVAR files from Pencil Code in proc0.
        
        If iqvar is not provided read qvar.dat, otherwise read specified QVAR.

        Parameters:
            qvarfile=''
            datadir='data/'
            iqvar=-1
            verbose=False
    """
    return DataCube(*args, **kwargs)

class DataCube(object):
    def __init__(self, qvarfile='', datadir='data/', iqvar=-1, verbose=False):
    	
    	"""
    	Description:
    	-----------
    	Read QVAR files from Pencil Code in proc0. 
    	
    	If iqvar is not provided read qvar.dat, otherwise read QVAR specified.
    	
    	Params:
    	------
    	qvarfile=''
    	datadir='data/'
    	iqvar=-1
    	verbose=False
    	
    	Example of usage
    	------
    	qvar=pc.read_qvar(iqvar=100)
    	"""
    	
    	keys,places=get_qvarnames(datadir=datadir)
    	if (not qvarfile):
    		if iqvar < 0:
    			qvarfile = 'qvar.dat'
    		else:
    			qvarfile = 'QVAR'+str(iqvar)
    	ipars, qvars, t = read_class_nqvar(datadir=datadir, pfile=qvarfile, verbose=verbose)
    	setattr(self,'ipars',ipars)
    	for i in places:
    		setattr(self, keys[int(i) - 1][1:], qvars[int(i) - 1, :])
    	setattr(self, 't', t)

# function obtains the names and Fortran indices of the parameters in qvarname.dat
def get_qvarnames(datadir='data'):
	with open(datadir+'/qvarname.dat') as file:
		lines = [line.rstrip('\n') for line in file]
	keys = []
	places = []
	for line in lines:
		place, key = tuple(line.split())
		keys.append(key)
		places.append(place)
	return keys,places

# function obtains number of particles, matrix composed of parameters by particle number, and time
def read_class_nqvar(datadir='data',pfile='',proc=0,verbose=False):
	dims = pc.read_dim(datadir,proc)
	if (dims.precision == 'S'):
		REAL = '<f4'
	else:
		REAL = '<f8'
	qdims = pc.read_qdim(datadir)
	nqpars=qdims.nqpar
	mqvars=qdims.mqvar
	ltot = nqpars*mqvars
	if (verbose):
		print(nqpars+'massqvarve particles on processor: '+proc)
	# MR modified array_shape
	array_shape = np.dtype([('header', '<i4'),
							('ipar', '<i4'),
							('footer', '<i4'),
							('header2', '<i4'),
							('fq', REAL,ltot),
							('footer2', '<i4'),
							('header3', '<i4'),
							('t', REAL),
							('footer3', '<i4')])
	p_data = np.fromfile(datadir+'/proc0/'+pfile,dtype=array_shape)
	ipar = np.squeeze(p_data[0]['ipar'].reshape(p_data[0]['ipar'].size))
	partpars = np.array(p_data[0]['fq'].reshape(mqvars,nqpars))
	t=p_data[0]['t']
	return ipar,partpars,t

