#!/usr/bin/env python

import pencil as pc
import numpy as np

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
        def __init__(self,varfile='pvar.dat',casedir='.',datadir='/data',proc=-1,verbose=False):
            keys,places=get_pvarnames(casedir=casedir,datadir=datadir)
            for i in places:
                setattr(self,keys[int(i)-1].replace('(','P').replace(')','P'),int(i)-1)
            if (proc==-1):
                procdirs = filter(lambda s:s.startswith('proc'),os.listdir(casedir+datadir))
                nprocs=len(procdirs)
                ipars,pvars = collect_class_pdata(casedir=casedir,pfile=varfile,datadir=datadir,nprocs=nprocs,verbose=verbose)
            else:
                ipars,pvars = read_class_npvar_red(casedir=casedir,datadir=datadir,proc=proc,verbose=verbose)
            
            setattr(self,'ipars',ipars)
            for i in places:
                setattr(self,keys[int(i)-1][1:].replace('(','P').replace(')','P'),pvars[int(i)-1,:])
                
                
                
def read_npar_loc(casedir='.',datadir='/data',pfile='pvar.dat',proc=0):
    array_npar_loc=np.dtype([('header','<i4'),('npar_loc','<i4')])
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

def read_class_npvar_red(casedir='.',datadir='/data',pfile='pvar.dat',proc=0,verbose=False):
    
    dims = pc.read_dim(casedir+datadir,proc)
    pdims = pc.read_pdim(casedir+datadir)
    npar_loc = read_npar_loc(casedir=casedir,datadir=datadir,pfile=pfile,proc=proc)
    if (verbose):
		print npar_loc,' particles on processor: ',proc
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
    
    p_data = np.fromfile(casedir+datadir+'/proc'+str(proc)+'/'+pfile,dtype=array_shape)
    partpars = np.array(p_data['fp'].reshape(mvars,npar_loc))
    ipar = np.squeeze(p_data['ipar'].reshape(p_data['ipar'].size))
    return ipar,partpars
    
def collect_class_pdata(casedir='.',pfile='pvar.dat',datadir='/data',nprocs='0',verbose=False):
    if (nprocs==0):
        print "this should be greater than zero"
    else:
        procs=range(nprocs)
    for i in procs:
		dom_ipar,dom_pvar = read_class_npvar_red(casedir=casedir,pfile=pfile,datadir=datadir,proc=i)
		if i == 0:
			ipars = dom_ipar
			pvars = dom_pvar
		else:
			ipars = np.hstack((ipars,dom_ipar))
			pvars = np.hstack((pvars,dom_pvar))
		if (verbose):
			print 'Reading processor '+ str(i)+'.'
    return ipars,pvars
    

    
    

