#!/usr/bin/env python

import pencil_old as pc
import numpy as np

from pencil_old.files.npfile import npfile
import os
import sys
#from string import maketrans

def read_qvar(*args, **kwargs):
    """ read qvar files from pencil code. if proc is not provided
        read all processors, if proc is provided read qvar.dat of 
        specified processor.
        
        Parameters:
            varfile=''
            datadir='data/'
            proc=-1
    """
    return pcqvar(*args, **kwargs)

class pcqvar(object):
        """
        a class to hold the qvar data
        """
        def __init__(self,varfile='qvar.dat',datadir='./data',proc=-1,verbose=False):
            keys,places=get_qvarnames(datadir=datadir)
            for i in places:
                setattr(self,keys[int(i)-1].replace('(','P').replace(')','P'),int(i)-1)
            if (proc==-1):
                procdirs = list(filter(lambda s:s.startswith('proc') and not s.endswith('.dat'),os.listdir(datadir)))
                nprocs=len(procdirs)
                ipars,pvars = collect_class_pdata(pfile=varfile,datadir=datadir,nprocs=nprocs,verbose=verbose)
            else:
                ipars,pvars = read_class_npvar_red(datadir=datadir,proc=proc,verbose=verbose)
            
            setattr(self,'ipars',ipars)
            for i in places:
                setattr(self,keys[int(i)-1][1:].replace('(','P').replace(')','P'),pvars[int(i)-1,:])
                                
def read_npar_loc(datadir='./data',pfile='qvar.dat',proc=0):
    array_nqpar=np.dtype([('header','<i4'),('nqpar','<i4')])
    nqpar = np.fromfile(datadir+'/proc'+str(proc)+'/'+pfile,dtype=array_nqpar)
    return npars['nqpar'][0]

def collect_pdata(datadir='./data'):
    procs = np.arange(32)
    for i in procs:
        dom_ipar,dom_part_pos = read_nqvar_red(datadir=datadir,proc=i)
        if i == 0:
            ipars = dom_ipar
            part_pos = dom_part_pos
        else:
            ipars = np.hstack((ipars,dom_ipar))
            part_pos = np.hstack((part_pos,dom_part_pos))

    return ipars,part_pos


def get_qvarnames(datadir='./data'):
    with open(datadir+'/qvarname.dat') as file:
        lines = [line.rstrip('\n') for line in file]
    keys = []
    places = []
    for line in lines:
        place, key = tuple(line.split())
        keys.append(key)
        places.append(place)
    return keys,places

def read_class_nqvar_red(datadir='./data',pfile='qvar.dat',proc=0,verbose=False):
    dims = pc.read_dim(datadir,proc)
    qdims = pc.read_qdim(datadir)
    nqpar = read_qpar(datadir=datadir,pfile=pfile,proc=proc)
    if (verbose):
        #print npar_loc,' particles on processor: ',proc # Python 2
        print(nqpar+'massive particles on processor: '+proc)
    mvars = pdims.mqaux+pdims.mqvar
    ltot = nqpar*mvars
    if (dims.precision == 'S'):
        REAL = '<f4'
    else:
        REAL = '<f8'

    array_shape= np.dtype([('header','<i4'),
                            ('nqpar','<i4'),
                            ('footer','<i4'),
                            ('header2','<i4'),
                            ('ipar','<i4',nqpar),
                            ('footer2','<i4'),
                            ('header3','<i4'),
                            ('fq',REAL,ltot),
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
    partpars = np.array(p_data['fq'].reshape(mvars,npar_loc))
    ipar = np.squeeze(p_data['ipar'].reshape(p_data['ipar'].size))
    return ipar,partpars
    
def collect_class_pdata(pfile='qvar.dat',datadir='/data',nprocs='0',verbose=False):
    if (nprocs==0):
        #print "this should be greater than zero" # Python 2
        print("this should be greater than zero")
    else:
        procs=range(nprocs)
    for i in procs:
        dom_ipar,dom_pvar = read_class_nqvar_red(pfile=pfile,datadir=datadir,proc=i)
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
    

    
    

