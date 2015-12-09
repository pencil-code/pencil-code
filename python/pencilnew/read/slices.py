# $Id$
#
# read slice files
#
# Author: J. Oishi (joishi@amnh.org). 
# 
#
import os
import numpy as N
import pylab as P
from pencilnew.io.npfile import npfile
from .param import param as read_param 
from .dim import dim as read_dim 
from time import sleep 

import sys

# slice file format is either
#   plane,t (old style)
#   plane,t,slice_z2pos (new style)


####### MISSING: if data/slice* are not available, ask to run pc_collectallmovie for the user
####### MISSING: default for extension and field should be to ask the user if not specified
####### MISSING: slice should be given back as object!

def slices(field='uu1',datadir='./data',proc=-1,extension='xz',format='native',oldfile=False):
    """
    read 2D slice files and return an array of (nslices,vsize,hsize).
    """
    datadir = os.path.expanduser(datadir)
    if proc < 0:
        filename = datadir+'/slice_'+field+'.'+extension
    else:
        filename = datadir+'/proc'+str(proc)+'/slice_'+field+'.'+extension

    # global dim
    param = read_param(datadir, quiet=True)

    dim = read_dim(datadir,proc) 
    if dim.precision == 'D':
        precision = 'd'
    else:
        precision = 'f'

    # set up slice plane
    if (extension == 'xy' or extension == 'Xy'):
        hsize = dim.nx
        vsize = dim.ny
    if (extension == 'xz'):
        hsize = dim.nx
        vsize = dim.nz
    if (extension == 'yz'):
        hsize = dim.ny
        vsize = dim.nz


    infile = npfile(filename,endian=format)

    ifirst = True
    islice = 0
    t = N.zeros(1,dtype=precision)
    slices = N.zeros(1,dtype=precision)

    while 1:
        try:
            raw_data = infile.fort_read(precision)
        except ValueError:
            break
        except TypeError:
            break
        
        if oldfile:
            t = N.concatenate((t,raw_data[-1:]))
            slices = N.concatenate((slices,raw_data[:-1]))
        else:
            t = N.concatenate((t,raw_data[-2:-1]))
            slices = N.concatenate((slices,raw_data[:-2]))
        islice += 1

    output = slices[1:].reshape(islice,vsize,hsize)

    return output,t[1:]
  