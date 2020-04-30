# $Id$
#
# read yaverages.dat file
#
# Author: J. Oishi (joishi@amnh.org). 
# 
#
__version__ = "$Id$"

import os
import re
import numpy as N

from ..files.npfile import npfile
from ..files.param import read_param 
from ..files.dim import read_dim 

def read_yaver(datadir='data/',format='native',point=(-1,-1)):

    """read 2D yaverage.dat file.

    point -- an array of 2-tuples (iz,ix) representing discrete
             points to be returned in an output array (not implemented yet)
    
    returns a tuple (yavg, t), yavg has shape (noutputs,nvars,nz,nx)

    """
    datadir = os.path.expanduser(datadir)
    datatopdir = re.sub('data\/*$','',datadir)
    filename = datadir+'/yaverages.dat'

    # which variables are averaged?
    infile = open(datatopdir+'yaver.in')
    variables = [line.strip() for line in infile.readlines()]
    infile.close()

    # global dim
    dim = read_dim(datadir) 
    if dim.precision == 'D':
        precision = 'd'
    else:
        precision = 'f'

    infile = npfile(filename,endian=format)
    
    t = N.zeros(1,dtype=precision)
    yaver = []
    yaver_shape = (len(variables),dim.nz,dim.nx)
    ntime = 0
    
    while 1:
        try:
            raw_data = infile.fort_read(precision,shape=1)
        except ValueError:
            break
        except TypeError:
            break

        t = N.concatenate((t,raw_data))

        try:
            raw_data = infile.fort_read(precision,shape=yaver_shape)
        except ValueError:
            #print "Problem: seems there is a t without corresponding data. yaverages.dat may be corrupted" # Python 2
            print("Problem: seems there is a t without corresponding data. yaverages.dat may be corrupted")
            break
        except TypeError:
            #print "Problem: seems there is a t without corresponding data. yaverages.dat may be corrupted" # Python 2
            print("Problem: seems there is a t without corresponding data. yaverages.dat may be corrupted")
            break
        yaver.append(raw_data)
        ntime += 1

    output = N.array(yaver)
    return output,t[1:]
