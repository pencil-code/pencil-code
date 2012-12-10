# $Id: zaver.py 10300 2009-01-21 22:57:47Z wdobler $
#
# read zaverages.dat file
#
# Author: J. Oishi (joishi@amnh.org). 
# 
#
__version__ = "$Id: zaver.py  $"

import os
import re
import numpy as N

from npfile import npfile
from param import read_param 
from dim import read_dim 

def read_zaver(datadir='data/',format='native',point=(-1,-1)):

    """read 2D zaverage.dat file.

    point -- an array of 2-tuples (iy,ix) representing discrete
             points to be returned in an output array (not implemented yet)
    
    returns a tuple (zavg, t), zavg has shape (noutputs,nvars,ny,nx)

    """
    datadir = os.path.expanduser(datadir)
    datatopdir = re.sub('data\/*$','',datadir)
    filename = datadir+'proc0/zaverages.dat'

    # which variables are averaged?
    infile = open(datatopdir+'zaver.in')
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
    zaver = []
    zaver_shape = (len(variables),dim.ny,dim.nx)
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
            raw_data = infile.fort_read(precision,shape=zaver_shape)
        except ValueError:
            print "Problem: seems there is a t without corresponding data. zaverages.dat may be corrupted"
            break
        except TypeError:
            print "Problem: seems there is a t without corresponding data. zaverages.dat may be corrupted"
            break
        zaver.append(raw_data)
        ntime += 1

    output = N.array(zaver)
    return output,t[1:]
