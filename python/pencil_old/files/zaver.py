# $Id: zaver.py 10300 2009-01-21 22:57:47Z wdobler $
#
# read zaverages.dat file
#
# Author: J. Oishi (joishi@amnh.org). 
# Fred: added options to select by time range and/or index ranges
#
#
__version__ = "$Id: zaver.py  $"

import os
import re
import numpy as N

from ..files.npfile import npfile
from ..files.param import read_param 
from ..files.dim import read_dim 

def read_zaver(datadir='data/',format='native',point=(-1,-1),proc=-1,
               trange=(0,None),tindex=(0,None,1)):

    """read 2D zaverage.dat file. If proc < 0, then load all data
    and assemble. Otherwise, load VAR file from specified processor.

    point -- an array of 2-tuples (iy,ix) representing discrete
             points to be returned in an output array (not implemented yet)
    
    proc -- Read data from proc if proc > -1, otherwise load all and assemble.

    trange -- read subset of data between time.min(), time.max() (0,None)

    tindex -- read every 1 data between indices 0 and None   
    
    returns a tuple (zavg, t), zavg has shape (noutputs,nvars,ny,nx)

    """
    datadir = os.path.expanduser(datadir)
    datatopdir = re.sub('data\/*$','',datadir)

    if len(trange) != 2:
        print("Warning: trange must be a real/integer tuple of length 2 with"+
              " start and end time specified")
    if len(tindex) != 3:
        print("Warning: tindex must be an integer tuple of length 3 with"+
              " start and end indices and step specified")

    # which variables are averaged?
    infile = open(datatopdir+'zaver.in')
    variables = [line.strip() for line in infile.readlines()]
    infile.close()

    # global dim
    dim = read_dim(datadir, proc=proc)
    if dim.precision == 'D':
        precision = 'd'
    else:
        precision = 'f'

    if proc < 0:
        procdirs = list(filter(lambda s:s.startswith('proc') and not s.endswith('.dat'),
                                 os.listdir(datadir)))
    else:
        procdirs = ['proc'+str(proc)]

    for directory in procdirs:
        ntime,ndx = 0,0
        # local dimensions
        core = int(directory[4:]) # SC: needed to rename proc to core to keep function argument
        procdim = read_dim(datadir,core)
        nxloc = procdim.nx
        nyloc = procdim.ny
        zaver_local = []
        zaver_loc_shape = (len(variables),procdim.ny,procdim.nx)

        #read data
        filename = os.path.join(datadir,directory,'zaverages.dat')
        try:
            infile = npfile(filename,endian=format)
            t = N.zeros(1,dtype=precision)
        except:
            continue
        while 1:
            try:
                raw_data = infile.fort_read(precision,shape=1)
            except ValueError:
                break
            except TypeError:
                break

            if ndx >= tindex[0] and N.mod(ndx,tindex[2])==0:
                if raw_data >= trange[0]:
                    t = N.concatenate((t,raw_data))
            try:
                raw_data = infile.fort_read(precision,shape=zaver_loc_shape)
            except ValueError:
                #print "Problem: seems there is a t without corresponding data. zaverages.dat may be corrupted" # Python 2                
                print("Problem: seems there is a t without corresponding data. zaverages.dat may be corrupted")
                break
            except TypeError:
                #print "Problem: seems there is a t without corresponding data. zaverages.dat may be corrupted" # Python
                print("Problem: seems there is a t without corresponding data. zaverages.dat may be corrupted")
                break
            if t.size-1 > ntime and t[-1] >= trange[0]:
                zaver_local.append(raw_data)
                ntime += 1
            ndx += 1
            if tindex[1] is not None and ndx > tindex[1]:
                break
            if trange[1] is not None and t[-1] > trange[1]:
                break

        try:
            zaver
            pass
        except:
            zaver = N.zeros((ntime,len(variables),dim.ny,dim.nx))
        
        if (proc < 0):
            # append to the global zaver
            for i in range(ntime):
                zaver[i,:,procdim.ipy*procdim.ny:(procdim.ipy+1)*procdim.ny,
                    procdim.ipx*procdim.nx:(procdim.ipx+1)*procdim.nx] = zaver_local[i]
        else:
            for i in range(ntime):
                zaver[i,:,:,:] = zaver_local[i]
        
    return zaver,t[1:]
