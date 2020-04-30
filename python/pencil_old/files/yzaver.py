# $Id$
#
# read yzaver files
#
# Author: J. Oishi (joishi@amnh.org).
#
#
import os
import re
# SC: commands i depracated since python 2.6
#from commands import getoutput
from subprocess import check_output
import numpy as N

#from param import read_param
from ..files.dim import read_dim

class YzAver:
    pass

def read_yzaver(varfile='yzaverages.dat',datadir='data/'):
    """read yz averaged data.

    returns a class yzaver containing a 1D time array and a 2D t-z array.

    """
    datadir = os.path.expanduser(datadir)
    dim = read_dim(datadir=datadir)
    nx = dim.nx
    datatopdir = re.sub('data\/*$','',datadir)

    infile = open(datatopdir+'yzaver.in')
    variables = [line.strip() for line in infile.readlines()]
    infile.close()

    #gotta be a better way to do this...
    n_lines = int(check_output(['wc', datadir+'/'+varfile]).split()[0])

    datafile = open(datadir+'/'+varfile)
    n_vars = len(variables)
    # each line of data has 8 values
    rec_length = 1 + int(n_vars*nx/8)
    if nx%8:
        rec_length += 1
    n_data_records = n_lines/rec_length
    #print "%s: reading %i records" % (__name__,n_data_records) # Python 2
    print(__name__ + ": reading {0} records".format(n_data_records))

    # change the hardcode dtype!
    t = []
    var_tmp = dict(zip(variables,[[] for var in variables]))
    for i in range(int(n_data_records)):
        t.extend(N.fromfile(datafile,dtype='float32',count=1,sep=' '))
        for var in variables:
            var_tmp[var].append(N.fromfile(datafile,
                                         dtype='float32',count=nx,sep=' '))

    datafile.close()
    #pack data into a class
    yzaver = YzAver()
    yzaver.t = N.array(t)

    for var in variables:
        setattr(yzaver,var,N.array(var_tmp[var]))
    return yzaver
