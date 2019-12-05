# $Id$

import numpy as N
from ..files.dim import read_dim

def read_power(file):
    """ 
    29-apr-2009/dintrans: coded
    t,dat=read_power(name_power_file)
    Read a power spectra file like 'data/poweru.dat'
    """ 
    infile = open(file, 'r')
    lines = infile.readlines()
    infile.close()
#
#  find the number of blocks (t,power) that should be read
#
    dim=read_dim()
    #nblock=len(lines)/int(N.ceil(dim.nxgrid/2/8.)+1) # Python 2 integer division
    nblock=int(len(lines)/int(N.ceil(int(dim.nxgrid/2)/8.)+1))
#
    infile = open(file, 'r')
    t=N.zeros(1, dtype='Float32')
    data=N.zeros(1, dtype='Float32')
    for i in range(nblock):
        st=infile.readline()
        t=N.append(t, float(st))
        #for ii in range(int(N.ceil(dim.nxgrid/2/8.))): # Python 2 integer division
        for ii in range(int(N.ceil(int(dim.nxgrid/2)/8.))):
            st=infile.readline()
            data=N.append(data, N.asarray(st.split()).astype('f'))
    infile.close()

    t=t[1:] ; data=data[1:]
    #nt=len(t) ; nk=len(data)/nt # Python 2 integer division
    nt=len(t) ; nk=int(len(data)/nt)
    data=data.reshape(nt, nk)
    return t, data

