# $Id: qdim.py $
#
# read qdim.dat
#
# Author: Betsy Hernandez (bhern@gmail.com) based on idl pc_read_qdim.pro and python pdim.py.
# 
# 
import sys
import os
import numpy as N

class PcQdim:
    """
    a class to hold the qdim.dat data.
    """
    def __init__(self,nqpar,mqvar):
        #primative quantities read directly from file
        self.nqpar = nqpar
        self.mqvar = mqvar
        
def read_qdim(datadir='data'):
    """
    read the qdim.dat file
    """
    
    filename = datadir+'/qdim.dat' # global particle properties

    try:
        filename = os.path.expanduser(filename)
        file = open(filename,"r")
    except IOError:
        print("File "+filename+" could not be opened.")
        return -1
    else:
        lines = file.readlines()[0].split()
        file.close()
        if N.size(lines) == 2:
            # code - case: 'nqpar','mqvar'
            nqpar,mqvar = tuple(map(int,lines))

    qdim = PcQdim(nqpar,mqvar)

    return qdim

