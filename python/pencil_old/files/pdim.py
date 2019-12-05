# $Id: pdim.py  dhruba.mitra $
#
# read pdim.dat
#
# Author: D. Mitra (dhruba.mitra@gmail.com). based on idl pc_read_pdim.pro.
# 
# 
import sys
import os
import numpy as N

class PcPdim:
    """
    a class to hold the pdim.dat data.
    """
    def __init__(self,npar,mpvar,npar_stalk,mpaux):
        #primative quantities read directly from file
        self.npar = npar
        self.mpvar = mpvar
        self.npar_stalk = npar_stalk
        self.mpaux = mpaux

        
def read_pdim(datadir='data'):
    """
    read the pdim.dat file
    """
    
    filename = datadir+'/pdim.dat' # global particle properties

    try:
        filename = os.path.expanduser(filename)
        file = open(filename,"r")
    except IOError:
        #print "File",filename,"could not be opened." # Python 2
        print("File "+filename+" could not be opened.")
        return -1
    else:
        lines = file.readlines()[0].split()
        file.close()
        if N.size(lines) == 3:
            # old code - case: 'npar','mpvar','mpaux'	(before fall 2014)
            npar,mpvar,mpaux = tuple(map(int,lines))
            npar_stalk = 0
        if N.size(lines) == 4:
            # new code - case: 'npar','mpvar','npar_stalk','mpaux'
            npar,mpvar,npar_stalk,mpaux = tuple(map(int,lines))

    pdim = PcPdim(npar,mpvar,npar_stalk,mpaux)

    return pdim
