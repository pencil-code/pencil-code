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
    def __init__(self,npar,mpvar,npar_stalk):
        #primative quantities read directly from file
        self.npar = npar
        self.mpvar = mpvar
        self.npar_stalk = npar_stalk

        
def read_pdim(datadir='data',proc=-1):
    """
    read the pdim.dat file. if proc is -1, then read the 'global'
    dimensions. if proc is >=0, then read the pdim.dat in the
    corresponding processor directory.
    """
    if (proc < 0 ):
        filename = datadir+'/pdim.dat' # global box dimensions
    else:
        filename = datadir+'/proc'+str(proc)+'/pdim.dat' # local proc. dimensions

    try:
	filename = os.path.expanduser(filename)
        file = open(filename,"r")
    except IOError:
        print "File",filename,"could not be opened."
        return -1
    else:
        lines = file.readlines()
        file.close()
        # 'npar','mpvar','npar_stalk'
	npar,mpvar,npar_stalk = tuple(map(int,lines[0].split()))

    pdim = PcPdim(npar,mpvar,npar_stalk)

    return pdim
