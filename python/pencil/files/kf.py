# $Id: kf.py 13777 2010-05-02 22:12:43Z iomsn $
#
# Read the forcing wavenumber kf from the k.dat file.
#
# Author: Simon Candelaresi (iomsn@physto.se, iomsn1@googlemail.com).
# 
#

import numpy as np


def read_kf(datadir = '.'):
    """
    Read the forcing wavenumber kf from the k.dat file and return its value.
    
    call signature::
    
      kf = read_kf(datadir = '.')
    
    Keyword arguments:

      *datadir*:
        The directory where the k.dat file resides.
    """
    
    
    # default value in case kf is not set
    kf = 0.

    kfile = open(datadir+'/k.dat', 'r')
    line = kfile.readline()
    kf = np.float(str.split(line)[1])
    
    kfile.close()
    
    return kf
    