# $Id$
#
# read dim.dat
#
# Author: J. Oishi (joishi@amnh.org). based on idl pc_read_dim.pro.
# 
# 
import sys
import os
import numpy as N

class PcDim:
    """
    a class to hold the dim.dat data.
    """
    def __init__(self,mx,my,mz,mvar,maux,precision,nghostx,nghosty,nghostz,nprocx,nprocy,nprocz,iprocz_slowest,ipx,ipy,ipz,mglobal):
        #primative quantities read directly from file
        self.mx = mx
        self.my = my
        self.mz = mz
        self.mvar = mvar
        self.maux = maux
        self.precision = precision
        self.nghostx = nghostx
        self.nghosty = nghosty
        self.nghostz = nghostz
        self.nprocx = nprocx
        self.nprocy = nprocy
        self.nprocz = nprocz
        self.iprocz_slowest = iprocz_slowest
        self.ipx = ipx
        self.ipy = ipy
        self.ipz = ipz
        self.mglobal = mglobal	
        
        #derived quantities
        self.nx = mx - (2 * nghostx)
        self.ny = my - (2 * nghosty)
        self.nz = mz - (2 * nghostz)
        self.mw = mx * my * mz      
        self.l1 = nghostx           
        self.l2 = mx-nghostx-1      
        self.m1 = nghosty           
        self.m2 = my-nghosty-1      
        self.n1 = nghostz           
        self.n2 = mz-nghostz-1      
        if (ipx == ipy == ipz == -1):
            # global
            self.nxgrid = self.nx
            self.nygrid = self.ny
            self.nzgrid = self.nz
            self.mxgrid = self.nxgrid + (2 * nghostx)
            self.mygrid = self.nygrid + (2 * nghosty)
            self.mzgrid = self.nzgrid + (2 * nghostz)
        else:
            # local
            self.nxgrid = self.nygrid = self.nzgrid = 0
            self.mxgrid = self.mygrid = self.mzgrid = 0

        
def read_dim(datadir='data',proc=-1,down=False):
    """
    read the dim.dat file. if proc is -1, then read the 'global'
    dimensions. if proc is >=0, then read the dim.dat in the
    corresponding processor directory.
    """
    if down:
        dimdat='dim_down.dat'
    else:
        dimdat='dim.dat'
    if (proc < 0 ):
        filename = datadir+'/'+dimdat # global box dimensions
    else:
        filename = datadir+'/proc'+str(proc)+'/'+dimdat # local proc. dimensions

    try:
        filename = os.path.expanduser(filename)
        file = open(filename,"r")
    except IOError:
        #print "File",filename,"could not be opened." # Python 2
        print("File " + filename + " could not be opened.")
        return -1
    else:
        lines = file.readlines()
        file.close()

    if len(lines[0].split()) == 6:
        mx,my,mz,mvar,maux,mglobal = tuple(map(int,lines[0].split()))
    else:
        mx,my,mz,mvar,maux = tuple(map(int,lines[0].split()))
        mglobal = 0

    precision = lines[1].strip("\n")
    nghostx,nghosty,nghostz = tuple(map(int,lines[2].split()))
    if (proc < 0):
        #global
        nprocx,nprocy,nprocz,iprocz_slowest = tuple(map(int,lines[3].split()))
        ipx = ipy = ipz = -1
    else:
        #local processor
        ipx,ipy,ipz = tuple(map(int,lines[3].split()))
        nprocx = nprocy = nprocz = iprocz_slowest = -1

    dim = PcDim(mx,my,mz,mvar,maux,precision,nghostx,nghosty,nghostz,nprocx,nprocy,nprocz,iprocz_slowest,ipx,ipy,ipz,mglobal)

    return dim
