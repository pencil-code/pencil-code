# $Id: var.py,v 1.3 2008-03-07 16:41:15 dintrans Exp $
#
# read VAR files. based on the read_var.pro IDL script.
#
# NB: the f array returned is C-ordered: f[nvar,nz,ny,nx]
#     NOT Fortran as in Pencil (& IDL):  f[nx,ny,nz,nvar]
#
# Author: J. Oishi (joishi@amnh.org). 
# 
# 
import numpy as N
from npfile import npfile
import os
from param import read_param 
from dim import read_dim 

class VarFileData:
    pass

# !!!  The file format written by output() (and used, e.g. in var.dat)
# !!!  consists of the followinig Fortran records:
# !!!    1. data(mx,my,mz,nvar)
# !!!    2. t(1), x(mx), y(my), z(mz), dx(1), dy(1), dz(1), deltay(1)
# !!!  Here nvar denotes the number of slots, i.e. 1 for one scalar field, 3
# !!!  for one vector field, 8 for var.dat in the case of MHD with entropy.

# but, deltay(1) is only there if lshear is on! need to know parameters...
def read_var(varfile='',datadir='data/',proc=-1,ivar=-1,quiet=False,trimall=False,format='native',param=None,dim=None, run2D=False):
    """
    read VAR files from pencil code. if proc < 0, then load all data
    and assemble. otherwise, load VAR file from specified processor.
    
    format -- one of (['native', 'n'], ['ieee-le', 'l'], ['ieee-be', 'B']) for byte-ordering
    """
    datadir = os.path.expanduser(datadir)
    if (dim == None):
      dim = read_dim(datadir,proc) 
    if (param == None):
      param = read_param(datadir)
    if dim.precision == 'D':
        precision = 'd'
    else:
        precision = 'f'
    
    if (param.lwrite_aux):
        totalvars = dim.mvar+dim.maux
    else:
        totalvars = dim.mvar

    # read index.pro to get positions and "names"
    # of variables in f(mx,my,mz,nvar)
    index = read_index(datadir)
    exec(index) # this loads the indicies.

    if (not varfile):
        if (ivar < 0):
            varfile='var.dat'
        else:
            varfile='VAR'+str(ivar)
    
    if (proc < 0):
        procdirs = filter(lambda s:s.startswith('proc'),os.listdir(datadir))
    else:
        procdirs = ['proc'+str(proc)]

    #global array
    if (not run2D):
      f = N.zeros((totalvars,dim.mz,dim.my,dim.mx),dtype=precision)
    else :
      if (dim.ny==1):
        f = N.zeros((totalvars,dim.mz,dim.mx),dtype=precision)
      else:
        f = N.zeros((totalvars,dim.my,dim.mx),dtype=precision)
    x = N.zeros(dim.mx,dtype=precision)
    y = N.zeros(dim.my,dtype=precision)
    z = N.zeros(dim.mz,dtype=precision)
    for directory in procdirs:

        proc = int(directory[4:])
        procdim = read_dim(datadir,proc)
        if(not quiet):
            print "reading data from processor",proc,"of",len(procdirs),"..."

        mxloc = procdim.mx
        myloc = procdim.my
        mzloc = procdim.mz

        #read data
        filename = datadir+directory+'/'+varfile
        infile = npfile(filename,endian=format)
        if (not run2D):
          f_loc = infile.fort_read(precision,shape=(-1,mzloc,myloc,mxloc))
        else:
          if (dim.ny==1):
            f_loc = infile.fort_read(precision,shape=(-1,mzloc,mxloc))
          else:
            f_loc = infile.fort_read(precision,shape=(-1,myloc,mxloc))
        raw_etc = infile.fort_read(precision)
        infile.close()

        t = raw_etc[0]
        x_loc = raw_etc[1:mxloc+1]
        y_loc = raw_etc[mxloc+1:mxloc+myloc+1]
        z_loc = raw_etc[mxloc+myloc+1:mxloc+myloc+mzloc+1]
        if (param.lshear):
            shear_offset = 1
            deltay = raw_etc[-1]
        else:
            shear_offset = 0
            
        dx = raw_etc[-3-shear_offset]
        dy = raw_etc[-2-shear_offset]
        dz = raw_etc[-1-shear_offset]

        if (len(procdirs) > 1):
            # calculate where the local processor will go in the global array

            #
            #  Don't overwrite ghost zones of processor to the left (and
            #  accordingly in y and z direction--makes a difference on the
            #  diagonals)
            #
            # recall that in NumPy, slicing is NON-INCLUSIVE on the right end
            # ie, x[0:4] will slice all of a 4-digit array, not produce
            # an error like in idl.
            
            if (procdim.ipx == 0): 
                i0x=0
                i1x=i0x+procdim.mx
                i0xloc=0 
                i1xloc=procdim.mx
            else:
                i0x=procdim.ipx*procdim.nx+procdim.nghostx 
                i1x=i0x+procdim.mx-procdim.nghostx
                i0xloc=procdim.nghostx
                i1xloc=procdim.mx
                
            if (procdim.ipy == 0):
                i0y=0
                i1y=i0y+procdim.my
                i0yloc=0 
                i1yloc=procdim.my
            else:
                i0y=procdim.ipy*procdim.ny+procdim.nghosty 
                i1y=i0y+procdim.my-procdim.nghosty
                i0yloc=procdim.nghosty 
                i1yloc=procdim.my
                    
            if (procdim.ipz == 0):
                i0z=0
                i1z=i0z+procdim.mz
                i0zloc=0 
                i1zloc=procdim.mz
            else:
                i0z=procdim.ipz*procdim.nz+procdim.nghostz 
                i1z=i0z+procdim.mz-procdim.nghostz
                i0zloc=procdim.nghostz 
                i1zloc=procdim.mz

            x[i0x:i1x] = x_loc[i0xloc:i1xloc]
            y[i0y:i1y] = y_loc[i0yloc:i1yloc]
            z[i0z:i1z] = z_loc[i0zloc:i1zloc]
            
            if (not run2D):
              f[:,i0z:i1z,i0y:i1y,i0x:i1x] = f_loc[:,i0zloc:i1zloc,i0yloc:i1yloc,i0xloc:i1xloc]
            else:
              if (dim.ny==1):
                f[:,i0z:i1z,i0x:i1x] = f_loc[:,i0zloc:i1zloc,i0xloc:i1xloc]
              else:
                f[:,i0y:i1y,i0x:i1x] = f_loc[:,i0yloc:i1yloc,i0xloc:i1xloc]
        else:
            f = f_loc
            x = x_loc
            y = y_loc
            z = z_loc
        #endif MPI run
        
    #endfor directories loop

        
    # create a VarFileData class object to return
    var = VarFileData()

    # trim ghost zones if asked
    if (trimall):
        var.x = x[dim.l1:dim.l2+1]
        var.y = y[dim.m1:dim.m2+1]
        var.z = z[dim.n1:dim.n2+1]
        if (not run2D):
          var.f = f[:,dim.n1:dim.n2+1,dim.m1:dim.m2+1,dim.l1:dim.l2+1]
        else:
         if (dim.ny==1):
           var.f = f[:,dim.n1:dim.n2+1,dim.l1:dim.l2+1]
         else:
           var.f = f[:,dim.m1:dim.m2+1,dim.l1:dim.l2+1]
    else:
        var.x = x
        var.y = y
        var.z = z
        var.f = f
        var.l1 = dim.l1
        var.l2 = dim.l2+1
        var.m1 = dim.m1
        var.m2 = dim.m2+1
        var.n1 = dim.n1
        var.n2 = dim.n2+1
        

    var.t = t
    var.dx = dx
    var.dy = dy
    var.dz = dz
    if param.lshear:
        var.deltay = deltay

    return var
            

def read_index(datadir='data/'):
    if datadir.endswith('/'):
        datadir += '/'

    index = ''
    f = open(datadir+'index.pro')
    for line in f.readlines():
        clean = line.strip()

        if (not clean.endswith('0') and not clean.startswith('i_') and clean.startswith('i')):
            index += clean+'\n'

    return index
