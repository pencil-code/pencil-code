# $Id$
#
# read VAR files. based on the read_var.pro IDL script.
#
# NB: the f array returned is C-ordered: f[nvar,nz,ny,nx]
#     NOT Fortran as in Pencil (& IDL):  f[nx,ny,nz,nvar]
#
# Author: J. Oishi (joishi@amnh.org). 
# 
# modify var.py in a more object-oriented way
# 03/08 : T. Gastine (tgastine@ast.obs-mip.fr)
# 
import numpy as N
from npfile import npfile
import os
import string
import sys
from param import read_param 
from index import read_index
from dim import read_dim 
from pencil.math.derivatives import curl, curl2

class read_var:
# !!!  The file format written by output() (and used, e.g. in var.dat)
# !!!  consists of the followinig Fortran records:
# !!!    1. data(mx,my,mz,nvar)
# !!!    2. t(1), x(mx), y(my), z(mz), dx(1), dy(1), dz(1), deltay(1)
# !!!  Here nvar denotes the number of slots, i.e. 1 for one scalar field, 3
# !!!  for one vector field, 8 for var.dat in the case of MHD with entropy.
# but, deltay(1) is only there if lshear is on! need to know parameters...

    def __init__(self, varfile='', datadir='data/', proc=-1, ivar=-1,
                 quiet=False, trimall=False, format='native',
                 param=None, dim=None, index=None, run2D=False, 
                 magic=None, setup=None):
        """
        Description:
        -----------
        read VAR files from pencil code. if proc < 0, then load all data
        and assemble. otherwise, load VAR file from specified processor.
        
        format -- one of (['native', 'n'], ['ieee-le', 'l'], 
        ['ieee-be', 'B']) for byte-ordering
        
        Params:
        ------
            varfile=''
            datadir='data/'
            proc=-1
            ivar=-1
            quiet=False
            trimall=False
            format='native'
            param=None
            dim=None
            index=None
            run2D=False
        """
        if (setup is not None):
            datadir=os.path.expanduser(setup.datadir)
            dim=setup.dim
            param=setup.param
            index=setup.index
            run2D=setup.run2D
        else:
            datadir = os.path.expanduser(datadir)
            if dim is None:
                dim = read_dim(datadir,proc) 
            if param is None:
                param = read_param(datadir=datadir, quiet=quiet)
            if index is None:
                index = read_index(datadir=datadir)

        if dim.precision == 'D':
            precision = 'd'
        else:
            precision = 'f'
        
        if param.lwrite_aux:
            totalvars = dim.mvar+dim.maux
        else:
            totalvars = dim.mvar

        # read index.pro to get positions and "names"
        # of variables in f(mx,my,mz,nvar)
        # Thomas: seems useless now ?
        #exec(index) # this loads the indicies.

        if (not varfile):
            if ivar < 0:
                varfile='var.dat'
            else:
                varfile='VAR'+str(ivar)
        
        if proc < 0:
            procdirs = filter(lambda s:s.startswith('proc'), 
                              os.listdir(datadir))
        else:
            procdirs = ['proc'+str(proc)]

        #global array
        if (not run2D):
            f = N.zeros((totalvars, dim.mz, dim.my, dim.mx),
                        dtype=precision)
        else:
            if dim.ny == 1:
                f = N.zeros((totalvars, dim.mz, dim.mx), dtype=precision)
            else:
                f = N.zeros((totalvars, dim.my, dim.mx), dtype=precision)
        x = N.zeros(dim.mx, dtype=precision)
        y = N.zeros(dim.my, dtype=precision)
        z = N.zeros(dim.mz, dtype=precision)
        for directory in procdirs:

            proc = int(directory[4:])
            procdim = read_dim(datadir, proc)
            if (not quiet):
                print "reading data from processor %i of %i ..." \
                      % (proc, len(procdirs))

            mxloc = procdim.mx
            myloc = procdim.my
            mzloc = procdim.mz

            #read data
            filename = os.path.join(datadir,directory,varfile)
            infile = npfile(filename, endian=format)
            if (not run2D):
                f_loc = infile.fort_read(precision, 
                               shape=(-1, mzloc, myloc, mxloc))
            else:
                if dim.ny == 1:
                    f_loc = infile.fort_read(precision,
                                   shape=(-1, mzloc, mxloc))
                else:
                    f_loc = infile.fort_read(precision,
                                   shape=(-1, myloc, mxloc))
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

            if len(procdirs) > 1:
                # calculate where the local processor will go in 
                # the global array

                #
                #  Don't overwrite ghost zones of processor to the left (and
                #  accordingly in y and z direction--makes a difference on the
                #  diagonals)
                #
                # recall that in NumPy, slicing is NON-INCLUSIVE on the right end
                # ie, x[0:4] will slice all of a 4-digit array, not produce
                # an error like in idl.
                
                if procdim.ipx == 0:
                    i0x = 0
                    i1x = i0x+procdim.mx
                    i0xloc = 0 
                    i1xloc = procdim.mx
                else:
                    i0x = procdim.ipx*procdim.nx+procdim.nghostx 
                    i1x = i0x+procdim.mx-procdim.nghostx
                    i0xloc = procdim.nghostx
                    i1xloc = procdim.mx
                    
                if procdim.ipy == 0:
                    i0y = 0
                    i1y = i0y+procdim.my
                    i0yloc = 0 
                    i1yloc = procdim.my
                else:
                    i0y = procdim.ipy*procdim.ny+procdim.nghosty 
                    i1y = i0y+procdim.my-procdim.nghosty
                    i0yloc = procdim.nghosty 
                    i1yloc = procdim.my
                        
                if procdim.ipz == 0:
                    i0z = 0
                    i1z = i0z+procdim.mz
                    i0zloc = 0 
                    i1zloc = procdim.mz
                else:
                    i0z = procdim.ipz*procdim.nz+procdim.nghostz 
                    i1z = i0z+procdim.mz-procdim.nghostz
                    i0zloc = procdim.nghostz 
                    i1zloc = procdim.mz

                x[i0x:i1x] = x_loc[i0xloc:i1xloc]
                y[i0y:i1y] = y_loc[i0yloc:i1yloc]
                z[i0z:i1z] = z_loc[i0zloc:i1zloc]
                
                if (not run2D):
                    f[:, i0z:i1z, i0y:i1y, i0x:i1x] = \
                        f_loc[:, i0zloc:i1zloc, i0yloc:i1yloc, i0xloc:i1xloc]
                else:
                    if dim.ny == 1:
                        f[:, i0z:i1z, i0x:i1x] = \
                              f_loc[:, i0zloc:i1zloc, i0xloc:i1xloc]
                    else:
                        f[:, i0y:i1y, i0x:i1x] = \
                              f_loc[:, i0yloc:i1yloc, i0xloc:i1xloc]
            else:
                f = f_loc
                x = x_loc
                y = y_loc
                z = z_loc
            #endif MPI run
        #endfor directories loop

        if (magic is not None):
            if ('bb' in magic):
                # compute the magnetic field before doing trimall
                aa = f[index['ax']-1:index['az'],...]
                self.bb = curl(aa,dx,dy,dz,run2D=param.lwrite_2d)
                if (trimall): self.bb=self.bb[:, dim.n1:dim.n2+1, 
                dim.m1:dim.m2+1, dim.l1:dim.l2+1]
            if ('jj' in magic):
                # compute the electric current field before doing trimall
                aa = f[index['ax']-1:index['az'],...]
                self.jj = curl2(aa,dx,dy,dz)
                if (trimall): self.jj=self.jj[:, dim.n1:dim.n2+1, 
                dim.m1:dim.m2+1, dim.l1:dim.l2+1]
            if ('vort' in magic):
                # compute the vorticity field before doing trimall
                uu = f[index['ux']-1:index['uz'],...]
                self.vort = curl(uu,dx,dy,dz,run2D=param.lwrite_2d)
                if (trimall):
                    if (param.lwrite_2d):
                        if (dim.nz == 1):
                            self.vort=self.vort[:, dim.m1:dim.m2+1, 
                            dim.l1:dim.l2+1]
                        else:
                            self.vort=self.vort[:, dim.n1:dim.n2+1, 
                            dim.l1:dim.l2+1]
                    else:
                        self.vort=self.vort[:, dim.n1:dim.n2+1, 
                        dim.m1:dim.m2+1, dim.l1:dim.l2+1]

        # trim the ghost zones of the global f-array if asked
        if trimall:
            self.x = x[dim.l1:dim.l2+1]
            self.y = y[dim.m1:dim.m2+1]
            self.z = z[dim.n1:dim.n2+1]
            if (not run2D):
                self.f = f[:, dim.n1:dim.n2+1, dim.m1:dim.m2+1, dim.l1:dim.l2+1]
            else:
               if dim.ny == 1:
                   self.f = f[:, dim.n1:dim.n2+1, dim.l1:dim.l2+1]
               else:
                   self.f = f[:, dim.m1:dim.m2+1, dim.l1:dim.l2+1]
        else:
            self.x = x
            self.y = y
            self.z = z
            self.f = f
            self.l1 = dim.l1
            self.l2 = dim.l2+1
            self.m1 = dim.m1
            self.m2 = dim.m2+1
            self.n1 = dim.n1
            self.n2 = dim.n2+1
            
        # Assign an attribute to self for each variable defined in
        # 'data/index.pro' so that e.g. self.ux is the x-velocity
        for key,value in index.items():
#          print key,value
          if key != 'global_gg':
            setattr(self,key,self.f[value-1,...])
        # Special treatment for vector quantities
        if index.has_key('uu'):
          self.uu = self.f[index['ux']-1:index['uz'],...]
        if index.has_key('aa'):
          self.aa = self.f[index['ax']-1:index['az'],...]

        self.t = t
        self.dx = dx
        self.dy = dy
        self.dz = dz
        if param.lshear:
            self.deltay = deltay

        # do the rest of magic after the trimall (i.e. no additional curl...)
        self.magic = magic
        if self.magic is not None:
            self.magicAttributes(param)

    def magicAttributes(self, param):
        for field in self.magic:
            if field == 'rho':
                if hasattr(self, 'lnrho'):
                    setattr(self, 'rho', N.exp(self.lnrho))
                else:
                    sys.exit("pb in magic!")

            if (field == 'tt' and not hasattr(self,'tt')):
                if hasattr(self, 'lnTT'):
                    tt = N.exp(self.lnTT)
                    setattr(self, 'tt', tt)
                else:
                    if (hasattr(self, 'lnrho') and hasattr(self, 'ss')):
                        cp = param.cp
                        gamma = param.gamma
                        cs20 = param.cs0**2
                        lnrho0 = N.log(param.rho0)
                        lnTT0 = N.log(cs20/(cp*(gamma-1.)))
                        lnTT = lnTT0+gamma/cp*self.ss+(gamma-1.)* \
                               (self.lnrho-lnrho0)
                        setattr(self, 'tt', N.exp(lnTT))
                    else:
                        sys.exit("pb in magic!")

            if (field == 'ss' and not hasattr(self, 'ss')):
                cp = param.cp
                gamma = param.gamma
                cs20 = param.cs0**2
                lnrho0 = N.log(param.rho0)
                lnTT0 = N.log(cs20/(cp*(gamma-1.)))
                if hasattr(self, 'lnTT'):
                    self.ss = cp/gamma*(self.lnTT-lnTT0- \
                              (gamma-1.)*(self.lnrho-lnrho0))
                elif hasattr(self, 'tt'):
                    self.ss = cp/gamma*(N.log(self.tt)- \
                              lnTT0-(gamma-1.)*(self.lnrho-lnrho0))
                else:
                    sys.exit("pb in magic!")

