
# ogvar.py
#
# Read the ogvar files.
# This routine is a simple adaptation of var.py used to read var files
# Necessary since dimension of var and ogvar are not the same
# NB: the f array returned is C-ordered: f[nvar, nz, ny, nx]
#     NOT Fortran as in Pencil (& IDL):  f[nx, ny, nz, nvar]
#
# Adapted from var.py by:
# J. Aarnes (jorgenaarnes@gmail.com)
"""
Contains the read class for the OGVAR file reading,
some simulation attributes and the data cube.
"""

def ogvar(*args, **kwargs):
    """
    Read OGVAR files from Pencil Code. If proc < 0, then load all data
    and assemble, otherwise load OGVAR file from specified processor.

    The file format written by output() (and used, e.g. in ogvar.dat)
    consists of the followinig Fortran records:
    1. data(mx, my, mz, nvar)
    2. t(1), x(mx), y(my), z(mz), dx(1), dy(1), dz(1), deltay(1)
    Here nvar denotes the number of slots, i.e. 1 for one scalar field, 3
    for one vector field, 8 for ogvar.dat in the case of MHD with entropy.
    but, deltay(1) is only there if lshear is on! need to know parameters.

    call signature:

    ogvar(ogvar_file='', datadir='data/', proc=-1, iogvar=-1,
          quiet=True, trimall=False,
          magic=None, sim=None, precision='f')

    Keyword arguments:
        ogvar_file: Name of the OGVAR file.
        sim:        Simulation sim object.
        magic:      Values to be computed from the data, e.g. B = curl(A).
        trimall:    Trim the data cube to exclude ghost zones.
        quiet:      Flag for switching off output.

        datadir:    Directory where the data is stored.
        proc:       Processor to be read. If -1 read all and assemble to one array.
        iogvar:     Index of the OGVAR file, if ogvar_file is not specified.
    """

    from ..sim import __Simulation__

    started = None

    for a in args:
        if type(a) == __Simulation__:
            started = a.started()
            break

    else:
        if 'sim' in kwargs.keys():
            started = kwargs['sim'].started()
        elif 'datadir' in kwargs.keys():
            from os.path import join, exists
            if exists(join(kwargs['datadir'], 'time_series.dat')): started = True
        #else:
        #    print('!! ERROR: No simulation of path specified..')

    if started == False:
        print('!! ERROR: Simulation has not jet started. There are not ogvar files.')
        return False

    ogvar_tmp = DataCube()
    ogvar_tmp.read(*args, **kwargs)
    return ogvar_tmp


class DataCube(object):
    """
    DataCube -- holds Pencil Code OGVAR file data.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.t = 0.
        self.dx = 1.
        self.dy = 1.
        self.dz = 1.


    def read(self, ogvar_file='', sim=None, datadir='data', proc=-1, iogvar=-1,
             quiet=True, trim_all=True, trimall=True,
             magic=None, transform=True, ogvarfile=''):
        """
        Read OGVAR files from pencil code. If proc < 0, then load all data
        and assemble. otherwise, load OGVAR file from specified processor.

        The file format written by output() (and used, e.g. in ogvar.dat)
        consists of the followinig Fortran records:
        1. data(mx, my, mz, nvar)
        2. t(1), x(mx), y(my), z(mz), dx(1), dy(1), dz(1), deltay(1)
        Here nvar denotes the number of slots, i.e. 1 for one scalar field, 3
        for one vector field, 8 for ogvar.dat in the case of MHD with entropy.
        but, deltay(1) is only there if lshear is on! need to know parameters.

        call signature:

        read(ogvar_file='', datadir='data/', proc=-1, iogvar=-1,
            quiet=True, trimall=False,
            magic=None, sim=None)

        Keyword arguments:
            ogvar_file/ogvarfile:
                        Name of the OGVAR file.
            sim:        Simulation sim object.
            magic:      Values to be computed from the data, e.g. B = curl(A).
            trimall:    Trim the data cube to exclude ghost zones.
            quiet:      Flag for switching off output.

            datadir:    Directory where the data is stored.
            proc:       Processor to be read. If -1 read all and assemble to one array.
            iogvar:     Index of the OGVAR file, if ogvar_file is not specified.
            transform:  By defauls, the velocities are given in polar cooridinates.
                        These are transformed to ux,uy,uz format (ur,uth,uz also in ogvar obj)

        """

        import numpy as np
        import os
        from scipy.io import FortranFile
        from pencilnew.math.derivatives import curl, curl2
        from pencilnew import read
        from ..sim import __Simulation__
        from pencilnew.math.transform import pospolar2cart, velpolar2cart

        ogdim = None; param = None; index = None

        if ogvarfile != '' and ogvar_file == '':
            ogvar_file = ogvarfile

        if type(ogvar_file) == __Simulation__:
            sim = ogvar_file
            ogvar_file = 'ogvar.dat'

        if sim is not None:
            datadir = os.path.expanduser(sim.datadir)
            ogdim = sim.ogdim
            param = read.param(datadir=sim.datadir, quiet=True)
            index = read.index(datadir=sim.datadir)
        else:
            datadir = os.path.expanduser(datadir)
            if ogdim is None:
                ogdim = read.ogdim(datadir, proc)
            if param is None:
                param = read.param(datadir=datadir, quiet=quiet)
            if index is None:
                index = read.index(datadir=datadir)

        run2D = param.lwrite_2d

        if ogdim.precision == 'D':
            precision = 'd'
        else:
            precision = 'f'

        #if param.lwrite_aux:
        #    total_ogvars = ogdim.mvar + ogdim.maux
        #else:
        #    total_ogvars = ogdim.mvar
        # Only mvar variables on overset grids at the moment
        total_ogvars = ogdim.mvar

        if not ogvar_file:
            if iogvar < 0:
                ogvar_file = 'ogvar.dat'
            else:
                ogvar_file = 'OGVAR' + str(iogvar)

        if proc < 0:
            proc_dirs = self.__natural_sort(filter(lambda s: s.startswith('proc'),
                                                   os.listdir(datadir)))
        else:
            proc_dirs = ['proc' + str(proc)]

        if trimall != False: trim_all = trimall

        # Set up the global array.
        if not run2D:
            f = np.zeros((total_ogvars, ogdim.mz, ogdim.my, ogdim.mx),
                         dtype=precision)
        else:
            if ogdim.ny == 1:
                f = np.zeros((total_ogvars, ogdim.mz, ogdim.mx), dtype=precision)
            else:
                f = np.zeros((total_ogvars, ogdim.my, ogdim.mx), dtype=precision)

        x = np.zeros(ogdim.mx, dtype=precision)
        y = np.zeros(ogdim.my, dtype=precision)
        z = np.zeros(ogdim.mz, dtype=precision)

        for directory in proc_dirs:
            proc = int(directory[4:])
            procdim = read.ogdim(datadir, proc)
            if not quiet:
                print("Reading data from processor {0} of {1} ...".format( \
                      proc, len(proc_dirs)))

            mxloc = procdim.mx
            myloc = procdim.my
            mzloc = procdim.mz

            # Read the data.
            file_name = os.path.join(datadir, directory, ogvar_file)
            infile = FortranFile(file_name)
            if not run2D:
                f_loc = infile.read_record(dtype=precision)
                f_loc = f_loc.reshape((-1, mzloc, myloc, mxloc))
            else:
                if ogdim.ny == 1:
                    f_loc = infile.read_record(dtype=precision)
                    f_loc = f_loc.reshape((-1, mzloc, mxloc))
                else:
                    f_loc = infile.read_record(dtype=precision)
                    f_loc = f_loc.reshape((-1, myloc, mxloc))
            raw_etc = infile.read_record(precision)
            infile.close()

            t = raw_etc[0]
            x_loc = raw_etc[1:mxloc+1]
            y_loc = raw_etc[mxloc+1:mxloc+myloc+1]
            z_loc = raw_etc[mxloc+myloc+1:mxloc+myloc+mzloc+1]
            if param.lshear:
                shear_offset = 1
                deltay = raw_etc[-1]
            else:
                shear_offset = 0

            dx = raw_etc[-3-shear_offset]
            dy = raw_etc[-2-shear_offset]
            dz = raw_etc[-1-shear_offset]

            if len(proc_dirs) > 1:
                # Calculate where the local processor will go in
                # the global array.
                #
                # Don't overwrite ghost zones of processor to the left (and
                # accordingly in y and z direction -- makes a difference on the
                # diagonals)
                #
                # Recall that in NumPy, slicing is NON-INCLUSIVE on the right end
                # ie, x[0:4] will slice all of a 4-digit array, not produce
                # an error like in idl.

                if procdim.ipx == 0:
                    i0x = 0
                    i1x = i0x + procdim.mx
                    i0xloc = 0
                    i1xloc = procdim.mx
                else:
                    i0x = procdim.ipx*procdim.nx + procdim.nghostx
                    i1x = i0x + procdim.mx - procdim.nghostx
                    i0xloc = procdim.nghostx
                    i1xloc = procdim.mx

                if procdim.ipy == 0:
                    i0y = 0
                    i1y = i0y + procdim.my
                    i0yloc = 0
                    i1yloc = procdim.my
                else:
                    i0y = procdim.ipy*procdim.ny + procdim.nghosty
                    i1y = i0y + procdim.my - procdim.nghosty
                    i0yloc = procdim.nghosty
                    i1yloc = procdim.my

                if procdim.ipz == 0:
                    i0z = 0
                    i1z = i0z+procdim.mz
                    i0zloc = 0
                    i1zloc = procdim.mz
                else:
                    i0z = procdim.ipz*procdim.nz + procdim.nghostz
                    i1z = i0z + procdim.mz - procdim.nghostz
                    i0zloc = procdim.nghostz
                    i1zloc = procdim.mz

                x[i0x:i1x] = x_loc[i0xloc:i1xloc]
                y[i0y:i1y] = y_loc[i0yloc:i1yloc]
                z[i0z:i1z] = z_loc[i0zloc:i1zloc]

                if not run2D:
                    f[:, i0z:i1z, i0y:i1y, i0x:i1x] = \
                        f_loc[:, i0zloc:i1zloc, i0yloc:i1yloc, i0xloc:i1xloc]
                else:
                    if ogdim.ny == 1:
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

        if magic is not None:
            if 'bb' in magic:
                # Compute the magnetic field before doing trim_all.
                aa = f[index.ax-1:index.az, ...]
                self.bb = curl(aa, dx, dy, dz, run2D=run2D)
                if trim_all:
                    self.bb = self.bb[:, ogdim.n1:ogdim.n2+1,
                                      ogdim.m1:ogdim.m2+1, ogdim.l1:ogdim.l2+1]
            if 'jj' in magic:
                # Compute the electric current field before doing trim_all.
                aa = f[index.ax-1:index.az, ...]
                self.jj = curl2(aa, dx, dy, dz)
                if trim_all:
                    self.jj = self.jj[:, ogdim.n1:ogdim.n2+1,
                                      ogdim.m1:ogdim.m2+1, ogdim.l1:ogdim.l2+1]
            if 'vort' in magic:
                # Compute the vorticity field before doing trim_all.
                uu = f[index.ux-1:index.uz, ...]
                self.vort = curl(uu, dx, dy, dz, run2D=run2D)
                if trim_all:
                    if run2D:
                        if ogdim.nz == 1:
                            self.vort = self.vort[:, ogdim.m1:ogdim.m2+1,
                                                  ogdim.l1:ogdim.l2+1]
                        else:
                            self.vort = self.vort[:, ogdim.n1:ogdim.n2+1,
                                                  ogdim.l1:ogdim.l2+1]
                    else:
                        self.vort = self.vort[:, ogdim.n1:ogdim.n2+1,
                                              ogdim.m1:ogdim.m2+1,
                                              ogdim.l1:ogdim.l2+1]

        # Trim the ghost zones of the global f-array if asked.
        if trim_all:
            self.x = x[ogdim.l1:ogdim.l2+1]
            self.y = y[ogdim.m1:ogdim.m2+1]
            self.z = z[ogdim.n1:ogdim.n2+1]
            if not run2D:
                self.f = f[:, ogdim.n1:ogdim.n2+1, ogdim.m1:ogdim.m2+1, ogdim.l1:ogdim.l2+1]
            else:
                if ogdim.ny == 1:
                    self.f = f[:, ogdim.n1:ogdim.n2+1, ogdim.l1:ogdim.l2+1]
                else:
                    self.f = f[:, ogdim.m1:ogdim.m2+1, ogdim.l1:ogdim.l2+1]
        else:
            self.x = x
            self.y = y
            self.z = z
            self.f = f
            self.l1 = ogdim.l1
            self.l2 = ogdim.l2 + 1
            self.m1 = ogdim.m1
            self.m2 = ogdim.m2 + 1
            self.n1 = ogdim.n1
            self.n2 = ogdim.n2 + 1

        # Assign an attribute to self for each variable defined in
        # 'data/index.pro' so that e.g. self.ux is the x-velocity
        for key in index.__dict__.keys():
            if key != 'global_gg' and key != 'keys':
                value = index.__dict__[key]
                setattr(self, key, self.f[value-1, ...])
        # Special treatment for vector quantities.
        if hasattr(index, 'uu'):
            self.uu = self.f[index.ux-1:index.uz, ...]
        if hasattr(index, 'aa'):
            self.aa = self.f[index.ax-1:index.az, ...]

        self.t = t
        self.dx = dx
        self.dy = dy
        self.dz = dz
        if param.lshear:
            self.deltay = deltay

        # Do the rest of magic after the trim_all (i.e. no additional curl...).
        self.magic = magic
        if self.magic is not None:
            self.magic_attributes(param)

        if transform: 
            if trim_all:
                zcoord=0
            else:
                zcoord=3

            self.r = self.x
            self.th = self.y
            self.x, self.y = pospolar2cart(self.r,self.th)
            self.ur = self.ux
            self.uth = self.uy
            self.ux, self.uy = velpolar2cart(self.ur,self.uth,self.r,self.th,zcoord)


    def __natural_sort(self, procs_list):
        """
        Sort array in a more natural way, e.g. 9OGVAR < 10OGVAR
        """

        import re

        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(procs_list, key=alphanum_key)


    def magic_attributes(self, param):
        """
        Compute some additional 'magic' quantities.
        """

        import numpy as np
        import sys

        for field in self.magic:
            if field == 'rho' and not hasattr(self, 'rho'):
                if hasattr(self, 'lnrho'):
                    setattr(self, 'rho', np.exp(self.lnrho))
                else:
                    sys.exit("pb in magic!")

            if field == 'tt' and not hasattr(self, 'tt'):
                if hasattr(self, 'lnTT'):
                    tt = np.exp(self.lnTT)
                    setattr(self, 'tt', tt)
                else:
                    if hasattr(self, 'ss'):
                        if hasattr(self, 'lnrho'):
                            lnrho = self.lnrho
                        elif hasattr(self, 'rho'):
                            lnrho = np.log(self.rho)
                        else:
                            sys.exit("pb in magic: missing rho or lnrho variable")
                        cp = param.cp
                        gamma = param.gamma
                        cs20 = param.cs0**2
                        lnrho0 = np.log(param.rho0)
                        lnTT0 = np.log(cs20/(cp*(gamma-1.)))
                        lnTT = lnTT0+gamma/cp*self.ss+(gamma-1.)* \
                               (lnrho-lnrho0)
                        setattr(self, 'tt', np.exp(lnTT))
                    else:
                        sys.exit("pb in magic!")

            if field == 'ss' and not hasattr(self, 'ss'):
                cp = param.cp
                gamma = param.gamma
                cs20 = param.cs0**2
                lnrho0 = np.log(param.rho0)
                lnTT0 = np.log(cs20/(cp*(gamma-1.)))
                if hasattr(self, 'lnTT'):
                    setattr(self, 'ss', cp/gamma*(self.lnTT-lnTT0- \
                            (gamma-1.)*(self.lnrho-lnrho0)))
                elif hasattr(self, 'tt'):
                    setattr(self, 'ss', cp/gamma*(np.log(self.tt)- \
                            lnTT0-(gamma-1.)*(self.lnrho-lnrho0)))
                else:
                    sys.exit("pb in magic!")

            if field == 'pp' and not hasattr(self, 'pp'):
                cp = param.cp
                gamma = param.gamma
                cv = cp/gamma
                if hasattr(self, 'lnrho'):
                    lnrho = self.lnrho
                elif hasattr(self, 'rho'):
                    lnrho = np.log(self.rho)
                else:
                    sys.exit("pb in magic: missing rho or lnrho variable")
                if hasattr(self, 'ss'):
                    setattr(self, 'pp', np.exp(gamma*(self.ss+lnrho)))
                elif hasattr(self, 'lntt'):
                    setattr(self, 'pp', (cp-cv)*np.exp(self.lntt+lnrho))
                elif hasattr(self, 'tt'):
                    setattr(self, 'pp', (cp-cv)*self.tt*np.exp(lnrho))
                else:
                    sys.exit("pb in magic!")
