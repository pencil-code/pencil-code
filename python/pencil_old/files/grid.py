# $Id$
#
# read grid
#
# Author: J. Oishi (joishi@amnh.org).
#
#
import numpy as N
from ..files.npfile import npfile
import os
from ..files.param import read_param
from ..files.dim import read_dim


def read_grid(*args, **kwargs):
    """Read grid from pencil code.

    If proc < 0, then load all data and assemble. otherwise, load grid
    from specified processor.
    """
    return Grid(*args, **kwargs)


class Grid(object):

    def __init__(self, datadir='data/', proc=-1, ivar=-1, quiet=False,
                 trim=False, format='native', param=None, down=False):
        """
        Read grid from pencil code. if proc < 0, then load all data
        and assemble. otherwise, load grid from specified processor.
        """
        datadir = os.path.expanduser(datadir)
        if param is None:
            param = read_param(datadir, quiet=quiet)
        dim = read_dim(datadir, proc, down=down)
        if dim.precision == 'D':
            precision = 'd'
        else:
            precision = 'f'

        if proc < 0:
            procdirs = list(filter(lambda s:s.startswith('proc'),
                                   os.listdir(datadir)))
            if (param.io_strategy == 'MPI-IO|collect' or param.io_strategy == 'collect'):
            	procdirs = ['allprocs']
        else:
            procdirs = ['proc'+str(proc)]

        #global array
        x = N.zeros(dim.mx, dtype=precision)
        y = N.zeros(dim.my, dtype=precision)
        z = N.zeros(dim.mz, dtype=precision)
        dx_1 = N.zeros(dim.mx, dtype=precision)
        dy_1 = N.zeros(dim.my, dtype=precision)
        dz_1 = N.zeros(dim.mz, dtype=precision)
        dx_tilde = N.zeros(dim.mx, dtype=precision)
        dy_tilde = N.zeros(dim.my, dtype=precision)
        dz_tilde = N.zeros(dim.mz, dtype=precision)

        for directory in procdirs:
            if directory != 'allprocs':
                proc = int(directory[4:])
            procdim = read_dim(datadir, proc, down=down)
            if not quiet:
                #print "reading data from processor %i of %i ..." \ # Python 2
                      #% (proc, len(procdirs)) # Python 2
                print("reading data from processor {0} of {1} ...".format(proc+1, len(procdirs)))

            mxloc = procdim.mx
            myloc = procdim.my
            mzloc = procdim.mz

            #read data
            if down:
                griddat='grid_down.dat'
            else:
                griddat='grid.dat'
            filename = os.path.join(datadir, directory, griddat)
            infile = npfile(filename, endian=format)
            grid_raw = infile.fort_read(precision)
            dx, dy, dz = tuple(infile.fort_read(precision))
            Lx, Ly, Lz = tuple(infile.fort_read(precision))
            dx_1_raw = infile.fort_read(precision)
            dx_tilde_raw = infile.fort_read(precision)
            infile.close()

            #reshape
            t = grid_raw[0]
            x_loc = grid_raw[1:mxloc+1]
            y_loc = grid_raw[mxloc+1:mxloc+myloc+1]
            z_loc = grid_raw[mxloc+myloc+1:mxloc+myloc+mzloc+1]
            dx_1_loc = dx_1_raw[0:mxloc]
            dy_1_loc = dx_1_raw[mxloc:mxloc+myloc]
            dz_1_loc = dx_1_raw[mxloc+myloc:mxloc+myloc+mzloc]
            dx_tilde_loc = dx_tilde_raw[0:mxloc]
            dy_tilde_loc = dx_tilde_raw[mxloc:mxloc+myloc]
            dz_tilde_loc = dx_tilde_raw[mxloc+myloc:mxloc+myloc+mzloc]

            if len(procdirs) >1:
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
                dx_1[i0x:i1x] = dx_1_loc[i0xloc:i1xloc]
                dy_1[i0y:i1y] = dy_1_loc[i0yloc:i1yloc]
                dz_1[i0z:i1z] = dz_1_loc[i0zloc:i1zloc]
                dx_tilde[i0x:i1x] = dx_tilde_loc[i0xloc:i1xloc]
                dy_tilde[i0y:i1y] = dy_tilde_loc[i0yloc:i1yloc]
                dz_tilde[i0z:i1z] = dz_tilde_loc[i0zloc:i1zloc]

            else:
                x = x_loc
                y = y_loc
                z = z_loc
                dx_1 = dx_1_loc
                dy_1 = dy_1_loc
                dz_1 = dz_1_loc
                dx_tilde = dx_tilde_loc
                dy_tilde = dy_tilde_loc
                dz_tilde = dz_tilde_loc
            #endif MPI run

        # end directories loop
        if trim:
            self.x = x[dim.l1:dim.l2+1]
            self.y = y[dim.m1:dim.m2+1]
            self.z = z[dim.n1:dim.n2+1]
            self.dx_1 = dx_1[dim.l1:dim.l2+1]
            self.dy_1 = dy_1[dim.m1:dim.m2+1]
            self.dx_1 = dz_1[dim.n1:dim.n2+1]
            self.dx_tilde = dx_tilde[dim.l1:dim.l2+1]
            self.dy_tilde = dy_tilde[dim.m1:dim.m2+1]
            self.dx_tilde = dz_tilde[dim.n1:dim.n2+1]
        else:
            self.x = x
            self.y = y
            self.z = z
            self.dx_1 = dx_1
            self.dy_1 = dy_1
            self.dx_1 = dz_1
            self.dx_tilde = dx_tilde
            self.dy_tilde = dy_tilde
            self.dx_tilde = dz_tilde

        self.t = t
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
