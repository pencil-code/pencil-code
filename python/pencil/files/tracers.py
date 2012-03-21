# $Id: tracers.py,v 1.4 2012-03-21 09:55:49 iomsn Exp $
#
# Reads the tracer files, composes a color map
# and reads the fixed point values.
#
# Author: Simon Candelaresi (iomsn@physto.se, iomsn1@googlemail.com).
# 
#

import struct
import array
import numpy as np
import os
import pencil as pc


def read_tracers(datadir = 'data/', filename = 'tracers.dat', zlim = []):
    """
    Reads the tracer files, composes a color map.

    call signature::
    
      tracers, mapping, t = read_tracers(filename = 'tracers.dat', datadir = 'data/', zlim)
    
    Reads from the tracer files and computes the color map according to
    A R Yeates and G Hornig 2011 J. Phys. A: Math. Theor. 44 265501
    doi:10.1088/1751-8113/44/26/265501.
    Returns the tracer values, the color mapping and the times of the snapshots.
    The color mapping can be plotted with:
    pc.animate_interactive(mapping[:,::-1,:,:], t, dimOrder = (2,1,0,3))
    
    Keyword arguments:
    
      *filename*:
        Name of the tracer file.

      *datadir*:
        Data directory.
        
      *zlim*:
        The upper limit for the field line mapping at which a field line is considered
        to have reached the upper boundary.
    """
    
    class data_struct():
        def __init__(self):
            self.xi = []
            self.yi = []
            self.xf = []
            self.yf = []
            self.zf = []
            self.l = []
            self.q = []

    data = []
    data = data_struct()

    # read the cpu structure
    dim = pc.read_dim(datadir = datadir)
    if (dim.nprocz > 1):
        print "error: number of cores in z-direction > 1"

    # read the parameters
    params = pc.read_param(datadir = datadir, quiet = True)
    
    # read the grid
    grid = pc.read_grid(datadir = datadir, quiet = True)

    # determine the file structure
    n_proc = dim.nprocx*dim.nprocy
    # sub sapling of the tracers
    trace_sub = params.trace_sub
    n_times = os.path.getsize("data/proc0/tracers.dat")/(4*(3 + 7*dim.nx*dim.ny*trace_sub**2/dim.nprocx/dim.nprocy))

    # prepare the output arrays
    tracers = np.zeros((dim.nx*trace_sub, dim.ny*trace_sub, n_times, 7))
    mapping = np.zeros((dim.nx*trace_sub, dim.ny*trace_sub, n_times, 3))

    # temporary arrays for one core                                 
    tracers_core = np.zeros((dim.nx*trace_sub/dim.nprocx, dim.ny*trace_sub/dim.nprocy, n_times, 7))
    mapping_core = np.zeros((dim.nx*trace_sub/dim.nprocx, dim.ny*trace_sub/dim.nprocy, n_times, 3))

    # set the upper z-limit to the domain boundary
    if zlim == []:
        zlim = grid.z[-dim.nghostz-1]
        
    # read the data from all cores
    for i in range(n_proc):
        # read the cpu structure
        dim_core = pc.read_dim(datadir = datadir, proc = i)
        stride = dim_core.nx*dim_core.ny*trace_sub**2    
        llen = 3 + 7*stride
        
        tracer_file = open(datadir+'proc{0}/'.format(i)+filename, 'rb')
        tmp = array.array('f')
        tmp.read(tracer_file, (3 + 7*dim_core.nx*dim_core.ny*trace_sub**2)*n_times)
        tracer_file.close()
        
        t = []
        
        for j in range(n_times):
            t.append(tmp[1+j*llen])
            data.xi = tmp[2+j*llen:2+stride + j*llen]
            data.yi = tmp[2+stride+j*llen:2+2*stride+j*llen]
            data.xf = tmp[2+2*stride+j*llen:2+3*stride+j*llen]
            data.yf = tmp[2+3*stride+j*llen:2+4*stride+j*llen]
            data.zf = tmp[2+4*stride+j*llen:2+5*stride+j*llen]
            data.l = tmp[2+5*stride+j*llen:2+6*stride+j*llen]
            data.q = tmp[2+6*stride+j*llen:2+7*stride+j*llen]

            # Squeeze the data into 2d array. This make the visualization much faster.
            for l in range(len(data.xi)):
                tracers_core[l%(dim_core.nx*trace_sub),l/(dim_core.nx*trace_sub),j,:] = \
                [data.xi[l], data.yi[l], data.xf[l], data.yf[l], data.zf[l], data.l[l], data.q[l]]
                if data.zf[l] >= zlim:
                    if (data.xi[l] - data.xf[l]) > 0:
                        if (data.yi[l] - data.yf[l]) > 0:
                            mapping_core[l%(dim_core.nx*trace_sub),l/(dim_core.nx*trace_sub),j,:] = [0,1,0]
                        else:
                            mapping_core[l%(dim_core.nx*trace_sub),l/(dim_core.nx*trace_sub),j,:] = [1,1,0]
                    else:
                        if (data.yi[l] - data.yf[l]) > 0:
                            mapping_core[l%(dim_core.nx*trace_sub),l/(dim_core.nx*trace_sub),j,:] = [0,0,1]
                        else:
                            mapping_core[l%(dim_core.nx*trace_sub),l/(dim_core.nx*trace_sub),j,:] = [1,0,0]
                else:
                    mapping_core[l%(dim_core.nx*trace_sub),l/(dim_core.nx*trace_sub),j,:] = [1,1,1]

            # copy single core data into total data arrays
            tracers[dim_core.ipx*dim_core.nx:(dim_core.ipx+1)*dim_core.nx, \
                    dim_core.ipy*dim_core.ny:(dim_core.ipy+1)*dim_core.ny,j,:] = \
                    tracers_core[:,:,j,:]
            mapping[dim_core.ipx*dim_core.nx:(dim_core.ipx+1)*dim_core.nx, \
                    dim_core.ipy*dim_core.ny:(dim_core.ipy+1)*dim_core.ny,j,:] = \
                    mapping_core[:,:,j,:]
                    
    return tracers, mapping, t



def read_fixed_points(datadir = 'data/', filename = 'fixed_points.dat'):
    """
    Reads the fixed points files.

    call signature::
    
      fixed = read_tracers(filename = 'tracers.dat', datadir = 'data/')
    
    Reads from the fixed points files.
    #Returns the fixed points positions with and the times of the snapshots.
    
    Keyword arguments:
    
      *filename*:
        Name of the tracer file.

      *datadir*:
        Data directory.
    """
    

    class data_struct():
        def __init__(self):
            self.t = []
            self.fidx = [] # number of fixed points at this time
            self.x = []
            self.y = []
            self.q = []

    # read the cpu structure
    dim = pc.read_dim()
    if (dim.nprocz > 1):
        print "error: number of cores in z-direction > 1"

    # determine the file structure
    n_proc = dim.nprocx*dim.nprocy
    
    data = []

    # total number of fixed points
    n_fixed = 0
    
    # read the data from all cores
    for i in range(n_proc):      
        fixed_file = open(datadir+'proc{0}/'.format(i)+filename, 'rb')
        tmp = fixed_file.read(4)
        
        data.append(data_struct())
        eof = 0
        if tmp == '':
            eof = 1
        while (eof == 0):
            data[i].t.append(struct.unpack("<ff", fixed_file.read(8))[0])
            n_fixed_core = int(struct.unpack("<fff", fixed_file.read(12))[1])
            n_fixed += n_fixed_core
            data[-1].fidx.append(n_fixed_core)

            x = list(np.zeros(n_fixed_core))
            y = list(np.zeros(n_fixed_core))
            for j in range(n_fixed_core):
                x[j] = struct.unpack("<ff", fixed_file.read(8))[1]
                y[j] = struct.unpack("<ff", fixed_file.read(8))[0]
            data[i].x.append(x)
            data[i].y.append(y)

            tmp = fixed_file.read(4)
            if tmp == '':
                eof = 1

        fixed_file.close()
        
    fixed = data_struct()
    for i in range(len(data[0].t)):
        fixed.t.append(data[0].t[i])
        x = []; y = []
        for proc in range(n_proc):
            x = x + data[proc].x[i]
            y = y + data[proc].y[i]
        fixed.x.append(x)
        fixed.y.append(y)
    
    fixed.t = np.array(fixed.t)
    fixed.x = np.array(fixed.x)
    fixed.y = np.array(fixed.y)
    
    return fixed



