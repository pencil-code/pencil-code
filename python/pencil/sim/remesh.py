# remesh.py
#
# Read VAR files from existing mature simulation located in the source_path
# directory to be remeshed into newly initialised target simulation of new
# grid and processor layout located in the target_path directory.
# 
# First set up new run compile and start it.
# To interpolate the old data into new shape call interp_var with string
# arguments specifying source_path and target_path directories.
# Then use the resulting fnew array and target_path as arguments for
# distribute_fort to overwrite the var.dat in the process tree
# 
#
# NB: to include particles or additional farray variables the scripts will
#     need to be modified
#     the interpolation does not handle switching the system of coordinates
#     for periodic boundaries the arrays are stretched or shrunk to maintain
#     continuity on the boundaries. This is the default.
#     Not implemented: clipping the data or extending e.g. repeat the array in
#     the periodic directions
#     Not implemented reading or writing farray in hdf5 format
#
# Author: F. Gent (fred.gent.ncl@gmail.com).
#
import numpy as np
from scipy.interpolate import interp1d
from ..math.derivatives import grad
from ..io import open_h5, group_h5, dataset_h5

def local_remesh(var,
           xsrc, ysrc, zsrc, xdst, ydst, zdst, quiet=True
          ):
    """
    Keyword arguments:

    *var*:
      snapshot scalar numpy array of shape [mz,my,mx].

    *xsrc, ysrc, zsrc:
      hdf5 grid object from source simulation.

    *xdst, ydst, zdst*:
      hdf5 grid object for destination simulation.

    *quiet*
      Flag for switching of output.

    """
    tmp = var.copy()
    if not xsrc == xdst:
        interpx = interp1d(xsrc,tmp,axis=-1)
        tmp = interpx(xdst)
    if not ysrc == ydst:
        interpy = interp1d(ysrc,tmp,axis=-2)
        tmp=interpy(ydst)
    if not zsrc == zdst:
        interpz = interp1d(zsrc,tmp,axis=-3)
        tmp = interpz(zdst)

    return tmp

def get_dstgrid(srch5, srcpar, dsth5,
                multxyz=[2,2,2], fracxyz=[1,1,1], srcghost=3, dstghost=3,
                dstprecision=[b'D'], lsymmetric=True
               ):
    """
    Keyword arguments:

    *srcgrid*:
      hdf5 object from source simulation.

    *srcpar*:
      simulation param dictionary object from source simulation.

    *dsth5*:
      hdf5 object for destination simulation data

    *srcghost*
      Number of ghost zones from the source order of accuracy (mx-nx)/2

    *dstghost*
      Number of ghost zones for the destination order of accuracy (mx-nx)/2

    *quiet*
      Flag for switching of output.

    *dstprecision*
      Precision used in destination simulation. Default double.

    """
    # check prime factorization of the result and display for proc options
    # if using fft check options for grid and cpu layout
    # handling non-equidistant grids tba
    if dstprecision[0] == b'D':
        dtype = np.float64
    elif dstprecision[0] == b'S':
        dtype = np.float32
    else:
        print('precision '+dstprecision+' not valid')
        return 1

    srcsets=srch5['settings']
    if not dsth5.__contains__('settings'):
        sets = dsth5.create_group('settings')
        for key in srcsets.keys():
            sets.create_dataset(key, data=srcsets[key][()])
    else:
        sets = dsth5['settings']
        for key in srcsets.keys():
            if not sets.__contains__(key):
                sets.create_dataset(key, data=srcsets[key][()])
    #update grid dimensions
    sets['nx'][()] = int(srcsets['nx'][()]*multxyz[0]/fracxyz[0])
    sets['mx'][()] = sets['nx'][()] + 2*dstghost
    sets['ny'][()] = int(srcsets['ny'][()]*multxyz[1]/fracxyz[1])
    sets['my'][()] = sets['ny'][()] + 2*dstghost
    sets['nz'][()] = int(srcsets['nz'][()]*multxyz[2]/fracxyz[2])
    sets['mz'][()] = sets['nz'][()] + 2*dstghost
    sets['l1'][()] = dstghost
    sets['l2'][()] = sets['mx'][()] - 1-dstghost
    sets['m1'][()] = dstghost
    sets['m2'][()] = sets['my'][()] - 1-dstghost
    sets['n1'][()] = dstghost
    sets['n2'][()] = sets['mz'][()] - 1-dstghost
    #calculate the new grid arrays
    srcgrid=srch5['grid']
    if not dsth5.__contains__('grid'):
        grid = dsth5.create_group('grid')
        for key in srcgrid.keys():
            grid.create_dataset(key, data=srcgrid[key][()], dtype=dtype)
    else:
        grid = dsth5['grid']
        for key in srcgrid.keys():
            if not grid.__contains__(key):
                grid.create_dataset(key, data=srcgrid[key][()])
    for ii,mm in [[0,'mx'],[1,'my'],[2,'mz']]:
        if not srcpar['lequidist'][ii]:
            print('get_dstgrid WARNING: non-equidistant grid not implemented\n',
                  'continuing with equidistant grid. Please correct omission')
        if not sets[mm][()] == srcsets[mm][()]:
            #assuming for now par.lxyz is the same
            mstr = mm[1]
            grid['d'+mstr][()] = dtype(srcgrid['d'+mstr][()]/
                             multxyz[ii]*fracxyz[ii])
            grid.__delitem__(mstr)
            grid.create_dataset(mstr, (sets[mm][()],), dtype=dtype)
            if not srcpar['lperi'][ii] and lsymmetric:
                grid[mstr][dstghost:-dstghost] = np.linspace(
                                srcgrid[mstr][srcghost]-grid['d'+mstr][()],
                                srcgrid[mstr][-srcghost-1],
                                sets['n'+mstr][()],dtype=dtype
                                ) + dtype(0.5*grid['d'+mstr][()])
            else:
                grid[mstr][dstghost:-dstghost] = np.linspace(
                                srcgrid[mstr][srcghost]-grid['d'+mstr][()],
                                srcgrid[mstr][-srcghost-1],
                                sets['n'+mstr][()],dtype=dtype
                                )
            for jj in range(0,dstghost):
                grid[mstr][jj] = grid[mstr][dstghost] -\
                                    (dstghost-jj)*grid['d'+mstr][()]
                grid[mstr][jj-dstghost] = grid[mstr][-dstghost-1] +\
                                    (jj+1)*grid['d'+mstr][()]
            if not srcpar['lperi'][ii]:
                grid['L'+mstr][()] = srcgrid['L'+mstr][()] + grid['d'+mstr][()]
                grid['O'+mstr][()] = srcgrid['O'+mstr][()] -\
                                                         0.5*grid['d'+mstr][()]
            grid['d'+mstr+'_1'][()] = 1./grid['d'+mstr][()]
            grid['d'+mstr+'_tilde'][()] = np.gradient(grid['d'+mstr+'_1'][()],
                                               grid['d'+mstr][()])
 
def src2dst_remesh(src, dst,
                   h5in='var.h5', h5out='var.h5',
                   multxyz=[2,2,2], fracxyz=[1,1,1], srcghost=3, dstghost=3,
                   srcdatadir='data/allprocs', dstdatadir='data/allprocs',
                   dstprecision=[b'D'], lsymmetric=True, quiet=True
                  ):
    import h5py
    from .. import read
    from os.path import join
    import os
    from ..io import mkdir	
    from . import is_sim_dir, simulation

    if dstprecision[0] == b'D':
        dtype = np.float64
    elif dstprecision[0] == b'S':
        dtype = np.float32
    else:
        print('precision '+dstprecision+' not valid')
        return 1

     
    if is_sim_dir(src):
        srcsim = simulation(src,quiet=quiet)
    else: 
        print('src2dst_remesh ERROR: src"'+src+
              '" is not a valid simulation path')
        return 1
    if is_sim_dir(dst):
        dstsim = simulation(dst,quiet=quiet)
    else:
        dstname = str.split(dst,'/')[-1]
        dstpath = str.strip(dst,dstname)
        if len(dstpath) == 0:
            dstpath = str.strip(srcsim.path,srcsim.name)
        dstsim = srcsim.copy(path_root=dstpath, name=dstname)
    with open_h5(join(srcsim.path,srcdatadir,h5in),mode='r') as srch5:
        with open_h5(join(dstsim.path,dstdatadir,h5out), mode='w') as dsth5:
            get_dstgrid(srch5, srcsim.param, dsth5,
                        multxyz=multxyz, fracxyz=fracxyz, srcghost=srcghost,
                        dstghost=dstghost,
                        dstprecision=dstprecision, lsymmetric=lsymmetric )
            group = group_h5(dsth5, 'unit', mode='w')
            for key in srch5['unit'].keys():
                if type(srch5['unit'][key][()]) == np.float64 or\
                    type(srch5['unit'][key][()]) == np.float32:
                    dset = dataset_h5(group, key, mode='w',
                                      data=srch5['unit'][key][()],
                                      overwrite=True, dtype=dtype)
                else:
                    dset = dataset_h5(group, key, mode='w',
                                      data=srch5['unit'][key][()],
                                      overwrite=True)
            if 'persist' in srch5.keys():
                group = group_h5(dsth5, 'persist', mode='w')
                for key in srch5['persist'].keys():
                    if type(srch5['persist'][key][()]) == np.float64 or\
                        type(srch5['persist'][key][()]) == np.float32:
                        dset = dataset_h5(group, key, mode='w',
                                          data=srch5['persist'][key][()],
                                          overwrite=True, dtype=dtype)
                    else:
                        dset = dataset_h5(group, key, mode='w',
                                          data=srch5['persist'][key][()],
                                          overwrite=True)
            dst = dataset_h5(dsth5, 'time', mode='w',
                             data=srch5['time'][()], dtype=dtype)
            srcgrid = srch5['grid']
            dstgrid = dsth5['grid']
            group = group_h5(dsth5, 'data', mode='w')
            for key in srch5['data'].keys():
                var = local_remesh(srch5['data'][key][()], 
                                   srcgrid['x'], srcgrid['y'], srcgrid['z'],
                                   dstgrid['x'], dstgrid['y'], dstgrid['z'],
                                   quiet=quiet )
                dst = dataset_h5(dsth5['data'], key, mode='w', data=var,
                                 overwrite=True, dtype=dtype)
# remains to copy other files and edit param files

