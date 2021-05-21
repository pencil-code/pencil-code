# remesh.py
#
# 05-may-20
# Author: F. Gent (fred.gent.ncl@gmail.com).
#
""" Remesh mature simulation snapshot with [nx,ny,nz] dimensions onto new
    simulation with new grid dimensions and optionally alternate cpu layout
    copying the base simulation files, existing start output files.

    uses:
      local_remesh to apply the interpolation onto a variable array
      get_dstgrid to derive the new grid layout
      src2dst_remesh to create the new simulation object and files

"""

from fileinput import input
import numpy as np
from scipy.interpolate import interp1d
import subprocess as sub
from sys import stdout

from pencil.math.derivatives import grad
from ..io import open_h5, group_h5, dataset_h5


def local_remesh(
    var, xsrc, ysrc, zsrc, xdst, ydst, zdst, quiet=True
):
    """
    Call signature:

    local_remesh(var, xsrc, ysrc, zsrc, xdst, ydst, zdst, quiet=True)

    Keyword arguments:

    *var*:
      snapshot scalar numpy array of shape [mz,my,mx].

    *xsrc, ysrc, zsrc:
      grid x,y,z arrays from source simulation.

    *xdst, ydst, zdst*:
      grid x,y,z arrays for destination simulation.

    *quiet*:
      Flag for switching of output.

    """
    tmp = var.copy()
    if not quiet:
        print('x',tmp.shape, xsrc.min(), xsrc.max(), xdst.min(), xdst.max())
        print('x',tmp.shape, xsrc.shape, xdst.shape)
    if not xsrc.size == xdst.size:
        interpx = interp1d(xsrc,tmp,axis=-1)
        tmp = interpx(xdst)
    if not quiet:
        print('y',tmp.shape, ysrc.min(), ysrc.max(), ydst.min(), ydst.max())
        print('y',tmp.shape, ysrc.shape, ydst.shape)
    if not ysrc.size == ydst.size:
        interpy = interp1d(ysrc,tmp,axis=-2)
        tmp=interpy(ydst)
    if not quiet:
        print('z',tmp.shape, zsrc.min(), zsrc.max(), zdst.min(), zdst.max())
        print('z',tmp.shape, zsrc.shape, zdst.shape)
    if not zsrc.size == zdst.size:
        interpz = interp1d(zsrc,tmp,axis=-3)
        tmp = interpz(zdst)

    return tmp


def get_dstgrid(srch5, srcpar, dsth5, ncpus=[1,1,1],
                multxyz=[2,2,2], fracxyz=[1,1,1], srcghost=3, dstghost=3,
                dtype=np.float64, lsymmetric=True, quiet=True
               ):
    """
    Call signature:

    get_dstgrid(srch5, srcpar, dsth5, ncpus=[1,1,1], multxyz=[2,2,2],
               fracxyz=[1,1,1], srcghost=3, dstghost=3, dtype=np.float64,
               lsymmetric=True, quiet=True)

    Keyword arguments:

    *srch5*:
      hdf5 object from source simulation.

    *srcpar*:
      simulation param dictionary object from source simulation.

    *dsth5*:
      hdf5 object for destination simulation data.

    *ncpus*:
      array of nprocx, nprocy, and nprocz to apply for new simulation.

    *multxyz*:
      factors by which to multiply old sim dimensions yxz order.

    *fracxyz*:
      factors by which to divide old sim dimensions yxz order.

    *srcghost*:
      Number of ghost zones from the source order of accuracy (mx-nx)/2

    *dstghost*:
      Number of ghost zones for the destination order of accuracy (mx-nx)/2

    *dtype*:
      Precision used in destination simulation. Default double.

    *lsymmetric*:
      Option to make non-periodic grid symmetric about old sim centre.
      Otherwise the lower boundary is retained from old sim grid.

    *quiet*:
      Flag for switching of output.

    """
    # TBA
    # check prime factorization of the result and display for proc options
    # if using fft check options for grid and cpu layout
    # handling non-equidistant grids tba

    # copy settings from srcsim and revise with changes to dstsim var.h5
    srcsets=srch5['settings']
    sets = group_h5(dsth5, 'settings', status='a')
    for key in srcsets.keys():
        dset = dataset_h5(sets, key, data=srcsets[key][()], status='a')
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
    if not ncpus == [1,1,1]:
        sets['nprocx'][()] = ncpus[0]
        sets['nprocy'][()] = ncpus[1]
        sets['nprocz'][()] = ncpus[2]
    #copy the grid from the srcsim to dstsim var.h5 and grid.h5
    srcgrid=srch5['grid']
    grid = group_h5(dsth5, 'grid', status='a')
    for key in srcgrid.keys():
        dset = dataset_h5(grid, key, data=srcgrid[key][()], status='a')
    #replace grid data changed for dstsim
    for ii,mm in [[0,'mx'],[1,'my'],[2,'mz']]:
        if not srcpar['lequidist'][ii]:
            print('get_dstgrid WARNING: non-equidistant grid not implemented\n',
                  'continuing with equidistant grid.\n',
                  'Please implement non-equidistant grid options.')
        if not sets[mm][()] == srcsets[mm][()]:
            #assuming for now par.lxyz is the same
            mstr = mm[1]
            grid['d'+mstr][()] = dtype((srcgrid[mstr][-srcghost]-
                                srcgrid[mstr][srcghost])/(sets['n'+mstr][()]-1))
            grid.__delitem__(mstr)
            grid.create_dataset(mstr, (sets[mm][()],), dtype=dtype)
            print('grid 161:',mstr,srcgrid[mstr][srcghost],grid['d'+mstr][()],srcgrid[mstr][-srcghost-1][()],sets['n'+mstr][()])
            grid[mstr][dstghost:-dstghost] = np.linspace(
                            srcgrid[mstr][srcghost]-grid['d'+mstr][()],
                            srcgrid[mstr][-srcghost-1][()],
                            sets['n'+mstr][0],dtype=dtype
                            )
            if srcpar['lshift_origin'][ii] or lsymmetric:
                grid[mstr][dstghost:-dstghost] += dtype(0.5*grid['d'+mstr][()])
            elif srcpar['lshift_origin_lower'][ii]:
                grid[mstr][dstghost:-dstghost] -= dtype(0.5*grid['d'+mstr][()])
            for jj in range(0,dstghost):
                grid[mstr][jj] = grid[mstr][dstghost] -\
                                    (dstghost-jj)*grid['d'+mstr][()]
                grid[mstr][jj-dstghost] = grid[mstr][-dstghost-1] +\
                                    (jj+1)*grid['d'+mstr][()]
            if not srcpar['lperi'][ii]:
                grid['L'+mstr][()] = srcgrid['L'+mstr][()] + grid['d'+mstr][()]
                grid['O'+mstr][()] = srcgrid['O'+mstr][()] -\
                                                         0.5*grid['d'+mstr][()]
            grid.__delitem__('d'+mstr+'_1')
            grid.create_dataset('d'+mstr+'_1',
                               data=1./np.gradient(grid[mstr][()]), dtype=dtype)
            grid.__delitem__('d'+mstr+'_tilde')
            grid.create_dataset('d'+mstr+'_tilde',
                         data=np.gradient(grid['d'+mstr+'_1'][()]), dtype=dtype)

def src2dst_remesh(src, dst,
                   h5in='var.h5', h5out='var.h5',
                   multxyz=[2,2,2], fracxyz=[1,1,1], srcghost=3, dstghost=3,
                   srcdatadir='data/allprocs', dstdatadir='data/allprocs',
                   dstprecision=[b'D'], lsymmetric=True, quiet=True,
                   check_grid=True, OVERWRITE=False, optionals=True, nmin=32,
                   rename_submit_script=False, MBmin=5.0, ncpus=[1,1,1],
                   start_optionals=False, hostfile=None, submit_new=False,
                   chunksize=1000.0, lfs=False, MB=1, count=1,
                   size=1, rank=0, comm=None
                  ):
    """
    Call signature:

    src2dst_remesh(src, dst, h5in='var.h5', h5out='var.h5', multxyz=[2,2,2],
                   fracxyz=[1,1,1], srcghost=3, dstghost=3,
                   srcdatadir='data/allprocs', dstdatadir='data/allprocs',
                   dstprecision=[b'D'], lsymmetric=True, quiet=True,
                   check_grid=True, OVERWRITE=False, optionals=True, nmin=32,
                   rename_submit_script=False, MBmin=5.0, ncpus=[1,1,1],
                   start_optionals=False, hostfile=None, submit_new=False)

    Keyword arguments:

    *src*:
      string relative or absolute path to source simulation.

    *dst*:
      string relative or absolute path to destination simulation.

    *h5in*:
      source simulation data file to be copied and remeshed.

    *h5out*:
      destination simulation file to be written.

    *multxyz*:
      factors by which to multiply old sim dimensions yxz order.

    *fracxyz*:
      factors by which to divide old sim dimensions yxz order.

    *srcghost*:
      Number of ghost zones from the source order of accuracy (mx-nx)/2

    *dstghost*:
      Number of ghost zones for the destination order of accuracy (mx-nx)/2

    *srcdatadir*:
      path from source simulation directory to data.

    *dstdatadir*:
      path from destination simulation directory to data.

    *dstprecision*:
      floating point precision settings [b'S'] or [b'D'].

    *lsymmetric*:
      Option to make non-periodic grid symmetric about old sim centre.
      Otherwise the lower boundary is retained from old sim grid.

    *quiet*:
      Flag for switching of output.

    *check_grid*:
      Flag to run check on grid and cpu layout before executing remesh.

    *OVERWRITE*:
      Flag to overwrite existing simulation directory and filesin dst.

    *optionals*:
      Copy simulation files with True or specify list of names (string) for
      additional files from src sim directory.

    *nmin*:
      Minimum length along coordinate after splitting by proc.

    *rename_submit_script:
      Edit lines in submission files vcopied from src to dst.
      Not yet operational.

    *MBmin*:
      Minimum size in MB of data on a sinlge proc pf ncpus total processes.

    *ncpus*:
      array of nprocx, nprocy, and nprocz to apply for new simulation.

    *start_optionals*
      Copy simulation files output by start.x with True or specify list of
      names (string) for additional files from src sim data directory.

    *hostfile:
      Specify name of host config file argument in pc_build.
      Not yet operational.

    *submit_new*:
      Execute changes to submission files, compile and run simulation.
      Not yet operational.

    *chunksize*:
      Size in megabytes of snapshot variable before chunked remesh is used.

    *lfs*:
      Flag to set the striping for large file sizes to imporve IO efficiency.

    *MB*:
      Size of data to write contiguously before moving to new OST on lustre.

    *count*:
      Number of OSTs across which the data will be shared for IO operations.

    *comm*:
      MPI library calls

    *rank*:
      Integer ID of processor

    *size*:
      Number of MPI processes

    """
    import h5py
    import os
    from os.path import join, abspath
    import time

    from .. import read
    from ..io import mkdir
    from . import simulation
    from pencil.math import cpu_optimal
    from pencil import is_sim_dir

    start_time = time.time()
    print('started at {}'.format(time.ctime(start_time)))
    # set dtype from precision
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
        dstsim = srcsim.copy(path_root=dstpath, name=dstname, quiet=quiet,
                             OVERWRITE=OVERWRITE, optionals=optionals,
                             start_optionals=start_optionals,
                             rename_submit_script=rename_submit_script)
    print('opening src file and dst file on rank{}'.format(rank))
    with open_h5(join(srcsim.path,srcdatadir,h5in),'r',rank=rank,comm=comm) as srch5:
        with open_h5(join(dstsim.path,dstdatadir,h5out),'w',lfs=lfs,MB=MB,count=count,rank=rank,comm=comm) as dsth5:
            #apply settings and grid to dst h5 files
            get_dstgrid(srch5, srcsim.param, dsth5, ncpus=ncpus,
                        multxyz=multxyz, fracxyz=fracxyz, srcghost=srcghost,
                        dstghost=dstghost, dtype=dtype, lsymmetric=lsymmetric,
                        quiet=quiet)
            print('get_dstgrid completed on rank {}'.format(rank))
            #use settings to determine available proc dist then set ncpus
            factors = cpu_optimal(
                   dsth5['settings/nx'][0],
                   dsth5['settings/ny'][0],
                   dsth5['settings/nz'][0],
                   mvar=dsth5['settings/mvar'][0],
                   maux=dsth5['settings/maux'][0],
                   par=srcsim.param, nmin=nmin, MBmin=MBmin)
            print('remesh check grid: optional cpus upto min grid of'+
                  'nmin={}\n'.format(nmin)+
                  'cpu options {}\n'.format(factors)+
                  'new mesh: {}, {}, {}\n'.format(dsth5['settings/nx'][0],
                  dsth5['settings/ny'][0], dsth5['settings/nz'][0])+
                 'To execute remesh set "check_grid=False".')
            if ncpus == [1,1,1]:
                ncpus = [factors[1][0],factors[1][1],factors[1][2]]
                dsth5['settings/nprocx'][0] = ncpus[0]
                dsth5['settings/nprocy'][0] = ncpus[1]
                dsth5['settings/nprocz'][0] = ncpus[2]
            nprocs = ncpus[0]*ncpus[1]*ncpus[2]
            srcprocs = srch5['settings/nprocx'][0]*\
                       srch5['settings/nprocy'][0]*\
                       srch5['settings/nprocz'][0]
            if srcprocs > nprocs:
                print(
               '\n**********************************************************\n'+
               'remesh WARNING: {} procs reduced from {}.\n'.format(
               nprocs, srcprocs)+
               'Review multxyz {} and fracxyz {} for more\n'.format(
               multxyz,fracxyz)+
               'efficient parallel processing options.'+
               '\n**********************************************************\n')
            if check_grid:
                return 1
            group = group_h5(dsth5, 'unit', status='w')
            for key in srch5['unit'].keys():
                if type(srch5['unit'][key][()]) == np.float64 or\
                    type(srch5['unit'][key][()]) == np.float32:
                    dset = dataset_h5(group, key, status='w',
                                      data=srch5['unit'][key][()],
                                      overwrite=True, dtype=dtype)
                else:
                    dset = dataset_h5(group, key, status='w',
                                      data=srch5['unit'][key][()],
                                      overwrite=True)
            gridh5 = open_h5(join(dstsim.datadir,'grid.h5'), status='w')
            dsth5.copy('settings', gridh5)
            dsth5.copy('grid', gridh5)
            dsth5.copy('unit', gridh5)
            gridh5.close()
            if 'persist' in srch5.keys():
                group = group_h5(dsth5, 'persist', status='w')
                for key in srch5['persist'].keys():
                    tmp = np.zeros(nprocs)
                    tmp[:] = srch5['persist'][key][0]
                    if type(srch5['persist'][key][()]) == np.float64 or\
                                 type(srch5['persist'][key][()]) == np.float32:
                        dset = dataset_h5(group, key, status='w',
                                          data=tmp, overwrite=True, dtype=dtype)
                    else:
                        dset = dataset_h5(group, key, status='w',
                                          data=tmp, overwrite=True)
            dset = dataset_h5(dsth5, 'time', status='w',
                             data=srch5['time'][()], dtype=dtype)
            nx, ny, nz = dsth5['settings']['nx'][0],\
                         dsth5['settings']['ny'][0],\
                         dsth5['settings']['nz'][0]
            dstchunksize = 8*nx*ny*nz/1024*1024
            lchunks = False
            if dstchunksize > chunksize:
                lchunks = True
                nchunks = cpu_optimal(nx,ny,nz,mvar=1,maux=0,MBmin=chunksize)[1]
                print('nchunks {}'.format(nchunks))
                indx = np.array_split(np.arange(nx)+dstghost,nchunks[0])
                indy = np.array_split(np.arange(ny)+dstghost,nchunks[1])
                indz = np.array_split(np.arange(nz)+dstghost,nchunks[2])
                mx, my, mz = dsth5['settings']['mx'][0],\
                             dsth5['settings']['my'][0],\
                             dsth5['settings']['mz'][0]
                if not quiet:
                    print('nx {}, ny {}, nz {}'.format(nx, ny, nz))
                    print('mx {}, my {}, mz {}'.format(mx, my, mz))
            group = group_h5(dsth5, 'data', status='w')
            for key in srch5['data'].keys():
                print('remeshing '+key)
                if not lchunks:
                    var = local_remesh(srch5['data'][key][()],
                                       srch5['grid']['x'],srch5['grid']['y'],
                                       srch5['grid']['z'],dsth5['grid']['x'],
                                       dsth5['grid']['y'], dsth5['grid']['z'],
                                       quiet=quiet)
                    print('writing '+key+' shape {}'.format(var.shape))
                    dset = dataset_h5(group, key, status='w', data=var,
                                      overwrite=True, dtype=dtype)
                else:
                    dset = dataset_h5(group, key, status='w', shape=[mz,my,mx],
                                      overwrite=True, dtype=dtype)
                    print('writing '+key+' shape {}'.format([mz,my,mx]))
                    for iz in range(nchunks[2]):
                        n1, n2 = indz[iz][ 0]-dstghost,\
                                 indz[iz][-1]+dstghost
                        srcn1 = np.max(np.where(srch5['grid/z'][()]<
                                                dsth5['grid/z'][n1]))
                        srcn2 = np.min(np.where(srch5['grid/z'][()]>
                                                dsth5['grid/z'][n2]))
                        n1out = n1+dstghost
                        n2out = n2-dstghost+1
                        varn1 =  dstghost
                        varn2 = -dstghost
                        if iz == 0:
                            n1out = 0
                            varn1 = 0
                        if iz == nchunks[2]-1:
                            n2out = n2+1
                            varn2 = n2+1
                        if not quiet:
                            print('n1 {}, n2 {}, srcn1 {}, srcn2 {}'.format(
                                   n1,    n2,    srcn1,    srcn2))
                        for iy in range(nchunks[1]):
                            m1, m2 = indy[iy][ 0]-dstghost,\
                                     indy[iy][-1]+dstghost
                            srcm1 = np.max(np.where(srch5['grid/y'][()]<
                                                    dsth5['grid/y'][m1]))
                            srcm2 = np.min(np.where(srch5['grid/y'][()]>
                                                    dsth5['grid/y'][m2]))
                            m1out = m1+dstghost
                            m2out = m2-dstghost+1
                            varm1 =  dstghost
                            varm2 = -dstghost
                            if iy == 0:
                                m1out = 0
                                varm1 = 0
                            if iy == nchunks[1]-1:
                                m2out = m2+1
                                varm2 = m2+1
                            if not quiet:
                                print('m1 {}, m2 {}, srcm1 {}, srcm2 {}'.format(
                                       m1,    m2,    srcm1,    srcm2))
                            for ix in range(nchunks[0]):
                                l1, l2 = indx[ix][ 0]-dstghost,\
                                         indx[ix][-1]+dstghost
                                srcl1 = np.max(np.where(srch5['grid/x'][()]<
                                                        dsth5['grid/x'][l1]))
                                srcl2 = np.min(np.where(srch5['grid/x'][()]>
                                                        dsth5['grid/x'][l2]))
                                l1out = l1+dstghost
                                l2out = l2-dstghost+1
                                varl1 =  dstghost
                                varl2 = -dstghost
                                if ix == 0:
                                    l1out = 0
                                    varl1 = 0
                                if ix == nchunks[0]-1:
                                    l2out = l2+1
                                    varl2 = l2+1
                                if not quiet:
                                    print(
                                      'l1 {}, l2 {}, srcl1 {}, srcl2 {}'.format(
                                       l1,    l2,    srcl1,    srcl2))
                                if not quiet:
                                    print('remeshing '+key+' chunk {}'.format(
                                           [iz,iy,ix]))

                                var = local_remesh(
                                             srch5['data'][key][srcn1:srcn2+1,
                                                                srcm1:srcm2+1,
                                                                srcl1:srcl2+1],
                                             srch5['grid']['x'][srcl1:srcl2+1],
                                             srch5['grid']['y'][srcm1:srcm2+1],
                                             srch5['grid']['z'][srcn1:srcn2+1],
                                             dsth5['grid']['x'][l1:l2+1],
                                             dsth5['grid']['y'][m1:m2+1],
                                             dsth5['grid']['z'][n1:n2+1],
                                             quiet=quiet )
                                if not quiet:
                                    print('writing '+key+
                                                   ' shape {} chunk {}'.format(
                                                         var.shape, [iz,iy,ix]))
                                dset[n1out:n2out,
                                     m1out:m2out,
                                     l1out:l2out] = dtype(var[
                                                            varn1:varn2,
                                                            varm1:varm2,
                                                            varl1:varl2])
    dstsim.update()
    dstsim.change_value_in_file('src/cparam.local','ncpus', str(nprocs))
    dstsim.change_value_in_file('src/cparam.local','nprocx',str(ncpus[0]))
    dstsim.change_value_in_file('src/cparam.local','nprocy',str(ncpus[1]))
    dstsim.change_value_in_file('src/cparam.local','nprocz',str(ncpus[2]))
    dstsim.change_value_in_file('src/cparam.local','nxgrid',
                                                         str(dstsim.dim.nxgrid))
    #dstsim.change_value_in_file('src/cparam.local','nygrid',
    #                                                    str(dstsim.dim.nygrid))
    dstsim.change_value_in_file('src/cparam.local','nzgrid',
                                                         str(dstsim.dim.nzgrid))

    #cmd = 'source '+join(srcsim.path,'src','.moduleinfo')
    #os.system(cmd)
    #os.chdir(dstsim.path)
    #cmd = 'pc_setupsrc; make cleann'
    #os.system(cmd)
    #cmd = 'pc_build'
    #if hostfile: cmd = cmd + ' -f '+hostfile
    #process = sub.Popen(cmd.split(),stdout=sub.PIPE)
    #process = sub.Popen(cmd.split(),stdout=sub.PIPE)
    #output, error = process.communicate()
    #print(cmd,output,error)
    if srcprocs > nprocs:
        print('\n**********************************************************\n'+
              'remesh WARNING: {} procs reduced from {}.\n'.format(
              nprocs, srcprocs)+
              'Review multxyz {} and fracxyz {} for more\n'.format(
              multxyz,fracxyz)+
              'efficient parallel processing options.'+
              '\n**********************************************************\n')
    end_time = time.time()
    print('end at {} after {} seconds'.format(
                                     time.ctime(end_time),end_time-start_time))
# remains to copy other files and edit param files
