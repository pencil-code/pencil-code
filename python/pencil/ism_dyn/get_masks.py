# get_masks.py
#
# 05-may-20
# Author: F. Gent (fred.gent.ncl@gmail.com).
#
""" Derive auxilliary data and other diagnostics from var.h5 file and 
    save to new h5 file
 
    uses:
      compute 'data' arrays of size [nz,ny,nx] as required
      store 'time' of snapshot
      compute 'masks' for example by temperature phase
      compute summary statistics 'stats' 
      compute 'structure' functions as required
"""
import numpy as np
from ..math import cpu_optimal
from ..io import open_h5, group_h5, dataset_h5
from .. import read 
import os

def thermal_decomposition(
    ss, pars, unit_key='unit_entropy', ent_cut=[2.32e9,]):
    """
    call signature:

    thermal_decomposition(ss, pars, unit='unit_entropy', ent_cut=[2.32e9,])
    
    Keyword arguments:
        ss:       dataset used for masks, default 'ss', alternate e.g.'tt'
        pars:     Param() object required for units rescaling
        unit_key: label of physical units in pars to apply to code values
        ent_cut:  list of boundary mask values, default see thesis
                  http://hdl.handle.net/10443/1755 Figure 5.10
                  may have multiple boundaries
    """
    temp_masks = list()
    hh = np.ma.array(np.copy(ss))
    for ent in ent_cut:
        hcut = ent/pars.__getattribute__(unit_key)
        hh[np.where(ss<hcut)] = np.ma.masked
        temp_masks.append(hh.mask)
    print('thermal_decomposition', temp_masks[0].shape,len(temp_masks))
    return temp_masks

def derive_masks(sim_path, src, dst, data_key='data/ss', par=[], comm=None,
                overwrite=False, rank=0, size=1, nghost=3, status='a',
                chunksize = 1000.0, quiet=True, nmin=32,
                ent_cuts=[2.32e9,], mask_keys=['hot',], unit_key='unit_entropy'
               ):
    if comm:
        overwrite = False
    if isinstance(par, list):
        os.chdir(sim_path)
        par = read.param(quiet=True,conflicts_quiet=True)
    #get data dimensions
    nx, ny, nz = src['settings']['nx'][0],\
                 src['settings']['ny'][0],\
                 src['settings']['nz'][0]
    mx, my, mz = src['settings']['mx'][0],\
                 src['settings']['my'][0],\
                 src['settings']['mz'][0]
    #split data into manageable memory chunks
    dstchunksize = 8*nx*ny*nz/1024*1024
    if dstchunksize > chunksize:
        nchunks = cpu_optimal(nx,ny,nz,quiet=quiet,
                              mvar=src['settings/mvar'][0],
                              maux=src['settings/maux'][0],
                              MBmin=chunksize,nmin=nmin,size=size)[1]
    else:
        nchunks = [1,1,1]
    print('nchunks {}'.format(nchunks)) 
    # for mpi split chunks across processes
    # for mpi split chunks across processes
    if size > 1:
        locindx = np.array_split(np.arange(nx)+nghost,nchunks[0]) 
        locindy = np.array_split(np.arange(ny)+nghost,nchunks[1]) 
        locindz = np.array_split(np.arange(nz)+nghost,nchunks[2])
        indx = [locindx[np.mod(rank+int(rank/nchunks[2])
                                   +int(rank/nchunks[1]),nchunks[0])]]
        indy = [locindy[np.mod(rank+int(rank/nchunks[2]),nchunks[1])]]
        indz = [locindz[np.mod(rank,nchunks[2])]]
        allchunks = 1
    else:
        locindx = np.array_split(np.arange(nx)+nghost,nchunks[0]) 
        locindy = np.array_split(np.arange(ny)+nghost,nchunks[1]) 
        locindz = np.array_split(np.arange(nz)+nghost,nchunks[2])
        indx = np.array_split(np.arange(nx)+nghost,nchunks[0]) 
        indy = np.array_split(np.arange(ny)+nghost,nchunks[1]) 
        indz = np.array_split(np.arange(nz)+nghost,nchunks[2])
        allchunks = nchunks[0]*nchunks[1]*nchunks[2]
    # ensure derived variables are in a list
    if isinstance(mask_keys, list):
        mask_keys = mask_keys
    else:
        mask_keys = [mask_keys]
    # initialise group 
    group = group_h5(dst, 'masks', status='a', overwrite=overwrite,
                     comm=comm, rank=rank, size=size)
    for key in mask_keys:
        ne = len(ent_cuts)
        dataset_h5(group, key, status=status, shape=[ne,mz,my,mx],
                      comm=comm, size=size, rank=rank,
                      overwrite=overwrite, dtype=np.bool_)
        print('writing '+key+' shape {}'.format([ne,mz,my,mx]))
        for ichunk in range(allchunks):
            for iz in [indz[np.mod(ichunk,nchunks[2])]]:
                n1, n2 = iz[ 0]-nghost,\
                         iz[-1]+nghost+1
                n1out = n1+nghost
                n2out = n2-nghost
                varn1 =  nghost
                varn2 = -nghost
                if iz[0] == locindz[0][0]:
                    n1out = 0
                    varn1 = 0
                if iz[-1] == locindz[-1][-1]:
                    n2out = n2
                    varn2 = n2
                for iy in [indy[np.mod(ichunk+
                                   int(ichunk/nchunks[2]),nchunks[1])]]:
                    m1, m2 = iy[ 0]-nghost,\
                             iy[-1]+nghost+1
                    m1out = m1+nghost
                    m2out = m2-nghost
                    varm1 =  nghost
                    varm2 = -nghost
                    if iy[0] == locindy[0][0]:
                        m1out = 0
                        varm1 = 0
                    if iy[-1] == locindy[-1][-1]:
                        m2out = m2
                        varm2 = m2
                    for ix in [indx[np.mod(ichunk+int(ichunk/nchunks[2])
                                   +int(ichunk/nchunks[1]),nchunks[0])]]:
                        l1, l2 = ix[ 0]-nghost,\
                                 ix[-1]+nghost+1
                        l1out = l1+nghost
                        l2out = l2-nghost
                        varl1 =  nghost
                        varl2 = -nghost
                        if ix[0] == locindx[0][0]:
                            l1out = 0
                            varl1 = 0
                        if ix[-1] == locindx[-1][-1]:
                            l2out = l2
                            varl2 = l2
                        if data_key in src.keys():
                            ss = src[data_key][n1:n2,m1:m2,l1:l2]
                        else:
                            if data_key in dst.keys():
                                ss = dst[data_key][n1:n2,m1:m2,l1:l2]
                            else:
                                print('masks: '+data_key+' does not exist in ',
                                        src,'or',dst)
                                return 1
                        masks = thermal_decomposition(ss,par,unit_key=unit_key,
                                                     ent_cut=ent_cuts
                                                    )
                        cut = 0
                        for mask in masks:
                            dst['masks'][key][cut,n1out:n2out,
                                         m1out:m2out,
                                         l1out:l2out] = mask[
                                                     varn1:varn2,
                                                     varm1:varm2,
                                                     varl1:varl2]
                            cut += 1
