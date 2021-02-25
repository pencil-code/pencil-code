# calc_tensors_zaver.py
#
# NB: the tensor array is returned C-ordered:
# f[1-rank,2-rank,3-rank,nz,ny,nx,nt] for writing efficiently to hdf5.
# The corresponding fortran array in the pencil mean field module is
# f[nt,nx,ny,nz,nvar]
#
# Authors: Fred Gent (fred.gent.ncl@gmail.com)
#          Simo Tuomisto (simo.tuomisto@aalto.fi)
#
#--------------------------------------------------------------------------
#pc.calc.zav2h5(folder='data', filename='emftensors.h5', timereducer='mean', hdf5dir='.', rmfzeros=15, rmbzeros=5, t_correction=0, l_correction=False)
#
#where:
#
# for 10 time steps between reseting, and do cut 30% first and 40% at the end (like idl routine)
# choose: rmfzeros=3, rmbzeros=4
#
# for 50 time steps between reseting, and do cut 30% first and 40% at the end (like idl routine)
# choose: rmfzeros=15, rmbzeros=20
#
#  zav2h5:   function included in aver2h5.py
#  folder:   location of the data you want to store
#  filename= 'location_of_the_new_file/filename.h5'
#  timereducer is the kind of set you want to calculate: it can be ?mean?, ?mean_last?, or ?none? (takes the full time-series)
#  rmfzeros and rmbzeros remove the chosen number of time points from the surrounding of the coefficients resetting.
#  l_correction is a correction for Fred's Millennium run.
#

def zav2h5(
           folder='.',
           dataset='',
           filename='emftensors.h5',
           timereducer='mean',
           hdf5dir='data/',
           l_correction=True,
           t_correction=8972.,
           rmfzeros=4,
           rmbzeros=2,
           dgroup='emftensor',
          ):
    """
    If large dataset MPI may be required.
    Loads Averages object and applies tensors calculation and reforms for
    efficient writing to hdf5 for mean field module simulations.
    MPI call needs to be improved to avoid MemoryError for large files
    with read.aver(plane_list=['z'])
    timereducers needs to be expanded to include various smoothing options
    """
    import numpy as np
    from .. import read
    from ..read import aver
    from ..export import create_h5, fvars, create_aver_sph
#    from ..export import create_h5.fvars as fvars
#    from ..export import create_aver_sph
    from ..calc import tensors_sph
    import h5py
    import copy

    timereducers= {
                  'mean': lambda x,args:  np.mean(x,axis=-3,keepdims=True),
                                           #np.std(x,axis=-3)),
                  'mean_last': lambda x,args:  np.mean(np.take(
            x,np.arange(-int(args[0]),0,1),axis=-3),axis=-3,keepdims=True),
                  'none': lambda x, args: x
                  }
    if not timereducer in timereducers:
        raise ValueError(
              'timereducer "{}" undefined in timereducers'.format(
               timereducer)+' options: {}'.format(timereducers.keys()))
    if len(dataset)==0:
        dataset=timereducer
    with open('zaver.in', 'r') as f:
        zavers = f.read().splitlines()

    """ Find out if the calculation is parallel and distribute the arrays
        according to y-index and ipz=0 processor layout
    """
    try:
        from mpi4py import MPI

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()    # rank of processor on which this script runs
        size = comm.Get_size()    # number of  ~  ~  ~  ~

        l_mpi = True
        l_mpi = l_mpi and (size != 1)
    except ImportError:
        l_mpi = False
        rank = 0
        size = 1
        comm=None
    dim=read.dim()
    nx, nny = dim.nx, dim.ny
    ayindex=np.arange(nny)
    if l_mpi:
        y_chunks = np.array_split(ayindex,size, axis=0)
        yindex   = y_chunks[rank]
        ny       = yindex.size
    else:
        yindex=ayindex            # vector 0 ... nygrid-1
        ny=nny

    ncpus=dim.nprocx*dim.nprocy
    aprocs = np.arange(ncpus)     # vector 0 ... nprocx*nprocy-1
    if np.mod(ncpus,size) > 0:
        print('number of processes must divide {} cpus'.format(ncpus))
        quit()
    if l_mpi:
        if size>aprocs.size:
            nsz=size/aprocs.size
            for ii in range(1,nsz):
                tmproc=np.append(aprocs,aprocs)
            aprocs=np.sort(tmproc)
        proc_chunks = np.array_split(aprocs,size,axis=0)
        proc=proc_chunks[rank]
    else:
        proc=aprocs

    """Set up hdf5 file and create datasets in which to save the tensors
    """
    lskip_zeros   = rmfzeros+rmbzeros > 0
    if rank==0:                                # if root processor
        grid=read.grid(trim=True,quiet=True)   # read grid
        zav=read.aver(proc=0,plane_list=['z']) # read zaverages of root proc of PC run
        tensor=tensors_sph(                    # decompose into individual effect tensors
                       zav,
                       proc=proc[0],
                       rank=0,
                       lskip_zeros=lskip_zeros,
                       iy=[int(ny/2/dim.nprocy),],
                       timereducer=timereducers[timereducer],
                       #trargs=trargs,
                       rmfzeros=rmfzeros,
                       rmbzeros=rmbzeros,
                       l_correction=l_correction,
                       t_correction=t_correction,
                       dim=dim,
                       #tindex=tindex
                       )
        if 'mean' in dataset:
            nt=1
        else:
            nt=tensor.t.size
        create_aver_sph(
        hdf5dir+filename,
        dataset,
        fvars,
        (1, ny, nx, nt),
        (0,grid.y,grid.x,tensor.t),
        hdf5dir=hdf5dir,
        dgroup=dgroup
        )
    if l_mpi:
        imask=comm.bcast(tensor.imask, root=0)
    else:
        imask=tensor.imask

    import os

    if os.path.exists('data/averages/z.h5'):
        zav=aver(plane_list=['z'])     # read all averages
        tensor_buf=tensors_sph(        # calculate tensors
                           aver=zav,
                           rank=rank,
                           lskip_zeros=lskip_zeros,
                           timereducer=timereducers[timereducer],
                           #trargs=trargs,
                           rmfzeros=rmfzeros,
                           rmbzeros=rmbzeros,
                           l_correction=l_correction,
                           t_correction=t_correction,
                           dim=dim,
                           #tindex=tindex,
                           imask=imask
                           )
    else:
        yndx_tmp = np.array_split(yindex, dim.nprocy)
        # list of vectors ipy*ny/nprocy ... (ipy+1)*ny/nprocy - 1, ipy=0,nprocy-1

        for ipy in range(dim.nprocy):            # over all y processors of the PC run
            for ipx in range(dim.nprocx):        # over all x processors of the PC run

                iproc=dim.nprocx*ipy+ipx         # proc rank of the PC run (0 ... nprocx*nprocy-1)
                yndx=yndx_tmp[ipy]-ipy*int(dim.nygrid/dim.nprocy)

                zav=aver(proc=iproc,plane_list=['z'])     # read averages from proc iproc

                print('calculating tensors on proc {0} rank {1}'.format(iproc,rank))
                """
                if iproc==1:             # as there is corrupted data on proc 1
                    with open('zaver.in', 'r') as f:
                        zavers = f.read().splitlines()
                    for zaver in  zavers:
                        zav.z.__setattr__(zaver,np.insert(
                                zav.z.__getattribute__(zaver),3766,
                                0.5*(zav.z.__getattribute__(zaver)[3766]+
                                zav.z.__getattribute__(zaver)[3767]),axis=0))
                        zav.t=np.insert(zav.t,3766,0.5*(zav.t[3766]+zav.t[3767]),axis=0)
                """
                tensor_buf=tensors_sph(   # calculate tensors
                                   aver=zav,
                                   proc=iproc,
                                   rank=rank,
                                   lskip_zeros=lskip_zeros,
                                   iy=yndx,
                                   timereducer=timereducers[timereducer],
                                   #trargs=trargs,
                                   rmfzeros=rmfzeros,
                                   rmbzeros=rmbzeros,
                                   l_correction=l_correction,
                                   t_correction=t_correction,
                                   dim=dim,
                                   #tindex=tindex,
                                   imask=imask
                                   )
                if ipx==0:
                    tensor=copy.deepcopy(tensor_buf)
                else:
                    for field, comp in fvars:
                        setattr(tensor,field,np.concatenate((tensor.__getattribute__(field),
                                                             tensor_buf.__getattribute__(field)),
                                                             axis=len(comp)+2))

        if l_mpi:
            comm.barrier()
            ds=h5py.File(hdf5dir+filename, 'a', driver='mpio', comm=comm)
        else:
            ds=h5py.File(hdf5dir+filename, 'a')     # open HDF5 file

        for field, comp in fvars:

            print('writing {0} from rank {1} for proc {2}'.format(field, rank, iproc))

            dsname='{0}/{1}/{2}'.format(dgroup,field,dataset)
            if len(comp)==1:
                ds[dsname][:,:,yndx_tmp[ipy],:]=tensor.__getattribute__(field)
            elif len(comp)==2:
                ds[dsname][:,:,:,yndx_tmp[ipy],:]=tensor.__getattribute__(field)
            else:
                ds[dsname][:,:,:,:,yndx_tmp[ipy],:]=tensor.__getattribute__(field)
        ds.close()
