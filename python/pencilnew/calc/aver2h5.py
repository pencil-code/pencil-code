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

def zav2h5(
           folder='.',
           dataset='',
           filename='data/emftensors.h5',
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
    from ..export import fvars, create_aver_sph 
    from ..calc import tensors_sph 
    import h5py
    timereducers= { 
                  'mean': lambda x,args:  np.mean(x,axis=-3,keepdims=True),
                  'mean_last': lambda x,args:  np.mean(np.take(
            x,np.arange(-int(args[0]),0,1),axis=-3),axis=-3,keepdims=True),
                  'none': lambda x, args: x,
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
        rank = comm.Get_rank()
        size = comm.Get_size()

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
        yindex=ayindex
        ny=nny
    ncpus=dim.nprocx*dim.nprocy
    aprocs = np.arange(ncpus)
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
    print('proc.size',proc.size)
    """Set up hdf5 file and create datasets in which to save the tensors
    """
    lskip_zeros   = rmfzeros+rmbzeros > 0
    if rank==0:
        grid=read.grid(trim=True,quiet=True)
        zav=read.aver(proc=0,plane_list=['z'])
        tensor=tensors_sph(
                       zav,
                       proc=proc[0],
                       rank=rank,
                       lskip_zeros=lskip_zeros,
                       iy=[ny/2/dim.nprocy,],
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
        filename,
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
    yndx_tmp = np.array_split(yindex, proc.size, axis=0)
    for iproc in proc:
        yndx=yndx_tmp[iproc]-iproc*(dim.nygrid/dim.nprocy)
        print('yndx {} from yindex {}'.format(yndx,yndx_tmp[iproc]))
        zav=read.aver(proc=iproc,plane_list=['z'])
        print('calculating tensors on proc {0} rank {1}'.format(iproc,rank))
        if iproc==1:
            with open('zaver.in', 'r') as f: 
                zavers = f.read().splitlines() 
            for zaver in  zavers:
                zav.z.__setattr__(zaver,np.insert(
                         zav.z.__getattribute__(zaver),3766,
                         0.5*(zav.z.__getattribute__(zaver)[3766]+
                              zav.z.__getattribute__(zaver)[3767]),axis=0))
                zav.t=np.insert(zav.t,3766,0.5*(zav.t[3766]+zav.t[3767]),axis=0)
        tensor=tensors_sph(
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
        if l_mpi:
            comm.barrier()
            ds=h5py.File(filename, 'a', driver='mpio', comm=comm)
        else:
            ds=h5py.File(filename, 'a')
        for field, comp in fvars:
            print('writing {0} from rank {1} for proc {2}'.format(
                   field, rank, proc[iproc]))
            if len(comp)==1:
                ds['{0}/{1}/{2}'.format(dgroup,field,dataset)][:,:,yndx_tmp[iproc],:]=\
                    tensor.__getattribute__(field)
            elif len(comp)==2:
                ds['{0}/{1}/{2}'.format(dgroup,field,dataset)][:,:,:,yndx_tmp[iproc],:]=\
                    tensor.__getattribute__(field)
            else:
                ds['{0}/{1}/{2}'.format(dgroup,field,dataset)][:,:,:,:,yndx_tmp[iproc],:]=\
                    tensor.__getattribute__(field)
        ds.close()
