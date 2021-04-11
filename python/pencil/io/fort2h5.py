# fort2h5.py
#
# Read existing Fortran unformatted simulation data and write as hdf5.
#
#
# Author: F. Gent (fred.gent.ncl@gmail.com).
#
"""
Contains the functions to read old Fortran binary simulation data and
write snapshots in hdf5 (data/allprocs/VAR*.h5 files), video slices to
data/slices/uu1_xy.h5, etc. and averages to data/averages/xy.h5, etc
"""

import time

def var2h5(newdir, olddir, allfile_names, todatadir, fromdatadir, snap_by_proc,
           precision, lpersist, quiet, nghost, settings, param, grid,
           x, y, z, lshear, lremove_old_snapshots, indx,
           trimall=False, l_mpi=False, driver=None, comm=None, rank=0, size=1
          ):

    """
    Copy a simulation snapshot set written in Fortran binary to hdf5.

    call signature:

    var2h5(newdir, olddir, allfile_names, todatadir, fromdatadir, snap_by_proc,
           precision, lpersist, quiet, nghost, settings, param, grid,
           x, y, z, lshear, lremove_old_snapshots, indx,
           trimall=False, l_mpi=False, driver=None, comm=None, rank=0, size=1
          )

    Keyword arguments:

    *newdir*:
      String path to simulation destination directory.

    *olddir*:
      String path to simulation destination directory.

    *allfile_names*:
      A list of names of the snapshot files to be written, e.g. VAR0.

    *todatadir*:
      Directory to which the data is stored.

    *fromdatadir*:
      Directory from which the data is collected.

    *snap_by_proc*:
      Read and write snapshots by procdir of the fortran binary tree

    *precision*:
      Single 'f' or double 'd' precision for new data.

    *lpersist*:
      option to include persistent variables from snapshots.

    *quiet*
      Option not to print output.

    *nghost*:
      Number of ghost zones.

    *settings*
      simulation properties.

    *param*
      simulation Param object.

    *grid*
      simulation Grid object.

    *xyz*:
      xyz arrays of the domain with ghost zones.

    *lshear*:
      Flag for the shear.

    *lremove_old_snapshots*:
      If True the old snapshots will be deleted once the new snapshot has
      been saved.

    *indx*:
      List of variable indices in the f-array.

    *trimall*:
      Strip ghost zones from snapshots

    *l_mpi*:
      Applying MPI parallel process

    *driver*:
      HDF5 file io driver either None or mpio

    *comm*:
      MPI library calls

    *rank*:
      Integer ID of processor

    *size*:
      Number of MPI processes
    """
    import os
    from os.path import exists, join
    import numpy as np
    import glob
    from .. import read
    from .. import sim
    from . import write_h5_snapshot
    import sys
    import time
    import subprocess as sub

    if isinstance(allfile_names, list):
        allfile_names = allfile_names
    else:
        allfile_names = [allfile_names]
    #proceed to copy each snapshot in varfile_names
    nprocs = settings['nprocx']*settings['nprocy']*settings['nprocz']
    if l_mpi:
        if not snap_by_proc:
            file_names = np.array_split(allfile_names, size)
            if 'VARd1' in allfile_names:
                varfile_names = file_names[size-rank-1]
            else:
                varfile_names = file_names[rank]
        else:
            os.chdir(olddir)
            if size > nprocs:
                nnames = len(allfile_names)
                if size > nnames*nprocs:
                    file_names = np.array_split(allfile_names, nnames)
                    varfile_names = file_names[np.mod(rank,nnames)]
                    nprocsplit = int(size/nnames)
                    iprocs = np.array_split(np.arange(nprocs),nprocs)
                    procs = iprocs[np.mod(rank,nprocs)]
                else:
                    file_names = np.array_split(allfile_names, nnames)
                    varfile_names = file_names[np.mod(rank,nnames)]
                    if nnames > size:
                        procs = np.arange(nprocs)
                    else:
                        nproc_per_fname = int(size/nnames)
                        isize = np.int(np.mod(rank,nnames)/nproc_per_fname)
                        if np.mod(isize,nproc_per_fname+1) == 0:
                            npf = nproc_per_fname+1
                            iprocs = np.array(np.array_split(np.arange(nprocs),
                                              npf)).T
                        else:
                            npf = nproc_per_fname
                            iprocs = np.array(np.array_split(np.arange(nprocs),
                                              npf)).T
                    procs = iprocs[np.mod(int((rank*nnames)/size),npf)]
            #if 'VARd1' in allfile_names:
            #    file_names = np.array_split(allfile_names, size)
            #    varfile_names = file_names[size-rank-1]
            else:
                if np.mod(nprocs,size) > 0:
                    procs = np.arange(nprocs+size-np.mod(nprocs,size))
                    procs[-size+np.mod(nprocs,size):] = np.arange(size-np.mod(nprocs,size))
                else:
                    procs = np.arange(nprocs)
                iprocs = np.array_split(procs,size)
                procs = iprocs[rank]
                varfile_names = allfile_names
            print('rank {} procs:'.format(rank),procs)
            sys.stdout.flush()
    else:
        varfile_names = allfile_names
        procs = np.arange(nprocs)
    if len(varfile_names) > 0:
        for file_name in varfile_names:
            #load Fortran binary snapshot
            if not quiet:
                print('rank {}:'.format(rank)+'saving '+file_name)
                sys.stdout.flush()
            if snap_by_proc:
                if len(procs) > 0:
                    proctime = time.time()
                    for proc in procs:
                        os.chdir(olddir)
                        if np.mod(proc,size) == size-1:
                            print('rank {}:'.format(rank)+'saving '+file_name+
                                          ' on proc{}\n'.format(proc),
                                   time.ctime())
                            sys.stdout.flush()
                        procdim = read.dim(proc=proc)
                        var = read.var(file_name, datadir=fromdatadir,
                                       quiet=quiet, lpersist=lpersist,
                                       trimall=trimall, proc=proc)
                        try:
                            var.deltay
                            lshear = True
                        except:
                            lshear = False

                        if lpersist:
                            persist = {}
                            for key in read.record_types.keys():
                                try:
                                    persist[key] = var.__getattribute__(key)[()]
                                    if (type(persist[key][0])==str):
                                        persist[key][0] = \
                                           var.__getattribute__(key)[0].encode()
                                except:
                                    pass
                        else:
                            persist = None
                        if np.mod(proc,size) == size-1:
                            print('rank {}:'.format(rank)+'loaded '+file_name+
                                  ' on proc{} in {} seconds'.format(
                                  proc,time.time()-proctime))
                            sys.stdout.flush()
                        #write data to h5
                        os.chdir(newdir)
                        write_h5_snapshot(var.f, file_name=file_name, state='a',
                                    datadir=todatadir, precision=precision,
                                    nghost=nghost, persist=persist, proc=proc,
                                    procdim=procdim, settings=settings,
                                    param=param, grid=grid, lghosts=True,
                                    indx=indx, t=var.t, x=x, y=y, z=z,
                                    quiet=quiet, rank=rank, size=size,
                                    lshear=lshear, driver=driver, comm=comm)
                        if np.mod(proc,size) == size-1:
                            print('rank {}:'.format(rank)+'written '+file_name+
                                  ' on proc{} in {} seconds'.format(
                                  proc,time.time()-proctime))
                            sys.stdout.flush()
                        proctime = time.time()
            else:
                var = read.var(file_name, datadir=fromdatadir, quiet=quiet,
                               lpersist=lpersist, trimall=trimall)
                try:
                    var.deltay
                    lshear = True
                except:
                    lshear = False

                if lpersist:
                    persist = {}
                    for key in read.record_types.keys():
                        try:
                            persist[key] = var.__getattribute__(key)[()]
                            if (type(persist[key][0])==str):
                                persist[key][0] = \
                                           var.__getattribute__(key)[0].encode()
                        except:
                            pass
                else:
                    persist = None
                #write data to h5
                os.chdir(newdir)
                write_h5_snapshot(var.f, file_name=file_name, datadir=todatadir,
                                  precision=precision, nghost=nghost,
                                  persist=persist, settings=settings,
                                  param=param, grid=grid, lghosts=True,
                                  indx=indx, t=var.t, x=x, y=y, z=z,
                                  lshear=lshear, driver=None, comm=None)
            if lremove_old_snapshots:
                os.chdir(olddir)
                cmd = "rm -f "+join(fromdatadir, 'proc*', file_name)
                process = sub.Popen(cmd.split(),stdout=sub.PIPE)
                output, error = process.communicate()
                print(cmd,output,error)
                #os.system(cmd)
            del(var)
            print('fred: rank {} end vars'.format(rank))

def slices2h5(newdir, olddir, grid,
              todatadir='data/slices', fromdatadir='data', precision='d',
              quiet=True, lremove_old_slices=False, lsplit_slices=False,
              l_mpi=False, driver=None, comm=None, rank=0, size=1,
              vlarge=1000000000):

    """
    Copy a simulation set of video slices written in Fortran binary to hdf5.

    call signature:

    slices2h5(newdir, olddir, grid,
              todatadir='data/slices', fromdatadir='data', precision='d',
              quiet=True, lremove_old_slices=False, lsplit_slices=False,
              l_mpi=False, driver=None, comm=None, rank=0, size=1,
              vlarge=1000000000):

    Keyword arguments:

    *newdir*:
      String path to simulation destination directory.

    *olddir*:
      String path to simulation destination directory.

    *grid*
      simulation Grid object.

    *todatadir*:
      Directory to which the data is stored.

    *fromdatadir*:
      Directory from which the data is collected.

    *precision*:
      Single 'f' or double 'd' precision for new data.

    *quiet*:
      Option not to print output.

    *lremove_old_slices*:
      If True the old video slices will be deleted once the new slices have
      been saved.

    *lsplit_slices*:
      Read the video slices in iterative chunks for large memory

    *l_mpi*:
      Applying MPI parallel process

    *driver*:
      HDF5 file io driver either None or mpio

    *comm*:
      MPI library calls

    *rank*:
      Integer ID of processor

    *size*:
      Number of MPI processes

    *vlarge*:
      Limit size of file without chunking
    """

    import os
    from os.path import exists, join
    import numpy as np
    from .. import read
    from .. import sim
    from . import write_h5_slices
    import sys
    import subprocess as sub

    #copy old video slices to new h5 sim
    os.chdir(olddir)
    #identify the coordinates and positions of the slices
    coordinates = {}
    positions = {}
    readlines1 = open('data/slice_position.dat','r').readlines()
    readlines2 = open('data/proc0/slice_position.dat','r').readlines()
    lines1, lines2 = [],[]
    for line in readlines1:
        lines1.append(int(line.split(' ')[-1].split('\n')[0]))
    """In newer binary sims lines2 obtains 7 strings below, but older sims may
    only yield the integer coordinates, so lines2 is hardcoded. The current
    version of the Pencil Code has 7 potential slices, but earlier versions
    may not. If your sim does not conform to this arrangement edit/copy this
    module and set lines1 and lines2 manually from data/slice_position.dat and
    the extensions present in your slice_*.xy  etc.
    """
    #check simulation includes the slice keys in data/proc*/slice_position.dat
    try:
        int(int(readlines2[0].split(' ')[-1].split('\n')[0]))
        lines2=['xy', 'xy2', 'xy3', 'xy4', 'xz', 'xz2', 'yz']
    except:
        for line in readlines2:
            lines2.append(line.split(' ')[-1].split('\n')[0].lower())
    #check if number of slice options as expected
    try:
        len(lines1)==7
    except ValueError:
        if rank == 0:
            print("ERROR: slice keys and positions must be set, "+
                  "see lines 212...")
            sys.stdout.flush()
        return -1
    for key, num in zip(lines2, lines1):
        if comm:
            key, num = comm.bcast([key,num],root=0)
        if num > 0:
            if 'xy' in key:
                positions[key] = grid.z[num-1]
            if 'xz' in key:
                positions[key] = grid.y[num-1]
            if 'yz' in key:
                positions[key] = grid.x[num-1]
            coordinates[key] = num
    if l_mpi:
        import glob
        slice_lists = glob.glob(join(fromdatadir,'slice_*'))
        slice_lists.remove(join(fromdatadir,'slice_position.dat'))
        slice_lists = np.array_split(slice_lists,size)
        slice_list = slice_lists[rank]
        print('rank {} slice_list {}'.format(rank,slice_list))
        if not quiet:
            print('rank', rank, 'slice_list', slice_list)
            sys.stdout.flush()
        if len(slice_list) > 0:
            for field_ext in slice_list:
                print('rank', rank, 'field_ext', field_ext)
                field=str.split(str.split(field_ext,'_')[-1],'.')[0]
                extension=str.split(str.split(field_ext,'_')[-1],'.')[1]
                os.chdir(olddir)
                vslice = read.slices(field=field,
                                     extension=extension, quiet=quiet)
                os.chdir(newdir)
                if not quiet:
                    print('rank', rank, 'writing', field_ext)
                    sys.stdout.flush()
                write_h5_slices(vslice, coordinates, positions,
                                datadir=todatadir, precision=precision,
                                quiet=quiet)
                if lremove_old_slices:
                    os.chdir(olddir)
                    cmd = "rm -f "+field_ext
                    process = sub.Popen(cmd.split(),stdout=sub.PIPE)
                    output, error = process.communicate()
                    print(cmd,output,error)
                    #os.system(cmd)
        else:
            vslice = None
    else:
        vslice = read.slices(quiet=quiet)
        #write new slices in hdf5
        os.chdir(newdir)
        write_h5_slices(vslice, coordinates, positions, datadir=todatadir,
                        precision=precision, quiet=quiet)
        if lremove_old_slices:
            os.chdir(olddir)
            cmd = "rm -f "+join(fromdatadir, 'proc*', 'slice_*')
            process = sub.Popen(cmd.split(),stdout=sub.PIPE)
            output, error = process.communicate()
            print(cmd,output,error)
            #os.system(cmd)
    print('fred: rank {} end slices'.format(rank))
    del(vslice)

def aver2h5(newdir, olddir,
            todatadir='data/averages', fromdatadir='data', l2D=True,
            precision='d', quiet=True, lremove_old_averages=False,
            aver_by_proc=False,
            laver2D=False, l_mpi=False, driver=None, comm=None, rank=0, size=1):

    """
    Copy a simulation set of video slices written in Fortran binary to hdf5.

    call signature:

    aver2h5(newdir, olddir,
            todatadir='data/averages', fromdatadir='data', l2D=True,
            precision='d', quiet=True, lremove_old_averages=False,
            aver_by_proc=False,
            laver2D=False, l_mpi=False, driver=None, comm=None, rank=0, size=1):

    Keyword arguments:

    *newdir*:
      String path to simulation destination directory.

    *olddir*:
      String path to simulation destination directory.

    *todatadir*:
      Directory to which the data is stored.

    *fromdatadir*:
      Directory from which the data is collected.

    *l2D*
     Option to include 2D averages if the file sizes are not too large

    *precision*:
      Single 'f' or double 'd' precision for new data.

    *quiet*
      Option not to print output.

    *lremove_old_averages*:
      If True the old averages data will be deleted once the new h5 data
      has been saved.

    *aver_by_proc*
      Option to read old binary files by processor and write in
      parallel

    *laver2D*
      If True apply to each plane_list 'y', 'z' and load each variable
      sequentially

    *l_mpi*:
      Applying MPI parallel process

    *driver*:
      HDF5 file io driver either None or mpio

    *comm*:
      MPI library calls

    *rank*:
      Integer ID of processor

    *size*:
      Number of MPI processes
    """

    import os
    from os.path import exists, join
    import numpy as np
    from .. import read
    from .. import sim
    from . import write_h5_averages
    import sys
    import subprocess as sub

    if laver2D:
        os.chdir(olddir)
        for xl in ['y','z']:
            if exists(xl+'aver.in'):
                if exists(join(fromdatadir,'t2davg.dat')):
                    f = open(join(fromdatadir,'t2davg.dat'))
                    niter = int(f.readline().split(' ')[-1].strip('\n'))-1
                else:
                    if not aver_by_proc:
                        av = read.aver(plane_list=xl, proc=0, var_index=0)
                        niter = av.t.size
                    else:
                        niter = None
                if aver_by_proc:
                    dim = read.dim()
                    if xl == 'y':
                        nproc = dim.nprocz
                    if xl == 'z':
                        nproc = dim.nprocy
                    all_list = np.array_split(np.arange(nproc), size)
                    proc_list = list(all_list[rank])
                    os.chdir(olddir)
                    if len(proc_list) > 0:
                        for proc in proc_list:
                            print('reading '+xl+'averages on proc', proc)
                            sys.stdout.flush()
                            av = read.aver(plane_list=xl, proc=proc)
                            procdim = read.dim(proc=proc)
                            write_h5_averages(av, file_name=xl, datadir=todatadir,
                                nt=niter, precision=precision, append=True,
                                aver_by_proc=True, nproc=nproc,
                                proc=proc, dim=dim, procdim=procdim, quiet=quiet,
                                driver=driver, comm=comm, rank=rank, size=size)
                    del(av)
                else:
                    all_list = np.array_split(np.arange(niter), size)
                    iter_list = list(all_list[rank])
                    os.chdir(olddir)
                    print('reading '+xl+'averages on rank', rank)
                    sys.stdout.flush()
                    av = read.aver(plane_list=xl, iter_list=iter_list)
                    os.chdir(newdir)
                    write_h5_averages(av, file_name=xl, datadir=todatadir,
                                      nt=niter, precision=precision,
                                      append=False, indx=iter_list, quiet=quiet,
                                      driver=driver, comm=comm, rank=rank,
                                      size=size)
                    del(av)
    else:
        #copy old 1D averages to new h5 sim
        os.chdir(olddir)
        plane_list = []
        for xl in ['xy','xz','yz']:
            if exists(xl+'aver.in'):
                plane_list.append(xl)
        if rank == size-1 or not l_mpi:
            if len(plane_list) > 0:
                av = read.aver(plane_list=plane_list)
                os.chdir(newdir)
                for key in av.__dict__.keys():
                    if not key in 't':
                        write_h5_averages(av, file_name=key, datadir=todatadir,
                                          precision=precision, quiet=quiet,
                                          driver=driver, comm=None, rank=None,
                                          size=size)
                del(av)
            if lremove_old_averages:
                os.chdir(olddir)
                cmd = "rm -f "+join(fromdatadir, '*averages.dat')
                process = sub.Popen(cmd.split(),stdout=sub.PIPE)
                output, error = process.communicate()
                print(cmd,output,error)
                #os.system(cmd)
            if l2D:
                plane_list = []
                os.chdir(olddir)
                for xl in ['x','y','z']:
                    if exists(xl+'aver.in'):
                        plane_list.append(xl)
                if len(plane_list) > 0:
                    for key in plane_list:
                        os.chdir(olddir)
                        av = read.aver(plane_list=key)
                        os.chdir(newdir)
                        write_h5_averages(av, file_name=key, datadir=todatadir,
                                          precision=precision, quiet=quiet,
                                          driver=None, comm=None)
                    del(av)
    if lremove_old_averages:
        if l_mpi:
            comm.Barrier()
        os.chdir(olddir)
        cmd = "rm -f "+join(fromdatadir, '*averages.dat')
        if rank == 0:
            process = sub.Popen(cmd.split(),stdout=sub.PIPE)
            output, error = process.communicate()
            print(cmd,output,error)
            #os.system(cmd)
    print('fred: rank {} end vars'.format(rank))

def sim2h5(newdir='.', olddir='.', varfile_names=None,
           todatadir='data/allprocs', fromdatadir='data',
           precision='d', nghost=3, lpersist=True,
           x=None, y=None, z=None, lshear=False,
           snap_by_proc=False, aver_by_proc=False,
           lremove_old_snapshots=False,
           lremove_old_slices=False, lread_all_videoslices=False,
           vlarge=100000000,
           lremove_old_averages=False, execute=False, quiet=True,
           l2D=True, lvars=True, lvids=True, laver=True, laver2D=False,
           lremove_deprecated_vids=False, lsplit_slices=False,
          ):
    """
    Copy a simulation object written in Fortran binary to hdf5.
    The default is to copy all snapshots from/to the current simulation
    directory. Optionally the old files can be removed to

    call signature:

    sim2h5(newdir='.', olddir='.', varfile_names=None,
           todatadir='data/allprocs', fromdatadir='data',
           precision='d', nghost=3, lpersist=False,
           x=None, y=None, z=None, lshear=False,
           snap_by_proc=False, aver_by_proc=False,
           lremove_old_snapshots=False,
           lremove_old_slices=False, lread_all_videoslices=True,
           lremove_old_averages=False, execute=False, quiet=True,
           l2D=True, lvars=True, lvids=True, laver=True)

    Keyword arguments:

    *olddir*:
      String path to simulation source directory.
      Path may be relative or absolute.

    *newdir*:
      String path to simulation destination directory.
      Path may be relative or absolute.

    *varfile_names*:
      A list of names of the snapshot files to be written, e.g. VAR0
      If None all varfiles in olddir+'/data/proc0/' will be converted

    *todatadir*:
      Directory to which the data is stored.

    *fromdatadir*:
      Directory from which the data is collected.

    *precision*:
      Single 'f' or double 'd' precision for new data.

    *nghost*:
      Number of ghost zones.
      TODO: handle switching size of ghost zones.

    *lpersist*:
      option to include persistent variables from snapshots.

    *xyz*:
      xyz arrays of the domain with ghost zones.
      This will normally be obtained from Grid object, but facility to
      redefine an alternative grid value.

    *lshear*:
      Flag for the shear.

    *execute*:
      optional confirmation required if lremove_old.

    *lremove_old_snapshots*:
      If True the old snapshot data will be deleted once the new h5 data
      has been saved.

    *lremove_old_slices*:
      If True the old video slice data will be deleted once the new h5 data
      has been saved.

    *lremove_old_averages*:
      If True the old averages data will be deleted once the new h5 data
      has been saved.

    *aver_by_proc*
      Option to read old binary files by processor and write in
      parallel

    *laver2D*
      If True apply to each plane_list 'y', 'z' and load each variable
      sequentially

    *l_mpi*:
      Applying MPI parallel process

    *driver*:
      HDF5 file io driver either None or mpio

    *comm*:
      MPI library calls

    *rank*:
      Integer ID of processor

    *size*:
      Number of MPI processes
    """

    import glob
    import numpy as np
    import os
    from os.path import exists, join
    import subprocess as sub
    import sys

    from .. import read
    from .. import sim
    from . import write_h5_grid
    from pencil.util import is_sim_dir

    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        driver='mpio'
        l_mpi = True
        l_mpi = l_mpi and (size != 1)
    except ImportError:
        comm = None
        driver=None
        rank = 0
        size = 1
        l_mpi = False
    if not l_mpi:
        comm = None
        driver=None
    print('rank {} and size {}'.format(rank,size))
    sys.stdout.flush()
    if rank == size-1:
        print('l_mpi',l_mpi)
        sys.stdout.flush()

    #test if simulation directories
    if newdir == '.':
        newdir = os.getcwd()
    if olddir == '.':
        olddir = os.getcwd()
    os.chdir(olddir)
    if not is_sim_dir():
        if rank == 0:
            print("ERROR: Directory ("+olddir+") needs to be a simulation")
            sys.stdout.flush()
        return -1
    if newdir != olddir:
        if not exists(newdir):
            cmd = 'pc_newrun -s '+newdir
            if rank == size-1:
                process = sub.Popen(cmd.split(),stdout=sub.PIPE)
                output, error = process.communicate()
                print(cmd,output,error)
                #os.system(cmd)
            if comm:
                comm.Barrier()
        os.chdir(newdir)
        if not is_sim_dir():
            if rank == 0:
                print("ERROR: Directory ("+newdir+") needs to be a simulation")
                sys.stdout.flush()
            return -1
    #
    lremove_old = lremove_old_snapshots or\
                  lremove_old_slices or lremove_old_averages
    if lremove_old:
        if not execute:
            os.chdir(olddir)
            if rank == 0:
                print("WARNING: Are you sure you wish to remove the Fortran"+
                      " binary files from \n"+
                      os.getcwd()+".\n"+
                      "Set execute=True to proceed.")
                sys.stdout.flush()
            return -1

    os.chdir(olddir)
    if lvars:
        if varfile_names == None:
            os.chdir(fromdatadir+'/proc0')
            lVARd = False
            varfiled_names = glob.glob('VARd*')
            if len(varfiled_names) > 0:
                varfile_names = glob.glob('VAR*')
                for iv in range(len(varfile_names)-1,-1,-1):
                    if 'VARd' in varfile_names[iv]:
                        varfile_names.remove(varfile_names[iv])
                lVARd = True
            else:
                varfile_names = glob.glob('VAR*')
            os.chdir(olddir)
        else:
            lVARd = False
            if isinstance(varfile_names, list):
                varfile_names = varfile_names
            else:
                varfile_names = [varfile_names]
            varfiled_names = []
            tmp_names = []
            for varfile_name in varfile_names:
                if 'VARd' in varfile_names:
                    varfiled_names.append(varfile_name)
                    lVARd = True
                else:
                    tmp_names.append(varfile_name)
            varfile_names = tmp_names
    gkeys = ['x', 'y', 'z', 'Lx', 'Ly', 'Lz', 'dx', 'dy', 'dz',
             'dx_1', 'dy_1', 'dz_1', 'dx_tilde', 'dy_tilde', 'dz_tilde',
            ]
    grid = None
    if rank == size-1:
        grid = read.grid(quiet=True)
    if l_mpi:
        grid=comm.bcast(grid, root=size-1)
    if not quiet:
        print(rank,grid)
        sys.stdout.flush()
    for key in gkeys:
        if not key in grid.__dict__.keys():
            if rank == 0:
                print("ERROR: key "+key+" missing from grid")
                sys.stdout.flush()
            return -1
    #obtain the settings from the old simulation
    settings={}
    skeys = ['l1', 'l2', 'm1', 'm2', 'n1', 'n2',
             'nx', 'ny', 'nz', 'mx', 'my', 'mz',
             'nprocx', 'nprocy', 'nprocz',
             'maux', 'mglobal', 'mvar', 'precision',
            ]
    if rank == 0:
        olddim = read.dim()
        for key in skeys:
            settings[key]=np.array(olddim.__getattribute__(key))
        olddim = None
        settings['nghost']=np.array(nghost)
        settings['precision']=precision.encode()
    if l_mpi:
        settings=comm.bcast(settings, root=0)
    if snap_by_proc:
        nprocs = settings['nprocx']*settings['nprocy']*settings['nprocz']
        if np.mod(nprocs,size) != 0:
            print("WARNING: efficiency requires cpus to divide ncpus")
            sys.stdout.flush()
    if not quiet:
        print(rank,grid)
        sys.stdout.flush()
    #obtain physical units from old simulation
    ukeys = ['length', 'velocity', 'density', 'magnetic', 'time',
                 'temperature', 'flux', 'energy', 'mass', 'system',
                ]
    param = read.param(quiet=True)
    param.__setattr__('unit_mass',param.unit_density*param.unit_length**3)
    param.__setattr__('unit_energy',param.unit_mass*param.unit_velocity**2)
    param.__setattr__('unit_time',param.unit_length/param.unit_velocity)
    param.__setattr__('unit_flux',param.unit_mass/param.unit_time**3)
    param.unit_system=param.unit_system.encode()
    #index list for variables in f-array
    if not quiet:
        print(rank,param)
        sys.stdout.flush()
    indx = None
    if rank == 0:
        indx = read.index()
    if l_mpi:
        indx=comm.bcast(indx, root=0)

    #check consistency between Fortran binary and h5 data
    os.chdir(newdir)
    dim = None
    if is_sim_dir():
        if rank == size-1:
            if exists(join(newdir,'data','dim.dat')):
                try:
                    dim = read.dim()
                except ValueError:
                    pass
        if l_mpi:
            dim=comm.bcast(dim, root=size-1)
        if dim:
            if not quiet:
                print(rank,dim)
                sys.stdout.flush()
            try:
                dim.mvar == settings['mvar']
                dim.mx   == settings['mx']
                dim.my   == settings['my']
                dim.mz   == settings['mz']
            except ValueError:
                if rank == size-1:
                    print("ERROR: new simulation dimensions do not match.")
                    sys.stdout.flush()
                return -1
            dim = None
    os.chdir(olddir)
    if rank == size-1:
        print('precision is ',precision)
        sys.stdout.flush()
    if laver2D:
        print('fred: rank {} laver2D'.format(rank))
        aver2h5(newdir, olddir,
                todatadir='data/averages', fromdatadir='data', l2D=False,
                precision=precision, quiet=quiet, laver2D=laver2D,
                lremove_old_averages=False, aver_by_proc=aver_by_proc,
                l_mpi=l_mpi, driver=driver, comm=comm, rank=rank, size=size)
        l2D = False
    #copy snapshots
    if lvars and len(varfile_names) > 0:
        var2h5(newdir, olddir, varfile_names, todatadir, fromdatadir,
               snap_by_proc, precision, lpersist, quiet, nghost, settings,
               param, grid, x, y, z, lshear, lremove_old_snapshots, indx,
               l_mpi=l_mpi, driver=driver, comm=comm, rank=rank, size=size)
    #copy downsampled snapshots if present
    if lvars and lVARd:
        print('fred: rank {} lVARd'.format(rank))
        var2h5(newdir, olddir, varfiled_names, todatadir, fromdatadir,
               False, precision, lpersist, quiet, nghost, settings,
               param, grid, x, y, z, lshear, lremove_old_snapshots, indx,
               trimall=True, l_mpi=l_mpi,
               driver=driver, comm=comm, rank=rank, size=size)
    if lvars:
        print('fred: rank {} lvars'.format(rank))
        var2h5(newdir, olddir, ['var.dat',], todatadir, fromdatadir,
               snap_by_proc,
               precision, lpersist, quiet, nghost, settings, param, grid,
               x, y, z, lshear, lremove_old_snapshots, indx,
               l_mpi=l_mpi, driver=driver, comm=comm, rank=rank, size=size)
    #copy old video slices to new h5 sim
    if lvids:
        print('fred: rank {} lvids'.format(rank))
        if lremove_deprecated_vids:
            for ext in ['bb.','uu.','ux.','uy.','uz.','bx.','by.','bz.']:
                cmd = 'rm -f '+join(fromdatadir,'proc*','slice_'+ext+'*')
                if rank == 0:
                    process = sub.Popen(cmd.split(),stdout=sub.PIPE)
                    output, error = process.communicate()
                    print(cmd,output,error)
                    #os.system(cmd)
                cmd = 'rm -f '+join(fromdatadir,'slice_'+ext+'*')
                if rank == 0:
                    process = sub.Popen(cmd.split(),stdout=sub.PIPE)
                    output, error = process.communicate()
                    print(cmd,output,error)
                    #os.system(cmd)
        if comm:
            comm.Barrier()
        cmd = 'src/read_all_videofiles.x'
        if rank == size-1 and lread_all_videoslices:
            process = sub.Popen(cmd.split(),stdout=sub.PIPE)
            output, error = process.communicate()
            print(cmd,output,error)
            #os.system(cmd)
        if comm:
            comm.Barrier()
        print('fred: rank {} slice'.format(rank))
        slices2h5(newdir, olddir, grid,
                  todatadir='data/slices', fromdatadir=fromdatadir,
                  precision=precision, quiet=quiet, vlarge=vlarge,
                  lsplit_slices=lsplit_slices,
                  lremove_old_slices=lremove_old_slices, l_mpi=l_mpi,
                  driver=driver, comm=comm, rank=rank, size=size)
    #copy old averages data to new h5 sim
    if laver:
        print('fred: rank {} aver 1D'.format(rank))
        aver2h5(newdir, olddir,
                todatadir='data/averages', fromdatadir=fromdatadir, l2D=l2D,
                precision=precision, quiet=quiet, aver_by_proc=False,
                lremove_old_averages=lremove_old_averages, l_mpi=l_mpi,
                driver=driver, comm=comm, rank=rank, size=size)
    #check some critical sim files are present for new sim without start
    #construct grid.h5 sim information if requied for new h5 sim
    os.chdir(newdir)
    if l_mpi:
        comm.Barrier()
    if rank == 0:
        write_h5_grid(file_name='grid', datadir='data', precision=precision,
                      nghost=nghost, settings=settings, param=param, grid=grid,
                      unit=None, quiet=quiet)
        source_file = join(olddir,fromdatadir,'proc0/varN.list')
        target_file = join(newdir,todatadir,'varN.list')
        if exists(source_file):
            cmd='cp '+source_file+' '+target_file
            process = sub.Popen(cmd.split(),stdout=sub.PIPE)
            output, error = process.communicate()
            print(cmd,output,error)
            #os.system(cmd)
        items=['def_var.pro', 'index.pro', 'jobid.dat', 'param.nml',
               'particle_index.pro', 'pc_constants.pro', 'pointmass_index.pro',
               'pt_positions.dat', 'sn_series.dat', 'svnid.dat',
               'time_series.dat', 'tsnap.dat', 'tspec.dat', 'tvid.dat',
               't2davg.dat', 'var.general', 'variables.pro', 'varname.dat']
        for item in items:
            source_file = join(olddir,fromdatadir,item)
            target_file = join(newdir,fromdatadir,item)
            if exists(source_file):
                if not exists(target_file):
                    cmd='cp '+source_file+' '+target_file
                    process = sub.Popen(cmd.split(),stdout=sub.PIPE)
                    output, error = process.communicate()
                    print(cmd,output,error)
                    #os.system(cmd)
    print('Simulation Fortran to h5 completed on rank {}.'.format(rank))
    sys.stdout.flush()
