import pencil as pc
import os
import glob

try: 
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    l_mpi = True
    l_mpi = l_mpi and (size != 1)
except ImportError:
    rank = 0
    size = 1

if comm:
    comm.Barrier()
    if os.path.exists('data/grid.h5'):
        cmd = 'mv data/grid.h5 data/grid.h5.bak'
        os.system(cmd)
pc.io.sim2h5(
             newdir='/users/fagent/pencil-code/python/tutorials/mps/ism_binary_h5',
             laver2D=False,
             lvars=True,
             lvids=True,
             laver=True, 
             lremove_old_snapshots=False,
             lremove_old_slices=False,
             lremove_old_averages=False,
             execute=False,
             lremove_deprecated_vids=True,
             lread_all_videoslices=False,
             vlarge=100000000,
             quiet=True,
            )
if comm:
    comm.Barrier()
if rank == 0:
    cmd = 'make distclean'
    os.system(cmd)
    src_list = glob.glob('./src/*')
    for src in src_list:
        if not 'local' in src:
            cmd = 'rm -rf '+src
            os.system(cmd)
