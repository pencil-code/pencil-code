"""
Read data and parameters from data directories.
"""

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
    comm = None
    l_mpi = False

try:
    import f90nml
    from .clonesimulations import clone_sims
    from .clonesimulations import clone_sims_from_obj
    from .makedict import parameter_table, make_sims_dict
except ModuleNotFoundError as e:
    if rank == 0:
        print("Warning: Module required for pipelines not found.", e)

