"""
Math functions and further calculations.
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

from .streamlines import *

####### working with particles and grid data
from .part_to_grid import *  # bin particle quantities to a grid
from .fill_gaps_in_grid import *
from .accuracy import *
from .draglift import *
from .tensors import *
from .Reynolds import *
from .shocktube import sod
from .Gaussian_averages import kernel_smooth, gauss_3Dsmooth 

try:
    from .aver2h5 import *
except:
    if rank == 0:
        print("Warning: Could not import calc.aver2h5. Try:")
        print("'pip3 install h5py' (Python 3) or 'pip install h5py' (Python 2).")
