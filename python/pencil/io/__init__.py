"""
Input und output functions, like saving data.
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

from .npfile import npfile
from .mkdir import mkdir
from .debug_breakpoint import debug_breakpoint
from .timestamp import timestamp
from .pc_hdf5 import *
from .snapshot import *

try:
    from .fort2h5 import *
except:
    if rank == 0:
        print("Warning: Could not import io.fort2h5. Try:")
        print("'pip3 install h5py' (Python 3) or 'pip install h5py' (Python 2).")

# io operation on cluster/computer
from .get_systemid import get_systemid
from .exists_file import exists_file
from .remove_files import remove_files

# io operation on simulation
from .get_value_from_file import get_value_from_file
from .change_value_in_file import change_value_in_file
from .rename_in_submit_script import rename_in_submit_script

# dill im-/exporter
try:
    from .dill_load import dill_load as load
    from .dill_save import dill_save as save
    from .dill_exists import dill_exists as exists
except:
    if rank == 0:
        print("Warning: Could not import io.dill*. Try:")
        print("'pip3 install dill' (Python 3) or 'pip install dill' (Python 2).")

# pkl im-/exporter
from .pkl_load import pkl_load  # as load
from .pkl_save import pkl_save  # as save
from .pkl_exists import pkl_exists  # as exists
from .walklevel import walklevel
