"""
Visualization routines.
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

from .animate_interactive import animate_interactive

from .animate_multislices import *
from .animate_slices_compareruns import *
from .animate_slices_maketomovie import *
from .animate_slices import *

## line integral convolution
from .lic import *

## general plotting and exporting
from . import internal

try:
    from . import pv_plotter
    from . import pv_plotter_utils
    from . import pv_volume_plotter
except Exception as e:
    if rank == 0:
        print(f"Exception while loading PyVista plotter tools. Exception: {e}")
        print(
            "Warning: Make sure you have all the required libraries! See pv_plotter.py"
            " docstrings for requirements."
        )

try:
    from . import rvid_box
except:
    if rank == 0:
        print("Warning: Could not import visu.rvid_box. Try:")
        print("$ conda install -c plotly plotly-orca psutil requests")
