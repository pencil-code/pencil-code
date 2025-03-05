"""
Visualization routines.
"""

from pencil.util import pc_print

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
#    from . import rvid_box_vista
except Exception as e:
    pc_print(f"Exception while loading PyVista plotter tools. Exception: {e}")
    pc_print(
        "Warning: Make sure you have all the required libraries! See pv_plotter.py"
        " docstrings for requirements."
    )

try:
    from . import rvid_box
except:
    pc_print("Warning: Could not import visu.rvid_box. Try:")
    pc_print("$ conda install -c plotly plotly-orca psutil requests")
