'''
Visualization routine.
'''

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
    from . import rvid_box
except:
    print('Warning: Could not import visu.rvid_box. Try:')
    print('$ conda install -c plotly plotly-orca psutil requests')
