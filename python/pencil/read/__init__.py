'''
Reading routines for the pencil code data.
'''

from .pdims import pdim, PDim
from .indices import index
from .timeseries import ts
from .powers import power
from .ogdims import ogdim
from .ogvar import ogvar
from .pencil_record_types import record_types
from .dims import dim
from .params import param
from .grids import grid
from .varfile import var
from .allslices import slices
from .averages import aver
from .pvarfile import pvar
from .phiaverages import phiaver

# idl workarounds
from .pstalk import pstalk
