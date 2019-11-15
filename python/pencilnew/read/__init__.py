'''
Reading routines for the pencil code data.
'''

from .dim import dim
from .pdim import pdim, __Pdim__
from .param import param
from .grid import grid
from .index import index
from .var import var
from .slices import slices
from .ts import ts
from .averages import aver
from .power import power
from .ogdim import ogdim
from .ogvar import ogvar
from .pencil_record_types import record_types

# idl workarounds
from .pstalk import pstalk
from .pvar import pvar
