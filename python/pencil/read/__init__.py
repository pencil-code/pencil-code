'''
Reading routines for the pencil code data.
'''

from .pdim import pdim, PDim
from .index import index
from .ts import ts
from .power import power
from .ogdim import ogdim
from .ogvar import ogvar
from .pencil_record_types import record_types
try:
    from .dim import dim
    from .param import param
    from .grid import grid
    from .var import var
    from .slices import slices
    from .averages import aver
    from .pvar import pvar
    from .phiaverages import phiaver
except:
    print('Warning: Could not import io.pc_hdf5. Try:')
    print('$ pip install h5py')

# idl workarounds
from .pstalk import pstalk
