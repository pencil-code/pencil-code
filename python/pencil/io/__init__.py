################
##
##	io
##
################

from .npfile import npfile
from .mkdir import mkdir
from .debug_breakpoint import debug_breakpoint
from .timestamp import timestamp
try:
    from .pc_hdf5 import *
    from .snapshot import *
    from .fort2h5 import *
except:
    print('Warning: Could not import io.pc_hdf5. Try:')
    print('$ pip install h5py')
    from .snapshot import write_snapshot

# io operation on cluster/computer
from .get_systemid import get_systemid
from .exists_file import exists_file
from .remove_files import remove_files

# io operation on simulation
from .get_value_from_file import get_value_from_file
from .change_value_in_file import change_value_in_file
from .rename_in_submit_script import rename_in_submit_script

# dill im-/exporter
from .dill_load import dill_load as load
from .dill_save import dill_save as save
from .dill_exists import dill_exists as exists

# pkl im-/exporter
from .pkl_load import pkl_load #as load
from .pkl_save import pkl_save #as save
from .pkl_exists import pkl_exists #as exists
from .walklevel import walklevel
