# ogdim.py
#
# Read the dimensions of the simulationp.
# This routine is a simple adaptation of dim.py used to read dim files
#
# Adapted from dim.py by:
# J. Aarnes (jorgenaarnes@gmail.com)
"""
Contains the classes and methods to read the simulation dimensions.
"""

from .dims import Dim

def ogdim(*args, **kwargs):
    """
    Read the ogdim.dat file.

    call signature:

    ogdim(datadir='data', proc=-1)

    Keyword arguments:

    *datadir*:
      Directory where the data is stored.

    *proc*
      Processor to be read. If proc is -1, then read the 'global'
      dimensions. If proc is >=0, then read the dim.dat in the
      corresponding processor directory.
    """

    ogdim_tmp = ogDim()
    ogdim_tmp.read(*args, ogrid=True, **kwargs)
    return ogdim_tmp


class ogDim(Dim):
    """
    ogDim -- holds pencil code dimension data.
    """

    def __init__(self):
        super(ogDim, self).__init__()
