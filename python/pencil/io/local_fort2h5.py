# local_fort2h5.py
#
# Call read existing Fortran unformatted simulation data and write as hdf5.
#
#
# Author: F. Gent (fred.gent.ncl@gmail.com).
#
"""
Example script to call from the comamnd line or from a batch file to execute the
copying of an existing simulation using Fortran binary data into the same or new
simulation using hdf5: farray data, etc to data/allprocs/VAR*.h5 files, video
slices to data/slices/uu1_xy.h5, etc. and averages to data/averages/xy.h5,
and necessary auxilliary file to continue the run in the new format.
Copy this file to the local drive and edit as required.

To apply from the command line:
> python local_fort2h5.py

"""

import pencil as pc

pc.io.sim2h5(
    newdir="full/path/to/fortran/binary/simulation",
    olddir="full/path/to/hdf5/copy/of/simulation",  # can be the same
    # if 2D averages files are very large an alternative reading method
    # is required so set True
    laver2D=False,
    # if var files not required or have already been converted set False
    lvars=True,
    # if videos not required or have already been converted set False
    lvids=True,
    # if averages not required or have already been converted set False
    laver=True,
    # if averages include 2D and not laver2D set True
    l2D=False,
    # CAUTION old files can be deleted during progress
    lremove_old_snapshots=True,
    lremove_old_slices=True,
    lremove_old_averages=True,
    # warning given before executing removal of data until set True
    execute=False,
)
