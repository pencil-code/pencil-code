# pc_hdf5.py
#
# pencil_code python wrappers for hdf5 operations.
#
# Authors:
# F. A. Gent (fred.gent.ncl@gmail.com & frederick.gent@aalto.fi)
# 18/06/2020
#
"""
Contains optional h5 operations with and without hdf5-parallel, for exampe
to open h5 file with/without MPI using common switches.
TODO open h5 pencil var object etc, as alternative to var object for 
large datasets with memory limits.
"""
import h5py

def open_h5(name, status, driver=None, comm=None):
    """This script opens file in serial or parallel.

    Keyword arguments:
        name:   relative or absolute path sting for name of hdf5 file.
        status: state of opened file 'w': write, 'r':read or 'a'/'r+': append.
        driver: 'mpio' required for parallel: version but absent for serial.
        comm:   only present for parallel version of h5py
    """

    if comm:
        if not driver:
            driver = 'mpio'
        dset = h5py.File(name, status, driver=driver, comm=comm)
    else:
        dset = h5py.File(name, status)

    return dset

#==============================================================================

