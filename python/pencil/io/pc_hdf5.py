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
import numpy as np
from pencil.io import mkdir
from os.path import exists, join
from os import remove
import subprocess as sub

#==============================================================================
def open_h5(filename, status, driver=None, comm=None, overwrite=False,
            size=1, rank=0, lfs=False, MB=1, count=1):
    """This function opens hdf5 file in serial or parallel.

    Keyword arguments:
        filename:  relative or absolute path string for name of hdf5 file.
        status:    open state of file 'w': write, 'r': read or 'a'/'r+': append.
        driver:    'mpio' required for parallel: version but absent for serial.
        comm:      only present for parallel version of h5py.
        overwrite: flag to replace existing file.
        rank:      processor rank with root = 0.
    """
    if '/' in filename:
        fname = filename.split('/')[-1]
        path = filename.split(fname)[0]
    else:
        fname = filename
        path = './'
    if not ('.h5' == filename[-3:] or '.hdf5' == filename[-5:]):
        if np.mod(rank,size) == 0:
            print('Relabelling h5 '+fname+' to '+
                   str.strip(fname,'.dat')+'.h5 on path '+path)
        fname = str.strip(fname,'.dat')+'.h5'
    mkdir(path, rank=rank, lfs=lfs, MB=MB, count=count)
    if status == 'w' and exists(join(path,fname)):
        if not overwrite:
            try:
                cmd = 'mv '+join(path,fname)+' '+join(path,fname+'.bak')
                process = sub.Popen(cmd.split(),stdout=sub.PIPE)
                output, error = process.communicate()
                print(cmd,output,error)
            except:
                while exists(join(path,fname)):
                    pass
    if comm:
        if not driver:
            driver = 'mpio'
        dset = h5py.File(join(path,fname), status, driver=driver, comm=comm)
    else:
        dset = h5py.File(join(path,fname), status)

    return dset

#==============================================================================
def group_h5(h5obj, groupname, status='r', delete=False, overwrite=False,
             rank=0, size=1, comm=None):
    """This function adds/removes hdf5 group objects.

    Keyword arguments:
        h5obj:     h5 object, may be the file or a sub group within the file
        groupname: string for name of the group.
        status:    open state of file 'w': write, 'r': read or 'a'/'r+': append.
        delete:    flag to remove existing group from h5 object.
        overwrite: flag to replace existing group from h5 object.
        rank:      processor rank with root = 0.
        comm:      only present for parallel version of h5py.
    """
    #if both overwrite and delete, delete is False
    if delete:
        delete = not overwrite
    if not h5obj.__contains__(groupname):
        if status == 'r':
            if np.mod(rank,size) == 0:
                print('group_h5: '+h5obj.filename+' does not contain '+groupname)
            return False
        else:
            h5obj.create_group(groupname)
    else:
        if not status == 'r' and (delete or overwrite):
            try:
                h5obj.__delitem__(groupname)
            except:
                while h5obj.__contains__(groupname):
                    pass
            if np.mod(rank,size) == 0:
                 print('group_h5: '+groupname+' deleted from '+h5obj.filename)
            if not delete:
                h5obj.create_group(groupname)
                if np.mod(rank,size) == 0:
                    print('group_h5: '+groupname+' replaced in '+h5obj.filename)
            else:
                return False
    return h5obj[groupname]

#==============================================================================
def dataset_h5(h5obj, dataname, status='r', data=None, shape=None, dtype=None,
               overwrite=False, delete=False, rank=0, size=1, comm=None):
    """This function adds/removes hdf5 dataset objects.

    Keyword arguments:
        h5obj:     h5 object, may be the file or a sub group within the file
        dataname:  string for name of the dataset.
        status:    open state of file 'w': write, 'r': read or 'a'/'r+': append.
        data:      h5 compatible data object; float, integer, string, array
        shape:     data shape tuple of length > 0
        dtype:     h5 compatible data type, eg. np.float64
        delete:    flag to remove existing group from h5 object.
        overwrite: flag to replace existing group from h5 object.
        rank:      processor rank with root = 0.
        comm:      only present for parallel version of h5py.
    """
    try:
        ldata = len(data)>0
    except:
        ldata = data is not None
    try:
        lshape = len(shape)>0
    except:
        lshape = shape is not None
    #if both overwrite and delete, delete is False
    if delete:
        delete = not overwrite
    if not h5obj.__contains__(dataname):
        if status == 'r':
            if np.mod(rank,size) == 0:
                print(h5obj.name+' does not contain '+dataname)
            return False
        else:
            if not ldata:
                if not lshape:
                    if np.mod(rank,size) == 0:
                        print('dataset_h5: data or shape must be provided')
                elif not dtype:
                    if np.mod(rank,size) == 0:
                        print('dataset_h5: data not present, provide dtype')
                else:
                    h5obj.create_dataset(dataname, shape, dtype=dtype)
            else:
                if not dtype:
                    h5obj.create_dataset(dataname, data=data)
                else:
                    h5obj.create_dataset(dataname, data=data, dtype=dtype)
    else:
        if not status == 'r' and (delete or overwrite):
            try:
                h5obj.__delitem__(dataname)
            except:
                while h5obj.__contains__(dataname):
                    pass
            if np.mod(rank,size) == 0:
                print('dataset_h5: '+dataname+' deleted from '+h5obj.name)
            if not delete:
                if not ldata:
                    if not lshape:
                        if np.mod(rank,size) == 0:
                            print('dataset_h5: data or shape must be provided')
                    elif not dtype:
                        if np.mod(rank,size) == 0:
                            print('dataset_h5: data not present, provide dtype')
                    else:
                        h5obj.create_dataset(dataname, shape, dtype=dtype)
                else:
                    if not dtype:
                        h5obj.create_dataset(dataname, data=data)
                    else:
                        h5obj.create_dataset(dataname, data=data, dtype=dtype)
            else:
                return False
    return h5obj[dataname]
