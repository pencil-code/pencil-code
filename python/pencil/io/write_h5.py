# write_h5.py
#
# Generate a snapshot or initial condition from a numpy array.
#
#
# Author: F. Gent (fred.gent.ncl@gmail.com).
#
"""wrappers for creating h5 files and groups or datasets within, to avoid 
   overwriting in error exisitng files or missing parents.

"""
from . import mkdir
from os.path import exists, join
import subprocess as sub
import h5py

def open_h5(filename, mode='r', overwrite=False):
    if '/' in filename:
        fname = str.split(filename,'/')[-1]
        path = str.strip(filename,fname)
    else:
        fname = filename
        path = './'
    if not ('.h5' == filename[-3:] or '.hdf5' == filename[-5:]):
        print('Relabelling h5 '+fname+' to '+fname+'.h5 on path '+path)
        fname = fname+'.h5'
    mkdir(path)
    if exists(join(path,fname)):
        if (mode == 'w' or mode=='a' or mode=='r+') and not overwrite:
            cmd = 'mv '+join(path,fname)+' '+join(path,fname+'.bak')
            process = sub.Popen(cmd.split(),stdout=sub.PIPE)
            output, error = process.communicate()
            print(cmd,output,error)
    return h5py.File(join(path,fname), mode)

def group_h5(h5obj, groupname, mode='r', delete=False):
    if not h5obj.__contains__(groupname):
        if mode == 'r':
            print('group_h5: '+h5obj.filename+' does not contain '+groupname)
            return 0
        else:
            h5obj.create_group(groupname)
    else:
        if not mode == 'r' and delete:
            h5obj.__delitem__(groupname)
            print('group_h5: '+groupname+' deleted from '+h5obj.filename)
            return 0
    return h5obj[groupname]

def dataset_h5(h5obj, dataname, mode='r', data=None, shape=None, dtype=None,
               overwrite=False, delete=False):
    try:
        ldata = len(data)>0
        lshape = len(shape)>0
    except:
        ldata = data is not None    
        lshape = shape is not None
    if not h5obj.__contains__(dataname):
        if mode == 'r':
            print(h5obj.name+' does not contain '+dataname)
            return 0
        else:
            
            if not ldata:
                if not lshape:
                    print('dataset_h5: data or shape must be provided')
                elif not dtype:
                    print('dataset_h5: data not present, provide dtype')
                else:
                    h5obj.create_dataset(dataname, (shape,), dtype=dtype)
            else:
                if not dtype:
                    h5obj.create_dataset(dataname, data=data)
                else:
                    h5obj.create_dataset(dataname, data=data, dtype=dtype)
    else:
        if not mode == 'r':
            if delete:
                h5obj.__delitem__(dataname)
                print('dataset_h5: '+dataname+' deleted from '+h5obj.name)
                return 0
            if overwrite:
                h5obj.__delitem__(dataname)
                if ldata:
                    if lshape:
                        print('dataset_h5: data or shape must be provided')
                    elif not dtype:
                        print('dataset_h5: data not present, provide dtype')
                    else:
                        h5obj.create_dataset(dataname, (shape,), dtype=dtype)
                else:
                    if not dtype:
                        h5obj.create_dataset(dataname, data=data)
                    else:
                        h5obj.create_dataset(dataname, data=data, dtype=dtype)
    return h5obj[dataname]
