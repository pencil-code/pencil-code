# allslices.py
#
# Read the slice files.
#
# Author: S. Candelaresi (iomsn1@gmail.com).
"""
Contains the classes and methods to read slice files.
"""

import numpy as np

def slices(*args, **kwargs):
    """
    Read Pencil Code slice data.

    Signature:

    read(field='', extension='', datadir='data', proc=-1, old_file=False,
         precision='f', iter_list=list(), vlarge=1000000000, quiet=True,
         tstart=0, tend=None)

    Parameters
    ----------
    *field*: Name of the field(s) to be read (string or string list).

    *extension*: Specifies the plane slice(s).

    *datadir*: Directory where the data is stored.

    *proc*: Processor to be read. If -1 read all and assemble to one array.

    *old_file*: Flag for reading old file format.

    *precision*: Precision of the data. Either float 'f' or double 'd'.

    *iter_list*: List of iteration indices for which to sample the slices.

    *vlarge*: Size of arrays over which to apply memory limits reshaping.

    *quiet*: Flag for switching off output.

    *tstart*: Start time interval from which to sample slices.

    *tend*: End time interval from which to sample slices.

    Returns 
    ------- 
    Class containing the fields and slices as attributes. 

    Notes
    -----
    Use the attribute keys to get a list of attributes

    Examples
    --------
    >>> vsl = pc.read.slices()
    >>> vsl.keys()
    t
    xy
    xy2
    xz
    yz
    position
    coordinate
    """

    slices_tmp = SliceSeries()
    slices_tmp.read(*args, **kwargs)
    return slices_tmp


class SliceSeries(object):
    """
    SliceSeries -- holds Pencil Code slices data and methods.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.t = np.array([])

    def keys(self):
        for i in self.__dict__.keys():
            print(i)

    def read(self, field='', extension='', datadir='data', proc=-1,
             old_file=False, precision='f',
             iter_list=list(), vlarge=1000000000,
             quiet=True, tstart=0, tend=None
            ):
        """
        Read Pencil Code slice data.

        call signature:

        read(field='', extension='', datadir='data', proc=-1, old_file=False,
             precision='f', iter_list=list(), vlarge=1000000000, quiet=True,
             tstart=0, tend=None)

        Keyword arguments:

        *field*: Name of the field(s) to be read (string or string list).

        *extension*: Specifies the plane slice(s).

        *datadir*: Directory where the data is stored.

        *proc*: Processor to be read. If -1 read all and assemble to one array.

        *old_file*: Flag for reading old file format.

        *precision*: Precision of the data. Either float 'f' or double 'd'.

        *iter_list*: List of iteration indices for which to sample the slices.

        *vlarge*: Size of arrays over which to apply memory limits reshaping.

        *quiet*: Flag for switching off output.

        *tstart*: Start time interval from which to sample slices.

        *tend*: End time interval from which to sample slices.

        """

        import os
        import sys
        import numpy as np
        from scipy.io import FortranFile
        from pencil import read

        if os.path.exists(os.path.join(datadir, 'grid.h5')):
            l_h5 = True
            import h5py
        else:
            l_h5 = False

        if l_h5:
            # Define the directory that contains the slice files.
            slice_dir = os.path.join(datadir, 'slices')
            # Initialize the fields list.
            if field:
                if isinstance(field, list):
                    field_list = field
                else:
                    field_list = [field]
            else:
                # Find the existing fields.
                field_list = []
                for file_name in os.listdir(slice_dir):
                    field_list.append(file_name.split('_')[0])
                # Remove duplicates.
                field_list = list(set(field_list))
            # Initialize the extensions list.
            if extension:
                if isinstance(extension, list):
                    extension_list = extension
                else:
                    extension_list = [extension]
            else:
                # Find the existing extensions.
                extension_list = []
                for file_name in os.listdir(slice_dir):
                    extension_list.append(
                              file_name.split('_')[1].split('.')[0])
                # Remove duplicates.
                extension_list = list(set(extension_list))

            class Foo(object):
                pass

            if len(iter_list) > 0:
                nt = len(iter_list)
                if tstart > 0 or tend:
                    print('read.slices: using iter_list.',
                          'If tstart or tend required set iter_list=None')
                tstart = 0; tend = None
            else:
                nt = None
            pos_object = Foo()
            ind_object = Foo()
            for extension in extension_list:
                if not quiet:
                    print('Extension: ' + str(extension))
                    sys.stdout.flush()
                # This one will store the data.
                ext_object = Foo()
                pos_list=[]
                ind_list=[]
                for field in field_list:
                    if not quiet:
                        print('  -> Field: ' + str(field))
                        sys.stdout.flush()
                    #Compose the file name according to field & extension.
                    file_name = os.path.join(slice_dir,
                                             field+'_'+extension+'.h5')
                    with h5py.File(file_name, 'r') as ds:
                        if tstart > 0 or tend:
                            iter_list = []; it = 0
                            while it < ds['last'][0]:
                                if ds[str(it+1)+'/time'][()] > tstart:
                                    iter_list.append(it)
                                    it += 1
                                if ds[str(it+1)+'/time'][()] > tend:
                                    it = ds['last'][0]
                                print('iter_list: it={}, time={}'.format(
                                      it,ds[str(it+1)+'/time'][()]))
                            nt = len(iter_list)
                        if not isinstance(nt, int):
                            nt = ds['last'][0]
                        vsize = ds['1/data'].shape[0]
                        hsize = ds['1/data'].shape[1]
                        slice_series = np.zeros([nt,vsize,hsize],dtype=precision)
                        if len(iter_list) == 0:
                            iter_list = list(range(nt))
                        for it in iter_list:
                            if ds.__contains__(str(it+1)):
                                slice_series[it] = ds[str(it+1)+'/data'][()]
                            else:
                                print('no data at {} in '.format(it+1)+
                                      file_name)
                        add_pos = len(pos_list)==0
                        if self.t.size == 0:
                            self.t = []
                            for it in iter_list:
                                self.t.append(ds[str(it+1)+'/time'][()])
                                if add_pos:
                                    ind_list.append(ds[str(it+1)+'/coordinate'][0])
                                    pos_list.append(ds[str(it+1)+'/position'][()])
                            self.t = np.array(self.t).astype(precision)
                            setattr(pos_object, extension, np.array(pos_list))
                            setattr(ind_object, extension, np.array(ind_list))
                        else:
                            if add_pos:
                                for it in iter_list:
                                    ind_list.append(ds[str(it+1)+'/coordinate'][0])
                                    pos_list.append(ds[str(it+1)+'/position'][()])
                                setattr(pos_object, extension, np.array(pos_list))
                                setattr(ind_object, extension, np.array(ind_list))
                    setattr(ext_object, field, slice_series)

                setattr(self, extension, ext_object)
                setattr(self, 'position', pos_object)
                setattr(self, 'coordinate', ind_object)
        else:
            # Define the directory that contains the slice files.
            if proc < 0:
                slice_dir = datadir
            else:
                slice_dir = os.path.join(datadir, 'proc{0}'.format(proc))

            # Initialize the fields list.
            if field:
                if isinstance(field, list):
                    field_list = field
                else:
                    field_list = [field]
            else:
                # Find the existing fields.
                field_list = []
                for file_name in os.listdir(slice_dir):
                    if file_name[:6] == 'slice_':
                        field_list.append(file_name.split('.')[0][6:])
                # Remove duplicates.
                field_list = list(set(field_list))
                try:
                    field_list.remove('position')
                except:
                    pass

            # Initialize the extensions list.
            if extension:
                if isinstance(extension, list):
                    extension_list = extension
                else:
                    extension_list = [extension]
            else:
                # Find the existing extensions.
                extension_list = []
                for file_name in os.listdir(slice_dir):
                    if file_name[:6] == 'slice_':
                        extension_list.append(file_name.split('.')[1])
                # Remove duplicates.
                extension_list = list(set(extension_list))
                try:
                    extension_list.remove('dat')
                except:
                    pass

            class Foo(object):
                pass

            for extension in extension_list:
                if not quiet:
                    print('Extension: ' + str(extension))
                    sys.stdout.flush()
                # This one will store the data.
                ext_object = Foo()

                for field in field_list:
                    if not quiet:
                        print('  -> Field: ' + str(field))
                        sys.stdout.flush()
                    # Compose the file name according to field and extension.
                    datadir = os.path.expanduser(datadir)
                    if proc < 0:
                        file_name = os.path.join(datadir,
                                                 'slice_'+field+'.'+extension)
                    else:
                        file_name = os.path.join(datadir,
                                                 'proc{0}'.format(proc),
                                                 'slice_'+field+'.'+extension)

                    dim = read.dim(datadir, proc)
                    if dim.precision == 'D':
                        read_precision = 'd'
                    else:
                        read_precision = 'f'

                    # Set up slice plane.
                    if extension == 'xy' or extension == 'Xy' or  extension == 'xy2':
                        hsize = dim.nx
                        vsize = dim.ny
                    if extension == 'xz':
                        hsize = dim.nx
                        vsize = dim.nz
                    if extension == 'yz':
                        hsize = dim.ny
                        vsize = dim.nz

                    try:
                        infile = FortranFile(file_name)
                    except:
                        continue

                    islice = 0
                    self.t = []
                    slice_series = []

                    if not quiet:
                        print('  -> Reading... ')
                        sys.stdout.flush()
                    while True:
                        try:
                            raw_data = infile.read_record(
                                       dtype=read_precision).astype(precision)
                        except ValueError:
                            break
                        except TypeError:
                            break

                        if old_file:
                            self.t.append(list(raw_data[-1]))
                            slice_series.append(raw_data[:-1])
                        else:
                            self.t.append(list(raw_data[-2:-1]))
                            slice_series.append(raw_data[:-2])
                        islice += 1
                    if not quiet:
                        print('  -> Done')
                        sys.stdout.flush()

                    # Remove first entry and reshape.
                    if not quiet:
                        print('Reshaping array')
                        sys.stdout.flush()
                    self.t = np.array(self.t, dtype=precision)[:, 0]
                    slice_series = np.array(slice_series, dtype=precision)
                    # Separate the data in chunks if too big
                    if slice_series.size > vlarge:
                        vscaled = (int(vlarge/vsize/hsize)+1)*vsize*hsize
                        nchunks = int(slice_series.size/vscaled)
                        endslice = vscaled*nchunks
                        subslice_series = np.array_split(slice_series[:endslice+1], nchunks, axis=0)
                        ichunk_size = int(subslice_series[0].size/vsize/hsize)
                        slice_tmp = subslice_series[0].reshape(ichunk_size, vsize, hsize)
                        if not quiet:
                            print("slice_series.shape",slice_series.shape)
                            print(field,"chunk_size",nchunks)
                            print("subslice_series[0].size", subslice_series[0].size)
                            sys.stdout.flush()
                        for subslice in subslice_series:
                            ichunk_size = int(subslice.size/vsize/hsize)
                            subslice = subslice.reshape(ichunk_size, vsize, hsize)
                            slice_tmp = np.concatenate([slice_tmp, subslice], axis=0)
                            if not quiet:
                                print("slice_tmp.shape",slice_tmp.shape)
                                sys.stdout.flush()
                        if not nchunks*vscaled==slice_series.size:
                            ichunk_size = int(slice_series[endslice+1:].size/vsize/hsize)
                            subslice = slice_series[endslice+1:].reshape(ichunk_size, vsize, hsize)
                            slice_tmp = np.concatenate([slice_tmp, subslice], axis=0)
                            if not quiet:
                                print("slice_tmp.shape",slice_tmp.shape)
                                sys.stdout.flush()
                        slice_series = slice_tmp
                    else:
                        slice_series = slice_series.reshape(islice, vsize, hsize)
                    setattr(ext_object, field, slice_series)

                setattr(self, extension, ext_object)
