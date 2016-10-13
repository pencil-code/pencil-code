# slices.py
#
# Read the slice files.
#
# Author: S. Candelaresi (iomsn1@gmail.com).
"""
Contains the classes and methods to read slice files.
"""


def slices(*args, **kwargs):
    """
    Read Pencil Code slices data.

    call signature:

    ts(file_name='time_series.dat', data_dir='data',
       quiet=False, comment_char='#')

    Keyword arguments:

    *file_name*:
      Name of the time series file.

    *data_dir*:
      Directory where the data is stored.

    *quiet*
      Flag for switching off output.

    *comment_char*
      Comment character in the time series file.
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

        import numpy as np

        self.t = np.array([])
        self.slices = []


    def read(self, field='uu1', data_dir='data', proc=-1, extension='xz',
             old_file=False, precision='f'):
        """
        Read Pencil Code time series data.

        call signature:

        read(self. field='uu1', data_dir='data', proc=-1, extension='xz',
             old_file=False)

        Keyword arguments:

        *field*:
          Name of the field to be read.

        *data_dir*:
          Directory where the data is stored.

        *proc*:
          Processor to be read. If -1 read all and assemble to one array.

        *extension*
          Specifies the slice.

        *old_file*
          Flag for reading old file format.

        *precision*:
          Precision of the data. Either float 'f' or double 'd'.
        """

        import os.path
        import numpy as np
        from scipy.io import FortranFile
        from pencilnew import read

        data_dir = os.path.expanduser(data_dir)
        if proc < 0:
            file_name = os.path.join(data_dir, 'slice_'+field+'.'+extension)
        else:
            file_name = os.path.join(data_dir, 'proc{0}'.format(proc),
                                     'slice_'+field+'.'+extension)

        dim = read.dim(data_dir, proc)
        if dim.precision == 'D':
            precision = 'd'
        else:
            precision = 'f'

        # Set up slice plane.
        if (extension == 'xy' or extension == 'Xy'):
            hsize = dim.nx
            vsize = dim.ny
        if (extension == 'xz'):
            hsize = dim.nx
            vsize = dim.nz
        if (extension == 'yz'):
            hsize = dim.ny
            vsize = dim.nz

        infile = FortranFile(file_name)

        islice = 0
        self.t = np.zeros(1, dtype=precision)
        self.slices = np.zeros(1, dtype=precision)

        while True:
            try:
                raw_data = infile.read_record(dtype=precision)
            except ValueError:
                break
            except TypeError:
                break

            if old_file:
                self.t = np.concatenate((self.t, raw_data[-1:]))
                self.slices = np.concatenate((self.slices, raw_data[:-1]))
            else:
                self.t = np.concatenate((self.t, raw_data[-2:-1]))
                self.slices = np.concatenate((self.slices, raw_data[:-2]))
            islice += 1

        # Reshape and remove first entry.
        self.t = self.t[1:]
        self.slices = self.slices[1:].reshape(islice, vsize, hsize)
