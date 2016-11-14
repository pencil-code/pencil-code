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
    Read Pencil Code time series data.

    call signature:

    read(self. field='uu1', extension='xz', data_dir='data', proc=-1,
         old_file=False, precision='f')

    Keyword arguments:

    *field*:
      Name of the field(s) to be read.

    *extension*
      Specifies the slice(s).

    *data_dir*:
      Directory where the data is stored.

    *proc*:
      Processor to be read. If -1 read all and assemble to one array.

    *old_file*
      Flag for reading old file format.

    *precision*:
      Precision of the data. Either float 'f' or double 'd'.
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


    def read(self, field='', extension='', data_dir='data', proc=-1,
             old_file=False, precision='f'):
        """
        Read Pencil Code time series data.

        call signature:

        read(self. field='', extension='', data_dir='data', proc=-1,
             old_file=False, precision='f')

        Keyword arguments:

        *field*:
          Name of the field(s) to be read.

        *extension*
          Specifies the slice(s).

        *data_dir*:
          Directory where the data is stored.

        *proc*:
          Processor to be read. If -1 read all and assemble to one array.

        *old_file*
          Flag for reading old file format.

        *precision*:
          Precision of the data. Either float 'f' or double 'd'.
        """

        import os
        import numpy as np
        from scipy.io import FortranFile
        from pencilnew import read

        # Define the directory that contains the slice files.
        if proc < 0:
            slice_dir = data_dir
        else:
            slice_dir = os.path.join(data_dir, 'proc{0}'.format(proc))

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
                if (file_name[:6] == 'slice_'):
                    field_list.append(file_name.split('.')[0][6:])
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
                if (file_name[:6] == 'slice_'):
                    extension_list.append(file_name.split('.')[1])
            # Remove duplicates.
            extension_list = list(set(extension_list))

        class Foo():
            pass

        for extension in extension_list:
            # This one will store the data.
            ext_object = Foo()

            for field in field_list:
                # Compose the file name according to field and extension.
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
                slice_series = np.zeros(1, dtype=precision)

                while True:
                    try:
                        raw_data = infile.read_record(dtype=precision)
                    except ValueError:
                        break
                    except TypeError:
                        break

                    if old_file:
                        self.t = np.concatenate((self.t, raw_data[-1:]))
                        slice_series = np.concatenate((slice_series, raw_data[:-1]))
                    else:
                        self.t = np.concatenate((self.t, raw_data[-2:-1]))
                        slice_series = np.concatenate((slice_series, raw_data[:-2]))
                    islice += 1

                # Reshape and remove first entry.
                self.t = self.t[1:]
                slice_series = slice_series[1:].reshape(islice, vsize, hsize)
                setattr(ext_object, field, slice_series)

            setattr(self, extension, ext_object)

