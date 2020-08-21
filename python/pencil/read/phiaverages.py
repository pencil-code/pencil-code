# phiaverages.py
#
# Read phi-average files.
#
# Adapted from averages.py by P. Kaepylae (pkaepyl@uni-goettingen.de)
"""
Contains the classes and methods to read phi-averaged files.
"""
import sys

def phiaver(*args, **kwargs):
    """
    Read Pencil Code phi-averaged data.

    call signature:

    read(avg_file='1'):

    Keyword arguments:

   *iter_list*
     list of iteration indices for which to sample the slices

   *var_index*:
     Index of single variable taken from among the 'phi' averages.
     Takes an integer value < len(phiaver.in).
    """

    averages_tmp = Averages()
    averages_tmp.read(*args, **kwargs)
    return averages_tmp


class Averages(object):
    """
    Averages -- holds Pencil Code averages data and methods.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        import numpy as np

        self.t = np.array([])


    def read(self, datadir='data', avg_file='1', var_index=-1, iter_list=None):
        """
        Read Pencil Code phi-average data.

        call signature:

        read(avf_file='1'):

        Keyword arguments:

        *iter_list*
          list of iteration indices for which to sample the slices

        *var_index*:
          Index of variable from among within the 'y' or 'z' averages.
          Takes an integer value < len(yaver.in or zaver.in).
        """

        import os

        l_h5 = False

        # Determine which average files to read.
        in_file_name_list = []
        aver_file_name_list = []
        if os.path.exists(os.path.join(datadir, 'grid.h5')):
            l_h5 = True
            #
            # Not implemented
            #
        else:
            if os.path.exists('data/averages/phiavg.list'):
                in_file_name_list.append('data/averages/phiavg.list')
                aver_file_name_list.append(os.path.join('averages', 'PHIAVG'+avg_file))
        if not in_file_name_list:
            print("error: invalid plane name")
            sys.stdout.flush()
            return -1

        class Foo(object):
            pass

        for in_file_name, aver_file_name in \
        zip(in_file_name_list, aver_file_name_list):
            # This one will store the data.
            ext_object = Foo()

            # Get the averaged quantities.
            file_id = open(os.path.join(os.path.dirname(datadir), in_file_name))
            variables = file_id.readlines()
            file_id.close()
            for i in range(sum(list(map(self.__equal_newline, variables)))):
                variables.remove('\n')
            n_vars = len(variables)

            t, r_cyl, z_cyl, raw_data = self.__read_phiaver(datadir, variables,
                                               aver_file_name, n_vars,
                                               var_index, iter_list, l_h5=l_h5)

            # Add the raw data to self.
            var_idx = 0
            for var in variables:
                if var_index >= 0:
                    if var_idx == var_index:
                        setattr(ext_object, var.strip(),
                                raw_data[:, ...])
                else:
                    setattr(ext_object, var.strip(),
                            raw_data[var_idx,:, ...])
                var_idx += 1

            self.t = t
            self.r_cyl = r_cyl
            self.z_cyl = z_cyl
            self.phiavg = ext_object

        return 0


    def __equal_newline(self, line):
        """
        Determine if string is equal new line.
        """

        return line == '\n'


    def __read_phiaver(self, datadir, variables, aver_file_name,
                       n_vars, var_index, iter_list, l_h5=False):
        """
        Read the PHIAVG file
        Return the time, cylindrical r and z and raw data.
        """

        import os
        import numpy as np
        from scipy.io import FortranFile
        from .. import read

        # Read the data
        if l_h5:
            import h5py
            #
            # Not implemented
            #
        else:
            glob_dim = read.dim(datadir)
            nu=glob_dim.nx/2
            nv=glob_dim.nz

            dim = read.dim(datadir)
            if dim.precision == 'S':
                dtype = np.float32
            if dim.precision == 'D':
                dtype = np.float64

            # Prepare the raw data.
            raw_data = []
            t = []

            # Read records
            #path=os.path.join(datadir, aver_file_name)
            #print(path)
            file_id = FortranFile(os.path.join(datadir, aver_file_name))

            data1 = file_id.read_record(dtype='i4')
            nr_phiavg = data1[0]
            nz_phiavg = data1[1]
            nvars = data1[2]
            nprocz = data1[3]

            data2 = file_id.read_record(dtype=np.float64)
            t = data2[0]
            r_cyl = data2[1:nr_phiavg+1]
            z_cyl = data2[nr_phiavg+1:nr_phiavg+nz_phiavg+1]

            data3 = file_id.read_record(dtype=np.float64)
            raw_data = data3.reshape(nvars,nz_phiavg,nr_phiavg)

            return t, r_cyl, z_cyl, raw_data
