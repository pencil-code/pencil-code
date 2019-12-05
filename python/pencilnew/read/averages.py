# averages.py
#
# Read the average files.
#
# Author: S. Candelaresi (iomsn1@gmail.com).
"""
Contains the classes and methods to read average files.
"""


def aver(*args, **kwargs):
    """
    Read Pencil Code average data.

    call signature:

    read(plane_list=['xy', 'xz', 'yz'], datadir='data', proc=-1):

    Keyword arguments:

    *plane_list*:
      A list of the 2d/1d planes over which the averages were taken.
      Takes 'xy', 'xz', 'yz', 'y', 'z'.

   *var_index*:
     Index of single variable taken from among the 'y' or 'z' averages.
     Takes an integer value < len(yaver.in or zaver.in).

    *datadir*:
      Directory where the data is stored.

    *proc*:
      Processor to be read. If -1 read all and assemble to one array.
      Only affects the reading of 'yaverages.dat' and 'zaverages.dat'.

    *millennium_bug*
      Correction required for missing data on proc 0 in millennium test
      field run - temp fix
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


    def read(self, plane_list=None, var_index=-1, datadir='data',
             proc=-1, millennium_bug=False):
        """
        Read Pencil Code average data.

        call signature:

        read(plane_list=['xy', 'xz', 'yz'], datadir='data', proc=-1):

        Keyword arguments:

        *plane_list*:
          A list of the 2d/1d planes over which the averages were taken.
          Takes 'xy', 'xz', 'yz', 'y', 'z'.

        *var_index*:
             Index of variable from among within the 'y' or 'z' averages.
             Takes an integer value < len(yaver.in or zaver.in).

        *datadir*:
          Directory where the data is stored.

        *proc*:
          Processor to be read. If -1 read all and assemble to one array.
          Only affects the reading of 'yaverages.dat' and 'zaverages.dat'.

        *millennium_bug*
          Correction required for missing data on proc 0 in millennium test
          field run - temp fix
        """

        import os

        l_h5 = False

        # Initialize the planes list.
        if plane_list:
            if isinstance(plane_list, list):
                plane_list = plane_list
            else:
                plane_list = [plane_list]
        else:
            plane_list = ['xy', 'xz', 'yz']

        if var_index >= 0:
            print("var_index {} requires plane_list = 'y' or 'z',".format(var_index))
        # Determine which average files to read.
        in_file_name_list = []
        aver_file_name_list = []
        if os.path.exists(os.path.join(datadir, 'grid.h5')):
            l_h5 = True
            if plane_list.count('xy') > 0:
                in_file_name_list.append('xyaver.in')
                aver_file_name_list.append(os.path.join('averages', 'xy.h5'))
            if plane_list.count('xz') > 0:
                in_file_name_list.append('xzaver.in')
                aver_file_name_list.append(os.path.join('averages', 'xz.h5'))
            if plane_list.count('yz') > 0:
                in_file_name_list.append('yzaver.in')
                aver_file_name_list.append(os.path.join('averages', 'yz.h5'))
            if plane_list.count('y') > 0:
                in_file_name_list.append('yaver.in')
                aver_file_name_list.append(os.path.join('averages', 'y.h5'))
            if plane_list.count('z') > 0:
                in_file_name_list.append('zaver.in')
                aver_file_name_list.append(os.path.join('averages', 'z.h5'))
        else:
            if plane_list.count('xy') > 0:
                in_file_name_list.append('xyaver.in')
                aver_file_name_list.append('xyaverages.dat')
            if plane_list.count('xz') > 0:
                in_file_name_list.append('xzaver.in')
                aver_file_name_list.append('xzaverages.dat')
            if plane_list.count('yz') > 0:
                in_file_name_list.append('yzaver.in')
                aver_file_name_list.append('yzaverages.dat')
            if plane_list.count('y') > 0:
                in_file_name_list.append('yaver.in')
                aver_file_name_list.append('yaverages.dat')
            if plane_list.count('z') > 0:
                in_file_name_list.append('zaver.in')
                aver_file_name_list.append('zaverages.dat')
        if not in_file_name_list:
            print("error: invalid plane name")
            return -1

        class Foo(object):
            pass

        for plane, in_file_name, aver_file_name in \
        zip(plane_list, in_file_name_list, aver_file_name_list):
            # This one will store the data.
            ext_object = Foo()

            # Get the averaged quantities.
            file_id = open(os.path.join(os.path.dirname(datadir), in_file_name))
            variables = file_id.readlines()
            file_id.close()
            for i in range(sum(list(map(self.__equal_newline, variables)))):
                variables.remove('\n')
            n_vars = len(variables)

            if plane == 'xy' or plane == 'xz' or plane == 'yz':
                t, raw_data = self.__read_2d_aver(plane, datadir, variables,
                                                  aver_file_name, n_vars, l_h5)
            if plane == 'y' or plane == 'z':
                t, raw_data = self.__read_1d_aver(plane, datadir, variables,
                                                  aver_file_name, n_vars,
                                                  var_index, proc, l_h5,
                                                  millennium_bug)

            # Add the raw data to self.
            var_idx = 0
            for var in variables:
                if var_index >= 0:
                    if var_idx == var_index:
                        setattr(ext_object, var.strip(),
                                raw_data[:, ...])
                else:
                    setattr(ext_object, var.strip(),
                            raw_data[:, var_idx, ...])
                var_idx += 1

            self.t = t
            setattr(self, plane, ext_object)

        return 0


    def __equal_newline(self, line):
        """
        Determine if string is equal new line.
        """

        return line == '\n'


    def __read_1d_aver(self, plane, datadir, variables, aver_file_name,
                       n_vars, var_index, proc, l_h5, millennium_bug):
        """
        Read the yaverages.dat, zaverages.dat.
        Return the raw data and the time array.
        """

        import os
        import numpy as np
        from scipy.io import FortranFile
        from .. import read

        # Read the data
        if l_h5:
            import h5py
            file_id = os.path.join(datadir, aver_file_name)
            print(file_id)
            with h5py.File(file_id, 'r') as tmp:
                n_times = len(tmp.keys()) - 1
                # Determine the structure of the xy/xz/yz averages.
                for var in variables:
                    nu = tmp[str(0) + '/' + var.strip()].shape[0]
                    nv = tmp[str(0) + '/' + var.strip()].shape[1]
                    break
            raw_data = np.zeros([n_times, n_vars, nv, nu])
            t = np.zeros(n_times, dtype=np.float32)
            with h5py.File(file_id, 'r') as tmp:
                for t_idx in range(0, n_times):
                    t[t_idx] = tmp[str(t_idx) + '/time'][()]
                    raw_idx = 0
                    for var in variables:
                        raw_data[t_idx, raw_idx] = \
                                    tmp[str(t_idx) + '/' +var.strip()][()]
                        raw_idx += 1
        else:
            glob_dim = read.dim(datadir)
            if plane == 'y':
                nu = glob_dim.nx
                nv = glob_dim.nz
            if plane == 'z':
                nu = glob_dim.nx
                nv = glob_dim.ny

            if proc < 0:
                offset = glob_dim.nprocx*glob_dim.nprocy
                if plane == 'z':
                    proc_list = range(offset)
                if plane == 'y':
                    proc_list = []
                    xr = range(glob_dim.nprocx)
                    for iz in range(glob_dim.nprocz):
                        proc_list.extend(xr)
                        xr = [x+offset for x in xr]
                all_procs = True
            else:
                proc_list = [proc]
                all_procs = False

            dim = read.dim(datadir, proc)
            if dim.precision == 'S':
                dtype = np.float32
            if dim.precision == 'D':
                dtype = np.float64

            # Prepare the raw data.
            # This will be reformatted at the end.
            raw_data = []
            for proc in proc_list:
                proc_dir = 'proc{0}'.format(proc)
                proc_dim = read.dim(datadir, proc)
                if plane == 'y':
                    pnu = proc_dim.nx
                    pnv = proc_dim.nz
                if plane == 'z':
                    pnu = proc_dim.nx
                    pnv = proc_dim.ny
                if var_index >= 0:
                    inx1 = var_index*pnu*pnv
                    inx2 = (var_index+1)*pnu*pnv
                # Read the data.
                t = []
                proc_data = []
                try:
                    file_id = FortranFile(os.path.join(datadir, proc_dir, aver_file_name))
                except:
                    # Not all proc dirs have a [yz]averages.dat.
                    print("Averages of processor {0} missing.".format(proc))
                    break
                while True:
                    try:
                        t.append(file_id.read_record(dtype=dtype)[0])
                        if var_index >= 0:
                            proc_data.append(file_id.read_record(dtype=dtype)[inx1:inx2])
                        else:
                            proc_data.append(file_id.read_record(dtype=dtype))
                    except:
                        # Finished reading.
                        break
                file_id.close()
                # Reshape the proc data into [len(t), pnu, pnv].
                proc_data = np.array(proc_data)
                if var_index >= 0:
                    proc_data = proc_data.reshape([len(t), 1, pnv, pnu])
                    if millennium_bug and proc==1:
                        proc_data=np.insert(proc_data,3766,0.5*(
                                            proc_data[3766]+
                                            proc_data[3767]),axis=0)
                else:
                    proc_data = proc_data.reshape([len(t), n_vars, pnv, pnu])

                if not all_procs:
                    return np.array(t), proc_data.swapaxes(2, 3)

                # Add the proc_data (one proc) to the raw_data (all procs)
                if plane == 'y':
                    if all_procs:
                        idx_u = proc_dim.ipx*proc_dim.nx
                        idx_v = proc_dim.ipz*proc_dim.nz
                    else:
                        idx_v = 0
                        idx_u = 0
                if plane == 'z':
                    if all_procs:
                        idx_u = proc_dim.ipx*proc_dim.nx
                        idx_v = proc_dim.ipy*proc_dim.ny
                    else:
                        idx_v = 0
                        idx_u = 0

                if not isinstance(raw_data, np.ndarray):
                    #Initialize the raw_data array with correct dimensions.
                    if var_index >= 0:
                        raw_data = np.zeros([len(t), 1, nv, nu])
                    else:
                        raw_data = np.zeros([len(t), n_vars, nv, nu])
                raw_data[:, :, idx_v:idx_v+pnv, idx_u:idx_u+pnu] = proc_data.copy()

            t = np.array(t)
            raw_data = np.swapaxes(raw_data, 2, 3)

        return t, raw_data


    def __read_2d_aver(self, plane, datadir, variables, aver_file_name, n_vars, l_h5):
        """
        Read the xyaverages.dat, xzaverages.dat, yzaverages.dat
        Return the raw data and the time array.
        """

        import os
        import numpy as np
        from .. import read

        if l_h5:
            import h5py
            file_id = os.path.join(datadir, aver_file_name)
            print(file_id)
            with h5py.File(file_id, 'r') as tmp:
                n_times = len(tmp.keys()) - 1
                # Determine the structure of the xy/xz/yz averages.
                for var in variables:
                    nw = tmp[str(0) + '/' + var.strip()].shape[0]
                    break
        else:
            # Determine the structure of the xy/xz/yz averages.
            if plane == 'xy':
                nw = getattr(read.dim(), 'nz')
            if plane == 'xz':
                nw = getattr(read.dim(), 'ny')
            if plane == 'yz':
                nw = getattr(read.dim(), 'nx')
            file_id = open(os.path.join(datadir, aver_file_name))
            aver_lines = file_id.readlines()
            file_id.close()
            entry_length = int(np.ceil(nw*n_vars/8.))
            n_times = int(len(aver_lines)/(1. + entry_length))

        # Prepare the data arrays.
        t = np.zeros(n_times, dtype=np.float32)

        # Read the data
        if l_h5:
            raw_data = np.zeros([n_times, n_vars, nw])
            with h5py.File(file_id, 'r') as tmp:
                for t_idx in range(0, n_times):
                    t[t_idx] = tmp[str(t_idx) + '/time'][()]
                    raw_idx = 0
                    for var in variables:
                        raw_data[t_idx, raw_idx] = tmp[str(t_idx) + '/' + var.strip()][()]
                        raw_idx += 1
        else:
            raw_data = np.zeros([n_times, n_vars*nw])
            line_idx = 0
            t_idx = -1
            for current_line in aver_lines:
                if line_idx % (entry_length+1) == 0:
                    t_idx += 1
                    t[t_idx] = np.float32(current_line)
                    raw_idx = 0
                else:
                    raw_data[t_idx, raw_idx*8:(raw_idx*8+8)] = \
                        list(map(np.float32, current_line.split()))
                    raw_idx += 1
                line_idx += 1

            # Restructure the raw data and add it to the Averages object.
            raw_data = np.reshape(raw_data, [n_times, n_vars, nw])

        return t, raw_data


    def __natural_sort(self, l):
        """
        Sort array in a more natural way, e.g. 9VAR < 10VAR
        """

        import re

        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(l, key=alphanum_key)
