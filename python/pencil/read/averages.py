# averages.py
#
# Read the average files.
#
# Author: S. Candelaresi (iomsn1@gmail.com).
"""
Contains the classes and methods to read average files.
"""
import sys


def aver(*args, **kwargs):
    """
    Read Pencil Code average data.

    call signature:

    read(plane_list=['xy', 'xz', 'yz'], datadir='data', proc=-1, var_index=-1, proc=-1):

    Keyword arguments:

    *plane_list*:
        A list of the 2d/1d planes over which the averages were taken.
        Takes 'xy', 'xz', 'yz', 'y', 'z'.
        By default, it is [p for p in ["xy", "xz", "yz"] if corresponding_dot_in_file_exists_in_simdir]

    *iter_list*
        list of iteration indices for which to sample the slices

    *avfile_list*, *infile_list*
        list of file names if alternative to standard files used

    *var_index*:
        Index of variable from among within the 'y' or 'z' averages.
        Takes an integer value < len(yaver.in or zaver.in).

    *datadir*:
        Directory where the data is stored. By default, "./data"

    *simdir*:
        Simulation directory containing the .in files.
        By default, parent directory of datadir.

    *proc*:
        Processor to be read. If -1 read all and assemble to one array.
        Only affects the reading of 'yaverages.dat' and 'zaverages.dat'.

    *precision*
        Float (f), double (d) or half (half).

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

    def keys(self):
        for i in self.__dict__.keys():
            print(i)

    def read(
        self,
        plane_list=None,
        var_index=-1,
        datadir="data",
        simdir=None,
        proc=-1,
        iter_list=None,
        avfile_list=None,
        infile_list=None,
        precision="f",
    ):
        """
        Read Pencil Code average data.

        call signature:

        read(plane_list=['xy', 'xz', 'yz'], datadir='data', proc=-1, var_index=-1, proc=-1):

        Keyword arguments:

        *plane_list*:
          A list of the 2d/1d planes over which the averages were taken.
          Takes 'xy', 'xz', 'yz', 'y', 'z'.
          By default, it is [p for p in ["xy", "xz", "yz"] if corresponding_dot_in_file_exists_in_simdir]

        *iter_list*
          list of iteration indices for which to sample the slices

        *avfile_list*, *infile_list*
          list of file names if alternative to standard files used

        *var_index*:
          Index of variable from among within the 'y' or 'z' averages.
          Takes an integer value < len(yaver.in or zaver.in).

        *datadir*:
          Directory where the data is stored.

        *simdir*:
          Simulation directory containing the .in files.
          By default, parent directory of datadir.

        *proc*:
          Processor to be read. If -1 read all and assemble to one array.
          Only affects the reading of 'yaverages.dat' and 'zaverages.dat'.

        *precision*
          Float (f), double (d) or half (half).

        """

        import os

        l_h5 = False

        if simdir is None:
            simdir = os.path.join(datadir, os.path.pardir)

        # Initialize the planes list.
        if plane_list:
            if isinstance(plane_list, list):
                plane_list = plane_list
            else:
                plane_list = [plane_list]
        else:
            plane_list = []
            for prefix in ["xy", "xz", "yz"]:
                if os.path.exists(os.path.join(simdir, prefix + "aver.in")):
                    plane_list.append(prefix)

        if var_index >= 0:
            print("var_index {} requires plane_list = 'y' or 'z',".format(var_index))
            sys.stdout.flush()
        # Determine which average files to read.
        if os.path.exists(os.path.join(datadir, "grid.h5")):
            l_h5 = True
        if avfile_list:
            if isinstance(avfile_list, list):
                aver_file_name_list = avfile_list
            else:
                aver_file_name_list = [avfile_list]
            if infile_list:
                if isinstance(infile_list, list):
                    in_file_name_list = infile_list
                else:
                    in_file_name_list = [infile_list]
            else:
                in_file_name_list = []
                if len(aver_file_name_list) == len(plane_list):
                    for plane in plane_list:
                        in_file_name_list.append(prefix + "aver.in")
                else:
                    print("plane_list and avfile_list must have matching length and order")
                    sys.stdout.flush()
        else:
            aver_file_name_list = []
            in_file_name_list = []
            if os.path.exists(os.path.join(datadir, "grid.h5")):
                l_h5 = True
                if plane_list.count("xy") > 0:
                    if os.path.exists(os.path.join(simdir, "xyaver.in")):
                        in_file_name_list.append("xyaver.in")
                        aver_file_name_list.append(os.path.join("averages", "xy.h5"))
                if plane_list.count("xz") > 0:
                    if os.path.exists(os.path.join(simdir, "xzaver.in")):
                        in_file_name_list.append("xzaver.in")
                        aver_file_name_list.append(os.path.join("averages", "xz.h5"))
                if plane_list.count("yz") > 0:
                    if os.path.exists(os.path.join(simdir, "yzaver.in")):
                        in_file_name_list.append("yzaver.in")
                        aver_file_name_list.append(os.path.join("averages", "yz.h5"))
                if plane_list.count("y") > 0:
                    if os.path.exists(os.path.join(simdir, "yaver.in")):
                        in_file_name_list.append("yaver.in")
                        aver_file_name_list.append(os.path.join("averages", "y.h5"))
                if plane_list.count("z") > 0:
                    if os.path.exists(os.path.join(simdir, "zaver.in")):
                        in_file_name_list.append("zaver.in")
                        aver_file_name_list.append(os.path.join("averages", "z.h5"))
            else:
                if plane_list.count("xy") > 0:
                    if os.path.exists(os.path.join(simdir, "xyaver.in")):
                        in_file_name_list.append("xyaver.in")
                        aver_file_name_list.append("xyaverages.dat")
                if plane_list.count("xz") > 0:
                    if os.path.exists(os.path.join(simdir, "xzaver.in")):
                        in_file_name_list.append("xzaver.in")
                        aver_file_name_list.append("xzaverages.dat")
                if plane_list.count("yz") > 0:
                    if os.path.exists(os.path.join(simdir, "yzaver.in")):
                        in_file_name_list.append("yzaver.in")
                        aver_file_name_list.append("yzaverages.dat")
                if plane_list.count("y") > 0:
                    if os.path.exists(os.path.join(simdir, "yaver.in")):
                        in_file_name_list.append("yaver.in")
                        aver_file_name_list.append("yaverages.dat")
                if plane_list.count("z") > 0:
                    if os.path.exists(os.path.join(simdir, "zaver.in")):
                        in_file_name_list.append("zaver.in")
                        aver_file_name_list.append("zaverages.dat")
            if not in_file_name_list:
                print("error: invalid plane name")
                sys.stdout.flush()
                return -1

        class Foo(object):
            pass

        print(plane_list,in_file_name_list, aver_file_name_list)
        for plane, in_file_name, aver_file_name in zip(
            plane_list, in_file_name_list, aver_file_name_list
        ):
            # This one will store the data.
            ext_object = Foo()

            # Get the averaged quantities.
            file_id = open(os.path.join(simdir, in_file_name))
            variables = file_id.readlines()
            file_id.close()
            variables = [v for v in variables if v[0] != '#' and not v.isspace()] #Ignore commented variables and blank lines in the .in file.
            n_vars = len(variables)

            if plane == "xy" or plane == "xz" or plane == "yz":
                t, raw_data = self.__read_2d_aver(
                    plane,
                    datadir,
                    variables,
                    aver_file_name,
                    n_vars,
                    l_h5=l_h5,
                    precision=precision,
                )
            if plane == "y" or plane == "z":
                t, raw_data = self.__read_1d_aver(
                    plane,
                    datadir,
                    variables,
                    aver_file_name,
                    n_vars,
                    var_index,
                    iter_list,
                    proc,
                    l_h5=l_h5,
                    precision=precision,
                )

            # Add the raw data to self.
            var_idx = 0
            for var in variables:
                if var_index >= 0:
                    if var_idx == var_index:
                        setattr(ext_object, var.strip(), raw_data[:, ...])
                else:
                    setattr(ext_object, var.strip(), raw_data[:, var_idx, ...])
                var_idx += 1

            self.t = t
            setattr(self, plane, ext_object)

        return 0

    def __equal_newline(self, line):
        """
        Determine if string is equal new line.
        """

        return line == "\n"

    def __read_1d_aver(
        self,
        plane,
        datadir,
        variables,
        aver_file_name,
        n_vars,
        var_index,
        iter_list,
        proc,
        l_h5=False,
        precision="f",
    ):
        """
        Read the yaverages.dat, zaverages.dat.
        Return the raw data and the time array.
        """

        import os
        import numpy as np
        from scipy.io import FortranFile
        from pencil import read

        # Read the data
        if l_h5:
            import h5py

            file_id = os.path.join(datadir, aver_file_name)
            print(file_id)
            sys.stdout.flush()
            with h5py.File(file_id, "r") as tmp:
                n_times = len(tmp.keys()) - 1
                # Determine the structure of the xy/xz/yz averages.
                for var in variables:
                    nu = tmp[str(0) + "/" + var.strip()].shape[0]
                    nv = tmp[str(0) + "/" + var.strip()].shape[1]
                    break
            raw_data = np.zeros([n_times, n_vars, nu, nv], dtype=precision)
            t = np.zeros(n_times, dtype=precision)
            with h5py.File(file_id, "r") as tmp:
                for t_idx in range(0, n_times):
                    t[t_idx] = tmp[str(t_idx) + "/time"][()]
                    raw_idx = 0
                    for var in variables:
                        raw_data[t_idx, raw_idx] = tmp[str(t_idx) + "/" + var.strip()][
                            ()
                        ]
                        raw_idx += 1
        else:
            glob_dim = read.dim(datadir)
            if plane == "y":
                nu = glob_dim.nx
                nv = glob_dim.nz
            if plane == "z":
                nu = glob_dim.nx
                nv = glob_dim.ny

            if proc < 0:
                offset = glob_dim.nprocx * glob_dim.nprocy
                if plane == "z":
                    proc_list = range(offset)
                if plane == "y":
                    proc_list = []
                    xr = range(glob_dim.nprocx)
                    for iz in range(glob_dim.nprocz):
                        proc_list.extend(xr)
                        xr = [x + offset for x in xr]
                all_procs = True
            else:
                proc_list = [proc]
                all_procs = False

            dim = read.dim(datadir, proc)
            if dim.precision == "S":
                read_precision = np.float32
            if dim.precision == "D":
                read_precision = np.float64

            # Prepare the raw data.
            # This will be reformatted at the end.
            raw_data = []
            for proc in proc_list:
                proc_dir = "proc{0}".format(proc)
                proc_dim = read.dim(datadir, proc)
                if plane == "y":
                    pnu = proc_dim.nx
                    pnv = proc_dim.nz
                if plane == "z":
                    pnu = proc_dim.nx
                    pnv = proc_dim.ny
                if var_index >= 0:
                    inx1 = var_index * pnu * pnv
                    inx2 = (var_index + 1) * pnu * pnv
                # Read the data.
                t = []
                proc_data = []
                try:
                    file_id = FortranFile(
                        os.path.join(datadir, proc_dir, aver_file_name)
                    )
                except:
                    # Not all proc dirs have a [yz]averages.dat.
                    print("Averages of processor {0} missing.".format(proc))
                    sys.stdout.flush()
                    break
                if iter_list:
                    if isinstance(iter_list, list):
                        iter_list = iter_list
                    else:
                        iter_list = [iter_list]
                    # split by iteration overrules split by variable
                    var_index = -1
                    iiter = 0
                    while True:
                        try:
                            if iiter in iter_list:
                                t.append(file_id.read_record(dtype=read_precision)[0])
                                proc_data.append(
                                    file_id.read_record(dtype=read_precision)
                                )
                                if iiter >= iter_list[-1]:
                                    # Finished reading.
                                    break
                                iiter += 1
                            else:
                                file_id.read_record(dtype=read_precision)[0]
                                file_id.read_record(dtype=read_precision)
                                iiter += 1
                        except:
                            # Finished reading.
                            break
                else:
                    while True:
                        try:
                            t.append(file_id.read_record(dtype=read_precision)[0])
                            if var_index >= 0:
                                proc_data.append(
                                    file_id.read_record(dtype=read_precision)[
                                        inx1:inx2
                                    ].astype(precision)
                                )
                            else:
                                proc_data.append(
                                    file_id.read_record(dtype=read_precision).astype(
                                        precision
                                    )
                                )
                        except:
                            # Finished reading.
                            break
                file_id.close()
                # Reshape the proc data into [len(t), pnu, pnv].
                proc_data = np.array(proc_data, dtype=precision)
                if var_index >= 0:
                    proc_data = proc_data.reshape([len(t), 1, pnv, pnu])
                else:
                    proc_data = proc_data.reshape([len(t), n_vars, pnv, pnu])

                if not all_procs:
                    return np.array(t, dtype=precision), proc_data.swapaxes(2, 3)

                # Add the proc_data (one proc) to the raw_data (all procs)
                if plane == "y":
                    if all_procs:
                        idx_u = proc_dim.ipx * proc_dim.nx
                        idx_v = proc_dim.ipz * proc_dim.nz
                    else:
                        idx_v = 0
                        idx_u = 0
                if plane == "z":
                    if all_procs:
                        idx_u = proc_dim.ipx * proc_dim.nx
                        idx_v = proc_dim.ipy * proc_dim.ny
                    else:
                        idx_v = 0
                        idx_u = 0

                if not isinstance(raw_data, np.ndarray):
                    # Initialize the raw_data array with correct dimensions.
                    if var_index >= 0:
                        raw_data = np.zeros([len(t), 1, nv, nu], dtype=precision)
                    else:
                        raw_data = np.zeros([len(t), n_vars, nv, nu], dtype=precision)
                raw_data[
                    :, :, idx_v : idx_v + pnv, idx_u : idx_u + pnu
                ] = proc_data.copy()

            t = np.array(t, dtype=precision)
            raw_data = np.swapaxes(raw_data, 2, 3)

        return t, raw_data

    def __read_2d_aver(
        self,
        plane,
        datadir,
        variables,
        aver_file_name,
        n_vars,
        l_h5=False,
        precision="f",
    ):
        """
        Read the xyaverages.dat, xzaverages.dat, yzaverages.dat
        Return the raw data and the time array.
        """

        import os
        import numpy as np
        from pencil import read

        if l_h5:
            import h5py

            file_id = os.path.join(datadir, aver_file_name)
            print(file_id)
            sys.stdout.flush()
            with h5py.File(file_id, "r") as tmp:
                n_times = len(tmp.keys()) - 1
                # Determine the structure of the xy/xz/yz averages.
                for var in variables:
                    nw = tmp[str(0) + "/" + var.strip()].shape[0]
                    break
        else:
            # Determine the structure of the xy/xz/yz averages.
            if plane == "xy":
                nw = getattr(read.dim(datadir=datadir), "nz")
            if plane == "xz":
                nw = getattr(read.dim(datadir=datadir), "ny")
            if plane == "yz":
                nw = getattr(read.dim(datadir=datadir), "nx")
            file_id = open(os.path.join(datadir, aver_file_name))
            aver_lines = file_id.readlines()
            file_id.close()
            entry_length = int(np.ceil(nw * n_vars / 8.0))
            n_times = int(len(aver_lines) / (1.0 + entry_length))

        # Prepare the data arrays.
        t = np.zeros(n_times, dtype=precision)

        # Read the data
        if l_h5:
            raw_data = np.zeros([n_times, n_vars, nw], dtype=precision)
            with h5py.File(file_id, "r") as tmp:
                for t_idx in range(0, n_times):
                    t[t_idx] = tmp[str(t_idx) + "/time"][()]
                    raw_idx = 0
                    for var in variables:
                        raw_data[t_idx, raw_idx] = tmp[str(t_idx) + "/" + var.strip()][
                            ()
                        ]
                        raw_idx += 1
        else:
            raw_data = np.zeros([n_times, n_vars * nw], dtype=precision)
            line_idx = 0
            t_idx = -1
            for current_line in aver_lines:
                if line_idx % (entry_length + 1) == 0:
                    t_idx += 1
                    t[t_idx] = current_line
                    raw_idx = 0
                else:
                    raw_data[t_idx, raw_idx * 8 : (raw_idx * 8 + 8)] = list(
                        map(np.float32, current_line.split())
                    )
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
        alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", key)]
        return sorted(l, key=alphanum_key)
