# averages.py
#
# Read the average files.
#
# Author: S. Candelaresi (iomsn1@gmail.com).
"""
Contains the classes and methods to read average files.
"""
import sys
from pencil import read
from pencil.math import natural_sort
import glob

def aver(*args, **kwargs):
    """
    aver(plane_list=None, datadir='data', proc=-1, var_index=-1, proc=-1):

    Read Pencil Code average data.

    Parameters
    ----------
    plane_list : list of string
        A list of the 2d/1d planes over which the averages were taken.
        Takes 'xy', 'xz', 'yz', 'y', 'z'.
        By default, it is [p for p in ["xy", "xz", "yz"] if corresponding_dot_in_file_exists_in_simdir].

    iter_list : list of int
        Iteration indices for which to sample the slices.

    avfile_list , infile_list : list of string
        File names if alternative to standard files used.

    var_index : int
        Index of variable from among within the 'y' or 'z' averages.
        Takes an integer value < len(yaver.in or zaver.in).

    datadir : string
        Directory where the data is stored. By default, "data"

    simdir : string
        Simulation directory containing the .in files.
        By default, current directory.

    proc : int
        Processor to be read. If -1 read all and assemble to one array.
        Only affects the reading of 'yaverages.dat' and 'zaverages.dat'.

    precision : string
        Float (f), double (d) or half (half).

    Returns
    -------
    Class containing the averages.

    Examples
    --------
    >>> aver = pc.read.aver()
    >>> avr.lnrho.shape
    (134, 134, 134)
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
        datadir="data",
        simdir=".",
        plane_list=None,
        avfile_list=None,
        infile_list=None,
        var_index=-1,
        var_names=list(),
        iter_list=None,
        iter_step=1,
        time_range=None,
        param=list(),
        proc=-1,
        precision="f",
    ):
        """
        read(plane_list=None, datadir='data', proc=-1, var_index=-1, proc=-1):

        Read Pencil Code average data.

        Parameters
        ----------
        datadir : string
            Directory where the data is stored. By default, "data"

        simdir : string
            Simulation directory containing the .in files.
            By default, current directory.

        plane_list : list of string
            A list of the 2d/1d planes over which the averages were taken.
            Takes 'xy', 'xz', 'yz', 'y', 'z'.
            By default, it is [p for p in ["xy", "xz", "yz"] if corresponding_dot_in_file_exists_in_simdir].

        avfile_list , infile_list : list of string
            File names if alternative to standard files used.

        var_names : string or list of strings matching format in *aver.in for given plane
            Variable names from among h5 averages.

        var_index : int
            Index of variable from among within the 'y' or 'z' averages.
            Takes an integer value < len(yaver.in or zaver.in).

        time_range : float or list of two floats.
            Single float indicates end time of sample taking 0 to be the start time.
            Two floats indicates start and end time to sample.

        iter_list : list of int, or index range
            Iteration indices for which to sample the slices.
            Single value only that index is sampled. Two indices indicate first and last index.
            More than two indicies specifies the precise list of indices to be sampled.

        iter_step : when sampling iterations number of steps, when iter_list length is two.

        param : Param object
            Param object belonging to the simulation.

        proc : int
            Processor to be read. If -1 read all and assemble to one array.
            Only affects the reading of 'yaverages.dat' and 'zaverages.dat'.

        precision : string
            Float (f), double (d) or half (half).

        Returns
        -------
        Class containing the averages.

        Examples
        --------
        >>> aver = pc.read.aver()
        >>> avr.lnrho.shape
        (134, 134, 134)
        """

        #check iter_list type is integer(s) if present
        if iter_list:
            if isinstance(iter_list, list):
                iter_list = iter_list
            else:
                iter_list = [iter_list]
            for it in iter_list:
                if not isinstance(it, int):
                    print("read.aver Error: iter_list contains {}, but must be integers".format(it))
                    exit()

        import os
        from os.path import join

        lh5 = False
        if isinstance(param, list):
            param = read.param(datadir=datadir)
        if param.io_strategy == "HDF5":
            lh5 = True
        # Keep this for sims that were converted from Fortran to hdf5
        if os.path.exists(os.path.join(datadir, "grid.h5")):
            lh5 = True

        if not lh5:
            #Reading averages from Fortran binary and ascii files
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
                        print(
                            "plane_list and avfile_list must have matching length and order"
                        )
                        sys.stdout.flush()
            else:
                aver_file_name_list = []
                in_file_name_list = []
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
                if plane_list.count("phi") > 0:
                    if os.path.exists(os.path.join(simdir, "phiaver.in")):
                        in_file_name_list.append("data/averages/phiavg.list")
                        aver_file_name_list.append(os.path.join("averages", "phi.h5"))
                if not in_file_name_list:
                    print("error: invalid plane name")
                    sys.stdout.flush()
                    return -1

            class Foo(object):
                pass

            print(plane_list, in_file_name_list, aver_file_name_list)
            for plane, in_file_name, aver_file_name in zip(
                plane_list, in_file_name_list, aver_file_name_list
            ):
                # This one will store the data.
                ext_object = Foo()

                # Get the averaged quantities.
                file_id = open(os.path.join(simdir, in_file_name))
                variables = file_id.readlines()
                file_id.close()
                variables = [
                    v for v in variables if v[0] != "#" and not v.isspace()
                ]  # Ignore commented variables and blank lines in the .in file.
                n_vars = len(variables)

                if plane == "xy" or plane == "xz" or plane == "yz":
                    t, raw_data = self.__read_2d_aver(
                        plane,
                        datadir,
                        variables,
                        aver_file_name,
                        n_vars,
                        precision=precision,
                    )
                if plane == "y" or plane == "z" or plane == "phi":
                    t, raw_data = self.__read_1d_aver(
                        plane,
                        datadir,
                        variables,
                        aver_file_name,
                        n_vars,
                        var_index,
                        iter_list,
                        proc,
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
                plane_keys = ext_object.__dict__.keys()
                setattr(ext_object, 'keys', plane_keys)

                self.t = t
                setattr(self, plane, ext_object)

        else:
            #Reading averages from hdf5 files
            dim = read.dim()
            av_files = glob.glob(join(datadir,"averages","*"))
            if len(av_files) > 0:
                for av_file in av_files:
                    if not ".h5" in av_file[-3:]:
                        av_files.remove(av_file)
                if len(av_files) == 0:
                    print('read.aver error: no averages files in '+join(datadir,"averages"))
                    sys.stdout.flush()
                    return -1
            else:
                print('read.aver error: no averages files in '+join(datadir,"averages"))
                sys.stdout.flush()
                return -1
            av_files_in = list()
            # Initialize the av_files_list of planes.
            if avfile_list:
                if isinstance(avfile_list, list):
                    avfile_list = avfile_list
                else:
                    avfile_list = [avfile_list]
            # Check plane_list if present matches data files and avfile_list if present
            if plane_list:
                if isinstance(plane_list, list):
                    plane_list = plane_list
                else:
                    plane_list = [plane_list]
                for prefix in plane_list:
                    # Check matches in avfile_list if present
                    if avfile_list:
                        if not prefix+".h5" in avfile_list:
                            print(prefix+" does not match in avfile_list and plane_list")
                            plane_list.remove(prefix)
                        else:
                            if join(datadir,"averages",prefix+".h5") in av_files:
                                av_files_in.append(join(datadir,"averages",prefix+".h5"))
                    else:
                        if join(datadir,"averages",prefix+".h5") in av_files:
                            av_files_in.append(join(datadir,"averages",prefix+".h5"))
                        else:
                            plane_list.remove(prefix)
                if len(av_files_in) == 0:
                    print('read.aver error: plane_lists and avlist or av_files have no match'.format(plane_list, av_files))
                    sys.stdout.flush()
                    return -1
            else:
                plane_list = list()
                for av_file in av_files:
                    if avfile_list:
                        if av_file.split('/')[-1] in avfile_list:
                            av_files_in.append(av_file)
                            plane_list.append(av_file.split('/')[-1].split('_')[-1][:-3])
                    else:
                        av_files_in.append(av_file)
                        plane_list.append(av_file.split('/')[-1].split('_')[-1][:-3])
                if len(av_files_in) == 0:
                    print('read.aver error: avfile_list has no match in av_files'.format(avlist, av_files))
                    sys.stdout.flush()
                    return -1

            print(av_files_in)
            for av_file, plane in zip(av_files_in, plane_list):

                # Get the averaged quantities
                t, ext_object = self.__read_h5_aver(
                        plane,
                        av_file,
                        var_names,
                        var_index,
                        iter_list,
                        iter_step,
                        time_range,
                        proc,
                        precision=precision,
                    )

                self.t = t
                setattr(self, plane, ext_object)
        return 0

    def __equal_newline(self, line):
        """
        Determine if string is equal new line.
        """

        return line == "\n"

    def __read_h5_aver(
        self,
        plane,
        av_file,
        var_names,
        var_index,
        iter_list,
        iter_step,
        time_range,
        proc,
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
        import h5py

        print(av_file)
        sys.stdout.flush()
        with h5py.File(av_file, "r") as tmp:
            n_times = len(tmp.keys()) - 1
            if tmp['last'][()].item() < n_times-2:
                n_times = tmp['last'][()].item() + 1
            start_time, end_time = 0, tmp[str(n_times-1)]['time'][()].item()
            if time_range:
                if isinstance(time_range, list):
                    time_range = time_range
                else:
                    time_range = [time_range]
                if len(time_range) == 1:
                    start_time = 0.
                    end_time = time_range[0]
                elif len(time_range) == 2:
                    start_time = time_range[0]
                    end_time = time_range[1]
            if iter_list:
                if len(iter_list) == 1:
                    n_times = min(n_times-1, iter_list[0])
                    itlist = range(n_times,n_times+1,iter_step)
                elif len(iter_list) == 2:
                    if iter_list[0] >= iter_list[1]:
                        print("Warning read.aver: iter_list pair must be list of integer - reading full series instead")
                    itlist = range(max(0,iter_list[0]),min(n_times,iter_list[1]),iter_step)
                else:
                    itlist = iter_list
                for it in itlist:
                    if not str(it) in tmp.keys():
                        itlist.remove(it)
                if not len(itlist) > 0:
                    print("Warning read.aver: iter_list has no match in {{av_file}} keys - reading only {{tmp['last'][0]}} instead")
                    itlist.append(tmp['last'][0])
            else:
                itlist = natural_sort(tmp.keys())[:n_times]

            if time_range:
                tmplist = list()
                for t_idx in itlist:
                    if tmp[str(t_idx) + "/time"][()].item() >= start_time:
                        if tmp[str(t_idx) + "/time"][()].item() <= end_time:
                            tmplist.append(t_idx)
                if len(tmplist) == 0:
                    print('read.aver error: no data in {{av_file}} within time range {{time_range}}.')
                    sys.stdout.flush()
                    return -1
                else:
                    itlist = tmplist
            # Determine the structure of the xy/xz/yz/y/z averages.
            if len(var_names) > 0:
                if isinstance(var_names, list):
                    var_names = var_names
                else:
                    var_names = [var_names]
                for var_name in var_names:
                    if not var_name in tmp[str(itlist[0])].keys():
                        var_names.remove(var_name)
                if len(var_names) < 1:
                    print("Warning read.aver: var_names has no match in {{av_file}} keys - reading all variables instead")
                    var_names = list(tmp[str(itlist[0])].keys())
            else:
                var_names = list(tmp[str(itlist[0])].keys())
            if "time" in var_names:
                var_names.remove('time')

            class Foo(object):
                pass

            # This one will store the data.
            if not hasattr(self, plane):
                ext_object = Foo()
            else:
                ext_object = self.__getattribute__(plane)

            for var in var_names:
                data_shape = list()
                data_shape.append(len(itlist))
                dshape = tmp[str(itlist[0])+ "/" + var.strip()].shape
                for dsh in dshape:
                    data_shape.append(dsh)
                raw_data = np.zeros(data_shape, dtype=precision)
                t = np.zeros(data_shape[0], dtype=precision)
                for t_idx, tmp_idx in zip(range(len(itlist)),itlist):
                    t[int(t_idx)] = tmp[str(tmp_idx) + "/time"][()]
                    if var.strip() in tmp[str(tmp_idx)].keys():
                        raw_data[int(t_idx)] = tmp[str(tmp_idx) + "/" + var.strip()][()]

                setattr(ext_object, var.strip(), raw_data)
            plane_keys = list(ext_object.__dict__.keys())
            if 'keys' in plane_keys:
                plane_keys.remove('keys')
            setattr(ext_object, 'keys', plane_keys)

        return t, ext_object

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
        precision="f",
    ):
        """
        Read the xyaverages.dat, xzaverages.dat, yzaverages.dat
        Return the raw data and the time array.
        """

        import os
        import numpy as np
        from pencil import read

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
        raw_data = np.zeros([n_times, n_vars * nw], dtype=precision)
        line_idx = 0
        t_idx = -1
        try:
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
        except:
            print("Error: There was a problem reading {} at line {}.\nCalculated values: n_vars = {}, nw = {}.\nAre these correct?".format(aver_file_name, line_idx, n_vars, nw))
            raise

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
