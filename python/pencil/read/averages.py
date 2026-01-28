# averages.py
#
# Read the average files.
#
# Author: S. Candelaresi (iomsn1@gmail.com).
"""
Contains the classes and methods to read average files.
"""
import warnings
import sys
import numpy as np
from pencil import read
from pencil.math import natural_sort
from pencil.util import copy_docstring
import glob
import time
import os
from collections.abc import Iterable

class Averages(object):
    """
    Averages -- holds Pencil Code averages data and methods.
    """

    def __init__(self):
        """
        Fill members with default values.
        """
        self._t = np.array([])

    @property
    def t(self):
        return self._t

    @t.setter
    def t(self, arr):
        if len(self.t) > 0 and (len(self.t) != len(arr) or np.any(self.t != arr)):
            warnings.warn("Mismatch between the times of different kinds of averages (usually happens when 1D and 2D averages are stored at different times). Please use the t attributes of the respective planes (e.g. av.xy.t, rather than av.t).")
        self._t = arr

    def keys(self):
        for i in self.__dict__.keys():
            if i == "_t":
               i = "t"
            print(i)

    def read(
        self,
        datadir="data",
        simdir=".",
        plane_list=None,
        avfile_list=None,
        infile_list=None,
        var_index=-1,
        var_names=None,
        iter_list=None,
        niter = 9999,
        iter_step=1,
        time_range=None,
        param=None,
        proc=-1,
        precision="f",
        comp_time=False,
        quiet=True
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
            Single value only that index is sampled.
            Two indices indicate first and last index (last value excluded).
            More than two indicies specifies the precise list of indices to be sampled.
            Note that the indices specified here start from 0.

        niter : int
            if iter_list is not used for fortran format a default index list will be used of size niter.

        iter_step : when sampling iterations number of steps, when iter_list length is two.

        param : Param object
            Param object belonging to the simulation.

        proc : int
            Processor to be read. If -1 read all and assemble to one array.
            Only affects the reading of 'yaverages.dat' and 'zaverages.dat'.

        precision : string
            Float (f), double (d) or half (half).

        quiet : bool
            Whether to suppress diagnostic output.
            Default: True

        Returns
        -------
        Class containing the averages.

        Examples
        --------
        >>> aver = pc.read.aver()
        >>> aver.lnrho.shape
        (134, 134, 134)

        Notes
        -----
        The axis ordering of the 2D averages (yaver and zaver) is [t,x,z] or [t,x,y]
        """


        #check iter_list type is integer(s) if present
        if iter_list is not None:
            if isinstance(iter_list, Iterable):
                iter_list = list(iter_list)
            else:
                iter_list = [iter_list]
            for it in iter_list:
                if not isinstance(it, int):
                    raise ValueError(f"read.aver Error: iter_list contains {it}, but must be integers")

        if time_range is not None:
            if isinstance(time_range, Iterable):
                time_range = list(time_range)
            else:
                time_range = [0, time_range]

        from os.path import join, abspath

        simdir = abspath(simdir)

        if param is None:
            param = read.param(datadir=datadir, quiet=True)

        if param.io_strategy != "HDF5":
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

            #if var_index >= 0:
            #    print("var_index {} requires plane_list = 'y' or 'z',".format(var_index))
            #    sys.stdout.flush()
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
                        raise ValueError(
                            "plane_list and avfile_list must have matching length and order"
                        )
            else:
                good_planes = []
                aver_file_name_list = []
                in_file_name_list = []
                for plane in plane_list:
                    if plane in ["xy", "xz", "yz", "y", "z"] and os.path.exists(os.path.join(simdir, f"{plane}aver.in")):
                        good_planes.append(plane)
                        in_file_name_list.append(f"{plane}aver.in")
                        aver_file_name_list.append(f"{plane}averages.dat")
                    elif plane == "phi" and os.path.exists(os.path.join(simdir, "phiaver.in")):
                        # Kishore (2025-Jun-09) is the use of a HDF5 file below correct? TODO
                        good_planes.append("phi")
                        in_file_name_list.append("data/averages/phiavg.list")
                        aver_file_name_list.append(os.path.join("averages", "phi.h5"))
                plane_list = good_planes

                if len(good_planes) == 0:
                    raise ValueError("All the specified planes are invalid")

            if not quiet:
                print(plane_list, in_file_name_list, aver_file_name_list)

            for plane, in_file_name, aver_file_name in zip(
                plane_list, in_file_name_list, aver_file_name_list
            ):
                # This one will store the data.
                ext_object = _Plane()

                # Get the averaged quantities.
                file_id = open(os.path.join(simdir, in_file_name))
                variables = file_id.readlines()
                file_id.close()
                variables = [
                    v.strip("\n") for v in variables if v[0] != "#" and not v.isspace()
                ]  # Ignore commented variables and blank lines in the .in file.
                n_vars = len(variables)
                if not quiet:
                    print(variables)
                if var_names is not None:
                    if isinstance(var_names, str):
                        plane_var_names = [var_names]
                    else:
                        plane_var_names = list(var_names)
                    for var_name in plane_var_names:
                        if not var_name in variables:
                            plane_var_names.remove(var_name)
                    if len(plane_var_names) < 1:
                        warnings.warn("read.aver: var_names has no match in {} - reading all variables instead".format(in_file_name))
                    else:
                        var_index = list()
                        for indx, var in zip(range(n_vars),variables):
                            if var in plane_var_names:
                                var_index.append(indx)
                if plane == "xy" or plane == "xz" or plane == "yz":
                    t, raw_data = self._read_2d_aver(
                        plane,
                        datadir,
                        aver_file_name,
                        n_vars,
                        var_names,
                        var_index,
                        iter_list,
                        iter_step,
                        time_range,
                        proc,
                        precision=precision,
                    )
                elif plane == "y" or plane == "z" or plane == "phi":
                    t, raw_data = self._read_1d_aver(
                        plane,
                        datadir,
                        aver_file_name,
                        n_vars,
                        var_names,
                        var_index,
                        iter_list,
                        iter_step,
                        time_range,
                        proc,
                        precision=precision,
                    )
                else:
                    raise ValueError(f"Unknown plane {plane}.")

                # Add the raw data to self.
                var_idx = 0
                if not isinstance(var_index,list):
                    for var in variables:
                        if not quiet:
                            print('var_index',var_index)
                        if var_index >= 0:
                            if var_idx == var_index:
                                setattr(ext_object, var.strip(), raw_data[:, ...])
                        else:
                            setattr(ext_object, var.strip(), raw_data[:, var_idx, ...])
                        var_idx += 1
                else:
                    for var, var_idx in zip(plane_var_names,
                                                range(len(plane_var_names))):
                        setattr(ext_object, var.strip(), raw_data[:, var_idx, ...])

                ext_object.t = t
                self.t = t
                setattr(self, plane, ext_object)

        else:
            #Reading averages from hdf5 files
                # Get the averaged quantities.
            dim = read.dim(datadir=datadir)

            av_files = glob.glob(join(datadir,"averages","*.h5"))
            if len(av_files) == 0:
                raise RuntimeError(f"read.aver error: no averages files in {join(datadir,'averages')}")

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
                            warnings.warn(prefix+" does not match in avfile_list and plane_list")
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
                    raise RuntimeError(f"read.aver error: {plane_list} and {av_files} have no match.")
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
                    raise RuntimeError(f"read.aver error: {avfile_list} has no match in {av_files}")

            if not quiet:
                print(av_files_in)
            for av_file, plane in zip(av_files_in, plane_list):

                # Get the averaged quantities
                t, ext_object = self._read_h5_aver(
                        plane,
                        av_file,
                        var_names,
                        var_index,
                        iter_list,
                        iter_step,
                        time_range,
                        proc,
                        simdir,
                        precision=precision,
                        comp_time=comp_time,
                        quiet=quiet,
                    )

                ext_object.t = t
                self.t = t
                setattr(self, plane, ext_object)

    def _read_h5_aver(
        self,
        plane,
        av_file,
        var_names,
        var_index,
        iter_list,
        iter_step,
        time_range,
        proc,
        simdir,
        precision="f",
        comp_time=False,
        quiet=True,
    ):
        """
        Return the raw data and the time array.
        """

        import h5py

        if not quiet:
            print(av_file)
        with open(os.path.join(simdir, plane+"aver.in")) as file_id:
            variables = file_id.readlines()
        variables = [
            v.strip() for v in variables if v[0] != "#" and not v.isspace()
        ]  # Ignore commented variables and blank lines in the .in file.
        with h5py.File(av_file, "r") as tmp:
            n_times = len(tmp.keys()) - 1
            n_times = min(n_times,tmp['last'][()].item() + 1)
            start_time, end_time = 0, tmp[str(n_times-1)]['time'][()].item()
            if time_range is not None:
                start_time, end_time = time_range
            if iter_list:
                if len(iter_list) == 1:
                    n_times = min(n_times-1, iter_list[0])
                    itlist = range(n_times,n_times+1,iter_step)
                elif len(iter_list) == 2:
                    if iter_list[0] >= iter_list[1]:
                        warnings.warn("read.aver: iter_list pair must be list of integer - reading full series instead")
                    itlist = range(max(0,iter_list[0]),min(n_times,iter_list[1]),iter_step)
                else:
                    itlist = iter_list
                for it in itlist:
                    if not str(it) in tmp.keys():
                        itlist.remove(it)
                if not len(itlist) > 0:
                    warnings.warn("read.aver: iter_list has no match in {} keys - reading only {} instead".format(av_file,tmp['last'][0]))
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
                    raise RuntimeError(f"read.aver error: no data in {av_file} within time range {time_range}.")
                else:
                    itlist = tmplist
            # Determine the structure of the xy/xz/yz/y/z averages.
            if var_names is not None:
                var_names = [v.strip() for v in var_names]
                if isinstance(var_names, list):
                    var_names = var_names
                else:
                    var_names = [var_names]
                for var_name in var_names:
                    if not var_name in variables:
                        var_names.remove(var_name)
                if len(var_names) < 1:
                    warnings.warn("read.aver: var_names has no match in {} keys - reading all variables instead".format(av_file))
                    var_names = variables
            else:
                var_names = variables
            if "time" in var_names:
                var_names.remove('time')

            # This one will store the data.
            ext_object = getattr(self, plane, _Plane())
            if comp_time:
                start_time = time.time()
                data_shape = None
                for var in var_names:
                    if not data_shape:
                        data_shape = list()
                        data_shape.append(len(itlist))
                        if var in tmp[str(itlist[0])].keys():
                            dshape = tmp[str(itlist[0])+ "/" + var].shape
                        for dsh in dshape:
                            data_shape.append(dsh)
                    raw_data = np.zeros(data_shape, dtype=precision)
                    t = np.zeros(data_shape[0], dtype=precision)
                    for t_idx, tmp_idx in zip(range(len(itlist)),itlist):
                        t[int(t_idx)] = tmp[str(tmp_idx) + "/time"][()]

                        if var in tmp[str(tmp_idx)].keys():
                            raw_data[int(t_idx)] = tmp[str(tmp_idx) + "/" + var][()]

                    setattr(ext_object, var, raw_data)
            else:
                start_time = time.time()

                # testkey = list(tmp[str(itlist[0])].keys())[0] #Fred
                # testkey = var_names[0] #Kishore
                # data_shape = [len(itlist), *tmp[str(itlist[0])][testkey].shape]
                # t = np.zeros(data_shape[0], dtype=precision)
                # Kishore: Fred, you replaced var_names[0] by the above, but this breaks the reading of 2D averages (e.g. yaver). Can you please check if what I have now done below addresses your concern?
                if len(var_names) > 0:
                    data_shape = [len(itlist), *tmp[f"{itlist[0]}/{var_names[0]}"].shape]
                else:
                    #TODO (Kishore): should this be an error? I can't think of a legitimate situation where var_names is empty (see the `if var_names is not None:` block at the beginning of this function, where a length-zero var_names is explicitly overridden, and entries of var_names not in `variables` are explicitly removed).
                    pass
                for var in var_names:
                    #TODO (Kishore): should this be np.nan, so that there is a distinction between the average being zero and the average not being present for that particular time?
                    setattr(ext_object, var, np.zeros(data_shape, dtype=precision))

                t = np.zeros(len(itlist), dtype=precision)

                for t_idx, tmp_idx in enumerate(itlist):
                    t[t_idx] = tmp[f"{tmp_idx}/time"][()]
                    for var in var_names:
                        if var in tmp[str(tmp_idx)].keys():
                            getattr(ext_object, var)[t_idx] = tmp[f"{tmp_idx}/{var}"][()]

        end_time = time.time()-start_time
        if not quiet:
            print("{} object reading time {:.0f} seconds".format(plane,end_time))

        for k in ext_object.__dict__.keys():
            if k != "t":
                val = getattr(ext_object, k)
                if val.ndim == 3:
                    """
                    To be consistent with what happens with io_dist.
                    Axis ordering is now [t,x,z] for yaver, or [t,x,y] for zaver.
                    """
                    setattr(ext_object, k, val.swapaxes(1,2))

                    # 2025-Nov-28/Kishore: I suppose we can remove this warning
                    # after a year or so.
                    warnings.warn(
                        "Axis ordering for yaver and zaver has recently been changed;"
                        "it is now [t,x,z] for yaver, or [t,x,y] for zaver."
                        )

        return t, ext_object

    def _read_1d_aver(
        self,
        plane,
        datadir,
        aver_file_name,
        n_vars,
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

        from scipy.io import FortranFile

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
            vindex = list()
            if not isinstance(var_index,list):
                if var_index >= 0:
                    inx1 = var_index * pnu * pnv
                    inx2 = (var_index + 1) * pnu * pnv
                else:
                    inx1 = 0
                    inx2 = (n_vars + 1) * pnu * pnv
                vindex.append([inx1, inx2])
            else:
                for var_ind in var_index:
                    inx1 = var_ind * pnu * pnv
                    inx2 = (var_ind + 1) * pnu * pnv
                    vindex.append([inx1,inx2])
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
            #indices can be specified for subset of time series
            if iter_list is not None:
                plane_iter_list = iter_list
                liter = True
            else:
                plane_iter_list = [0,1,2,3]
                liter = False
            if time_range is not None:
                start_time, end_time = time_range
                ltime = True
            else:
                ltime = False
            #    # split by iteration overrules split by variable
            #    var_index = -1
            iiter = 0
            while True:
                if not liter:
                    plane_iter_list.append(iiter+1)
                try:
                    if iiter in plane_iter_list:
                        t_tmp = file_id.read_record(dtype=read_precision)[0].astype(precision)
                        if ltime:
                            if t_tmp >= start_time:
                                if t_tmp <= end_time:
                                    t.append(t_tmp)
                                    v_tmp = file_id.read_record(dtype=read_precision).astype(precision)
                                    for inx1, inx2 in vindex:
                                        proc_data.append(v_tmp[inx1:inx2])
                        else:
                            t.append(t_tmp)
                            v_tmp = file_id.read_record(dtype=read_precision).astype(precision)
                            for inx1, inx2 in vindex:
                                proc_data.append(v_tmp[inx1:inx2])
                    if iiter >= plane_iter_list[-1]:
                        # Finished reading.
                        break
                    iiter += 1
                except:
                    # Finished reading.
                    break
            file_id.close()
            # Reshape the proc data into [len(t), pnu, pnv].
            proc_data = np.array(proc_data, dtype=precision)
            if not isinstance(var_index, list):
                if var_index >= 0:
                    proc_data = proc_data.reshape([len(t), 1, pnv, pnu])
                else:
                    proc_data = proc_data.reshape([len(t), n_vars, pnv, pnu])
            else:
                    proc_data = proc_data.reshape([len(t), len(vindex), pnv, pnu])

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
                if not isinstance(var_index, list):
                    if var_index >= 0:
                        raw_data = np.zeros([len(t), 1, nv, nu], dtype=precision)
                    else:
                        raw_data = np.zeros([len(t), n_vars, nv, nu], dtype=precision)
                else:
                    raw_data = np.zeros([len(t), len(vindex), nv, nu], dtype=precision)
            raw_data[
                :, :, idx_v : idx_v + pnv, idx_u : idx_u + pnu
            ] = proc_data.copy()

        t = np.array(t, dtype=precision)
        raw_data = np.swapaxes(raw_data, 2, 3)

        return t, raw_data

    def _read_2d_aver(
        self,
        plane,
        datadir,
        aver_file_name,
        n_vars,
        var_names,
        var_index,
        iter_list,
        iter_step,
        time_range,
        proc,
        precision="f",
    ):
        """
        Read the xyaverages.dat, xzaverages.dat, yzaverages.dat
        Return the raw data and the time array.
        """

        if time_range is not None:
            warnings.warn("Averages._read_2d_aver: time_range is not implemented")

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
            raise RuntimeError(f"Error: There was a problem reading {aver_file_name} at line {line_idx}.\nCalculated values: n_vars = {n_vars}, nw = {nw}.\nAre these correct?")

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

@copy_docstring(Averages.read)
def aver(*args, **kwargs):
    """
    Wrapper for :py:meth:`Averages.read`
    """
    averages_tmp = Averages()
    averages_tmp.read(*args, **kwargs)
    return averages_tmp

class _Plane():
    """
    Used to store the averages in a particular plane
    """
    def keys(self):
        for i in self.__dict__.keys():
            if not i == "keys":
               print(i)
