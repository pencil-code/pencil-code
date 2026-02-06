# power.py
#
# Read the power spectra of the simulation.
#
"""
Contains the classes and methods to read the power spectra.
"""

from collections.abc import Iterable
import os
import numpy as np
from pencil import read
from pencil.util import ffloat, copy_docstring
import re
import warnings
import functools
import pathlib

class Power(object):
    """
    Power -- holds power spectrum data.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.t = []

    def keys(self):
        return list(self.__dict__.keys())

    def read(
            self,
            datadir="data",
            file_name=None,
            quiet=False,
            time_range=None,
            lazy=False,
            ):
        """
        read(datadir='data', file_name='', quiet=False)

        Read the power spectra.

        Parameters
        ----------
        datadir : string
            Directory where the data is stored.

        file_name : string
            Filename to read.
            If a filename is given, only that power spectrum is read.
            By default it reads all the power spectrum files.

        quiet : bool
            Flag for switching off output.

        lazy: bool
            If True, the HDF5 files will not be completely read into memory;
            rather, just the relevant portions will be read into memory and
            returned when the user applies indices.
            Default: False

        Returns
        -------
        Class containing the different power spectrum as attributes.

        Notes
        -----
        Use the attribute keys to get a list of attributes

        Examples
        --------
        >>> pw = pc.read.power()
        >>> pw.keys()
        t
        kin
        krms
        hel_kin
        """

        if file_name is not None:
            if not quiet:
                print("Reading only ", file_name)

            if os.path.isfile(os.path.join(datadir, file_name)):
                power_list, file_list = self._parse_filelist([file_name], quiet)
            else:
                raise ValueError(f"File {file_name} does not exist.")

        else:
            power_list, file_list = self._parse_filelist(os.listdir(datadir), quiet)

        # Read the power spectra.
        for power_name, file_name in zip(power_list, file_list):
            if not quiet:
                print(file_name)

            if re.match("power.*_xy.dat", file_name):
                self._read_power2d(power_name, file_name, datadir)
            elif re.match("power.*_xy.h5", file_name):
                if lazy:
                    self._read_power2d_hdf5_lazy(power_name, file_name, datadir)
                else:
                    self._read_power2d_hdf5(power_name, file_name, datadir)
            elif (
                file_name == "poweruz_x.dat"
                or file_name == "powerux_x.dat"
                or file_name == "poweruy_x.dat"
                ):
                self._read_power_1d(power_name, file_name, datadir)
            elif file_name == "power_krms.dat":
                self._read_power_krms(power_name, file_name, datadir)
            else:
                self._read_power(power_name, file_name, datadir)

        # 11-Dec-2025/Kishore: It is unclear what the point of time_range is;
        # 11-Dec-2025/Kishore: the slicing only happens after everything has
        # 11-Dec-2025/Kishore: been parsed and read into memory! Fred, do you
        # 11-Dec-2025/Kishore: find this option useful?
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
            ilist = list()
            for i, time in zip(range(self.t.size),self.t):
                if time >= start_time:
                    if time <= end_time:
                        ilist.append(i)
            for key in self.__dict__.keys():
                if not key=="krms":
                    tmp = self.__getattribute__(key)[ilist]
                    self.__delattr__(key)
                    setattr(self, key, tmp)

    def _read_power2d(self, power_name, file_name, datadir):
        """
        Handles output of power_xy subroutine.
        Axis order of self.power_name will be [t,ivec,z,ky,kx] or [t,z,ky,kx] or [t,z,k].
        """
        dim = read.dim(datadir=datadir)
        param = read.param(datadir=datadir)

        with open(os.path.join(datadir, file_name), "r") as f:
            _ = f.readline()  # ignore first line
            header = f.readline()

            # Get k vectors:
            if param.lintegrate_shell:
                nk = int(
                    header
                    .split()[header.split().index("k") + 1]
                    .split(")")[0][1:]
                    )
                k = []
                for _ in range(int(np.ceil(nk / 8))):
                    line = f.readline()
                    k.extend([float(j) for j in line.split()])
                k = np.array(k)
                self.k = k
            else:
                nkx = int(
                    header
                    .split()[header.split().index("k_x") + 1]
                    .split(")")[0][1:]
                    )
                kx = []
                for _ in range(int(np.ceil(nkx / 8))):
                    line = f.readline()
                    kx.extend([float(j) for j in line.split()])
                kx = np.array(kx)
                self.kx = kx

                nky = int(
                    header
                    .split()[header.split().index("k_y") + 1]
                    .split(")")[0][1:]
                    )
                ky = []
                for _ in range(int(np.ceil(nky / 8))):
                    line = f.readline()
                    ky.extend([float(j) for j in line.split()])
                ky = np.array(ky)
                self.ky = ky

                nk = nkx * nky

            # Now read z-positions, if any
            if param.lintegrate_z:
                nzpos = 1
            else:
                ini = f.tell()
                line = f.readline()
                if "z-pos" in line:
                    nzpos = int(re.search(r"\((\d+)\)", line)[1])
                    block_size = int(np.ceil(nzpos / 8))
                    zpos = []
                    for _ in range(block_size):
                        line = f.readline()
                        zpos.extend([ffloat(j) for j in line.split()])
                    self.zpos = np.array(zpos)
                else:
                    # there was no list of z-positions, so reset the position of the reader.
                    f.seek(ini)

                    nzpos = dim.nzgrid
                    grid = read.grid(datadir=datadir, trim=True, quiet=True)
                    self.zpos = grid.z

            # Now read the rest of the file
            time = []
            if param.lcomplex:
                power_array = np.array([], dtype=np.csingle)
            else:
                power_array = np.array([], dtype=np.single)

            if param.lintegrate_shell:
                block_size = np.ceil(nk / 8) * nzpos + 1
            else:
                block_size = np.ceil(int(nk * nzpos) / 8) + 1

            for line_idx, line in enumerate(f):
                if line_idx % block_size == 0:
                    time.append(self._parse_time_line(line))

                    power_array.resize([len(time), nzpos*nk])
                    ik = 0
                else:
                    lsp = line.strip().split()

                    if param.lcomplex:
                        # complex power spectrum
                        real = lsp[0::2]
                        imag = lsp[1::2]
                        for a, b in zip(real, imag):
                            power_array[-1,ik] = ffloat(a) + 1j * ffloat(b)
                            ik += 1
                    else:
                        for value_string in lsp:
                            power_array[-1,ik] = ffloat(value_string)
                            ik += 1

        time = np.array(time)

        if param.lintegrate_shell or (dim.nxgrid == 1 or dim.nygrid == 1):
            power_array = power_array.reshape([len(time), nzpos, nk])
        else:
            power_array = power_array.reshape([len(time), nzpos, nky, nkx])

        self.t = time.astype(np.single)
        self.nzpos = nzpos
        setattr(self, power_name, power_array)

    def _read_power2d_hdf5(self, power_name, file_name, datadir):
        """
        Handles HDF5 output of power_xy subroutine.
        Axis order of self.power_name will be [t,ivec,z,ky,kx] or [t,z,ky,kx].
        """
        import h5py

        param = read.param(datadir=datadir)

        with h5py.File(os.path.join(datadir, file_name)) as f:
            if param.lintegrate_shell:
                raise NotImplementedError
            elif param.lintegrate_z:
                raise NotImplementedError
            else:
                try:
                    version = f['metadata/output_version'][()]
                except KeyError:
                    #We don't know how to read this file
                    return

                self.kx = f['metadata/kx'][()]
                self.ky = f['metadata/ky'][()]
                self.zpos = f['metadata/z'][()]
                self.nzpos = len(self.zpos)

                nt = int(np.squeeze(f['last'][()]))
                time = np.empty([nt])
                power_shape = f[f"{nt}"]['data_re'].shape
                power_re = np.empty([nt, *power_shape])
                power_im = np.empty_like(power_re)

                for it in range(nt):
                    time[it] = self._h5_get(f, "time", it)
                    power_re[it] = self._h5_get(f, "data_re", it)
                    power_im[it] = self._h5_get(f, "data_im", it)

        self.t = time

        power_array = power_re + 1j*power_im
        if power_array.shape[1] == 1:
            power_array = np.squeeze(power_array, axis=1)
        setattr(self, power_name, power_array)

    def _read_power2d_hdf5_lazy(self, power_name, file_name, datadir):
        """
        Replacement for :py:meth:`Power._read_power2d_hdf5` that does not read
        the entire data into memory (since it can be hundreds of GB). Rather,
        the relevant part of the data is read into memory when the user indexes
        `self.*_xy`.
        Axis order of self.power_name will be [t,ivec,z,ky,kx] or [t,z,ky,kx].
        """
        import h5py

        param = read.param(datadir=datadir)

        with h5py.File(os.path.join(datadir, file_name)) as f:
            if param.lintegrate_shell:
                raise NotImplementedError
            elif param.lintegrate_z:
                raise NotImplementedError
            else:
                try:
                    version = f['metadata/output_version'][()]
                except KeyError:
                    #We don't know how to read this file
                    return

                self.kx = f['metadata/kx'][()]
                self.ky = f['metadata/ky'][()]
                self.zpos = f['metadata/z'][()]
                self.nzpos = len(self.zpos)

                nt = int(np.squeeze(f['last'][()]))
                time = np.empty([nt])
                power_shape = f[f"{nt}"]['data_re'].shape

                data_re_list = []
                data_im_list = []
                for it in range(nt):
                    time[it] = self._h5_get(f, "time", it)

                setattr(
                    self,
                    power_name,
                    _LazyPowerArray(f, nt),
                    )

        self.t = time

    def _h5_get(self, f, key, it, lazy=False):
        """
        The only purpose of this function is to mention the iteration number
        in the error message.

        Arguments:
            f: h5py.File instance
            key: str. Dataset name.
            it: int. Iteration number (zero-based indexing).
        """
        try:
            val = f[f"{it+1}/{key}"]

            if not lazy:
                val = val[()]
        except KeyError as e:
            e.add_note(f"iteration = {it+1}")
            raise e
        return val

    def _read_power_1d(self, power_name, file_name, datadir):
        """
        Handle output of subroutine power_1d
        """
        dim = read.dim(datadir=datadir)

        block_size = np.ceil(int(dim.nxgrid / 2) / 8.0) + 1

        time = []
        power_array = []
        with open(os.path.join(datadir, file_name), "r") as f:
            for line_idx, line in enumerate(f):
                if line_idx % block_size == 0:
                    time.append(self._parse_time_line(line))
                elif line.find(",") == -1:
                    # if the line does not contain ',', assume it represents a series of real numbers.
                    for value_string in line.strip().split():
                        power_array.append(float(value_string))
                else:
                    # Assume we have complex numbers.
                    for value_string in line.strip().split("( ")[1:]:
                        value_string = (
                            value_string.replace(")", "j")
                            .strip()
                            .replace(", ", "")
                            .replace(" ", "+")
                        )
                        power_array.append(complex(value_string))

        time = np.array(time)
        power_array = np.array(power_array).reshape([len(time), int(dim.nxgrid / 2)])
        self.t = time
        setattr(self, power_name, power_array)

    def _read_power_krms(self, power_name, file_name, datadir):
        """
        Read power_krms.dat.
        """
        nk = self._get_nk_xyz(datadir)
        block_size = np.ceil(nk/8)

        power_array = []
        with open(os.path.join(datadir, file_name), "r") as f:
            for line_idx, line in enumerate(f):
                # KG: Is this file expected to contain extra lines? If not, the if below can be removed.
                if line_idx < block_size:
                    for value_string in line.strip().split():
                        power_array.append(float(value_string))
        power_array = (
            np.array(power_array).reshape([nk]).astype(np.float32)
        )
        setattr(self, power_name, power_array)

    def _read_power(self, power_name, file_name, datadir):
        """
        Handles output of power subroutine.
        """
        nk = self._get_nk_xyz(datadir)
        block_size = np.ceil(nk/8) + 1

        time = []
        power_array = []
        with open(os.path.join(datadir, file_name), "r") as f:
            for line_idx, line in enumerate(f):
                if line_idx % block_size == 0:
                    time.append(self._parse_time_line(line))
                else:
                    for value_string in line.strip().split():
                        power_array.append(ffloat(value_string))

        # Reformat into arrays.
        time = np.array(time)
        power_array = (
            np.array(power_array)
            .reshape([len(time), nk])
            .astype(np.float32)
        )
        self.t = time.astype(np.float32)
        setattr(self, power_name, power_array)

    @functools.lru_cache(maxsize=128)
    def _get_nk_xyz(self, datadir):
        """
        See variable nk_xyz in power_spectrum.f90.

        NOTE: If you want to read output from non-cubic-box simulations run using older versions of Pencil where the number of k-vectors was always taken as nxgrid/2, you can do
        ```
        >>> class Power_wrong(pc.read.powers.Power):
        ...     def _get_nk_xyz(self, datadir):
        ...         dim = read.dim(datadir=datadir)
        ...         return int(dim.nxgrid/2)

        >>> p = Power_wrong()
        >>> p.read()
        ```
        """
        dim = read.dim(datadir=datadir)
        try:
            grid = read.grid(datadir=datadir, quiet=True)
        except FileNotFoundError:
            # KG: Handling this case because there is no grid.dat in `tests/input/serial-1/proc0` and we don't want the test to fail. Should we just drop this and add a grid.dat in the test input?
            warnings.warn("grid.dat not found. Assuming the box is cubical.")
            return int(dim.nxgrid/2)

        Lx = grid.Lx
        Ly = grid.Ly
        Lz = grid.Lz

        nx = dim.nx
        ny = dim.ny
        nz = dim.nz

        L_min = np.inf

        if nx != 1:
            L_min = min(L_min, Lx)
        if ny != 1:
            L_min = min(L_min, Lz)
        if nz != 1:
            L_min = min(L_min, Lz)
        if L_min == np.inf:
            L_min = 2*np.pi

        nk = np.inf
        if nx != 1:
            nk = min(nk, np.round(nx*L_min/(2*Lx)))
        if ny != 1:
            nk = min(nk, np.round(ny*L_min/(2*Ly)))
        if nz != 1:
            nk = min(nk, np.round(nz*L_min/(2*Lz)))
        if nk == np.inf:
            nk = 1

        return int(nk)

    def _parse_filelist(self, file_list_in, quiet):
        power_list = []
        file_list = []

        for file_name in file_list_in:
            fileext = file_name.split('.')[-1]
            if file_name[:5] == "power" and fileext in ["dat", "h5"]:
                if file_name[:6] == "power_":
                    power_name = file_name.split(".")[0][6:]
                else:
                    power_name = file_name.split(".")[0][5:]

                if not quiet:
                    print("appending", power_name)

                power_list.append(power_name)
                file_list.append(file_name)

        return power_list, file_list

    def _parse_time_line(self, line):
        """
        Since commit c9911621738739ac6847274ff59f3cd8bf9ac254, the 'time' line in
        the power spectrum files contains two quantities: the first one is the
        triggering quantity (which may be, e.g., the scale factor), and the second
        is the time.

        It is important to be able to read both the old format and the new format.
        """
        entries = [float(t) for t in line.strip().split()]
        return entries[-1]

class _m_LazyPowerArray:
    """
    Mixin to populate properties of _LazyPowerArray
    """
    @property
    def ndim(self):
        return len(self.shape)

    @property
    def dtype(self):
        return np.cdouble

class _LazyPowerArrayVD(_m_LazyPowerArray):
    """
    A container that reads the HDF5 power_xy data only when indexed.

    Example:
    >>> a = _LazyPowerArray(hdf5_file_handle, [hdf5_dset_1, hdf5_dset_2, ...], [hdf5_dset_1, hdf5_dset_2, ...]) #will not read anything into memory
    >>> a[2,:,64,32] #will read only the requested values into memory

    2025-Dec-13/Kishore: while this works fine with the samples in the auto-test,
    it does not with my large production runs (all the read data is np.nan even
    though the correct values are in the HDF5 files). I have thus replaced it by
    _LazyPowerArrayNoVD.
    """

    def __init__(self, h5file, nt):
        import h5py

        self._memfile = h5py.File.in_memory() #will hold virtual datasets

        for key in ["data_re", "data_im"]:
            shape = None
            for it in range(nt):
                try:
                    dset = h5file[f"{it+1}/{key}"]
                except KeyError as e:
                    e.add_note(f"key = '{key}'; iteration = {it+1}")
                    raise e

                if shape is None:
                    shape = (nt, *dset.shape)
                    layout = h5py.VirtualLayout(shape=shape, dtype=dset.dtype)
                else:
                    assert dset.shape == shape[1:]

                layout[it] = h5py.VirtualSource(dset)

            self._memfile.create_virtual_dataset(key, layout, fillvalue=np.nan)

        self._dset_re = self._memfile['data_re']
        self._dset_im = self._memfile['data_im']

        dset_shape = self._dset_re.shape
        if self._dset_re.shape[1] == 1:
            self.shape = (dset_shape[0], *dset_shape[2:])
        else:
            self.shape = dset_shape

    def __getitem__(self, k):
        if isinstance(k, (list, np.ndarray)):
            raise NotImplementedError("fancy indexing")

        if isinstance(k, Iterable):

            k = list(k)

            if self._dset_re.shape[1] == 1:
                #scalar dataset
                #this mirrors the removal of size-1 axes in _read_power2d
                if len(k) == 0:
                    k = [np.s_[:], 0]
                else:
                    k = [k[0], 0, *k[1:]]
        else:
            k = [k]

        if Ellipsis in k:
            raise NotImplementedError("use of `...`")

        if len(k) < 5:
            k = k + [np.s_[:]]*(5-len(k))
        elif len(k) > 5:
            raise ValueError("number of indices specified is more than the number of axes")

        #convert back to tuple since we don't want fancy indexing
        k = tuple(k)
        return self._dset_re[k] + 1j*self._dset_im[k]

    def __del__(self):
        if self._memfile:
            self._memfile.close()

class _LazyPowerArrayNoVD(_m_LazyPowerArray):
    """
    Variant of _LazyPowerArray that does not rely on virtual datasets.
    """
    def __init__(self, h5file, nt):
        self._file_path = pathlib.Path(h5file.filename).absolute()
        self._nt = nt

        dset_shape = h5file[f"{nt}/data_re"].shape
        if dset_shape[0] == 1:
            self._is_vec = False
            self.shape = (nt, *dset_shape[1:])
        else:
            self._is_vec = True
            self.shape = (nt, *dset_shape)

    def __getitem__(self, k):
        import h5py

        if isinstance(k, (list, np.ndarray)):
            raise NotImplementedError("fancy indexing")

        if isinstance(k, Iterable):

            k = list(k)

            if not self._is_vec:
                #scalar dataset
                #this mirrors the removal of size-1 axes in _read_power2d
                if len(k) == 0:
                    k = [np.s_[:], 0]
                else:
                    k = [k[0], 0, *k[1:]]
        else:
            k = [k]

        if Ellipsis in k:
            raise NotImplementedError("use of `...`")

        if len(k) < 5:
            k = k + [np.s_[:]]*(5-len(k))
        elif len(k) > 5:
            raise ValueError("number of indices specified is more than the number of axes")

        #convert back to tuple since we don't want fancy indexing
        k = tuple(k)
        ret = []
        t_sl = k[0]
        inds = k[1:]
        with h5py.File(self._file_path, mode='r') as f:
            for it in np.atleast_1d(range(self._nt)[t_sl]):
                ret.append(f[f"{it+1}/data_re"][inds] + 1j*f[f"{it+1}/data_im"][inds])
        return np.array(ret, dtype=ret[0].dtype)

_LazyPowerArray = _LazyPowerArrayNoVD

@copy_docstring(Power.read)
def power(*args, **kwargs):
    """
    Wrapper for :py:meth:`Power.read`
    """
    power_tmp = Power()
    power_tmp.read(*args, **kwargs)
    return power_tmp
