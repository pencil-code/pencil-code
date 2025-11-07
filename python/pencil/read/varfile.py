# varfile.py
#
# Read the var files.
# NB: the f array returned is C-ordered: f[nvar, nz, ny, nx]
#     NOT Fortran as in Pencil (& IDL):  f[nx, ny, nz, nvar]
#
# Authors:
# J. Oishi (joishi@amnh.org)
# T. Gastine (tgastine@ast.obs-mip.fr)
# S. Candelaresi (iomsn1@gmail.com).
# 20.10.30 IL modified: added comments in the definition
"""
Contains the read class for the VAR file reading,
some simulation attributes and the data cube.
"""

import os
from os.path import join, exists
import numpy as np
import warnings
import time
import re
import sys

#from scipy.io import FortranFile
from .fortran_file import FortranFileExt

from pencil.math import natural_sort
from pencil.math.derivatives import curl, curl2
from pencil import read
from pencil.sim import __Simulation__

def var(*args, **kwargs):
    """
    var(var_file='', datadir='data', proc=-1, ivar=-1, quiet=True,
        trimall=False, magic=None, sim=None, precision='f', flist=None,
        timing=True, fbloc=True, lvec=True, lonlyvec=False, lpersist=False,
        range_x=None, range_y=None, range_z=None,
        irange_x=None, irange_y=None, irange_z=None)

    Read VAR files from Pencil Code. If proc < 0, then load all data
    and assemble, otherwise load VAR file from specified processor.

    The file format written by output() (and used, e.g. in var.dat)
    consists of the following Fortran records:
    1. data(mx, my, mz, nvar)
    2. t(1), x(mx), y(my), z(mz), dx(1), dy(1), dz(1), deltay(1)
    Here nvar denotes the number of slots, i.e. 1 for one scalar field, 3
    for one vector field, 8 for var.dat in the case of MHD with entropy.
    but, deltay(1) is only there if lshear is on! need to know parameters.


    Parameters
    ----------
     var_file : string
         Name of the VAR file.
         If not specified, use var.dat (which is the latest snapshot of the fields)

     datadir : string
         Directory where the data is stored.

     proc : int
         Processor to be read. If -1 read all and assemble to one array.

     ivar : int
       Index of the VAR file, if var_file is not specified.

     quiet : bool
         Flag for switching off output.

     trimall : bool
         Trim the data cube to exclude ghost zones.

     magic : bool
         If present list of derived values to be computed from the data, e.g. B = curl(A).

     sim : pencil code simulation object
         Contains information about the local simulation.

     precision : string
         Float 'f', double 'd' or half 'half'.

     flist : list
         If present list of exclusive basic farrays to include

     timing : bool
         Report the time taken to create the obbject

     fbloc : bool
         If memory is restricted omit duplicate farray copy

     lvec : bool
         Combine components to form a vector

     lonlyvec : bool
         If memory is restricted omit components and provide only the vector

     lpersist : bool
         Read the persistent variables if they exist

     range_[xyz] : 2-tuple of real
         coordinate range selection for subdomain

     irange_[xyz] : 2-tuple of integer
         index range selection for subdomain


    Returns
    -------
    DataCube
        Instance of the pencil.read.var.DataCube class.
        All of the computed fields are imported as class members.

    Examples
    --------
    Read the latest var.dat file and print the shape of the uu array:
    >>> var = pc.read.var()
    >>> print(var.uu.shape)

    Read the VAR2 file, compute the magnetic field B = curl(A),
    the vorticity omega = curl(u) and remove the ghost zones:
    >>> var = pc.read.var(var_file='VAR2', magic=['bb', 'vort'], trimall=True)
    >>> print(var.bb.shape)
    """

    started = None

    for a in args:
        if isinstance(a, __Simulation__):
            started = a.started()
            break

    if "sim" in kwargs.keys():
        # started = kwargs['sim'].started()

        started = True
    elif "datadir" in kwargs.keys():
        if exists(join(kwargs["datadir"], "time_series.dat")):
            started = True
    else:
        if exists(join("data", "time_series.dat")):
            started = True

    if not started:
        if "ivar" in kwargs:
            if kwargs["ivar"] != 0:
                print("ERROR: Simulation has not yet started. There are no var files.")
                return False

    var_tmp = DataCube()
    var_tmp.read(*args, **kwargs)
    return var_tmp


class DataCube(object):
    """
    DataCube -- holds Pencil Code VAR file data.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.t = 0.0
        self.dx = 1.0
        self.dy = 1.0
        self.dz = 1.0
        self.x = None
        self.y = None
        self.z = None
        self.f = None
        self.l1 = None
        self.l2 = None
        self.m1 = None
        self.m2 = None
        self.n1 = None
        self.n2 = None
        self.magic = None

    def keys(self):
        for i in self.__dict__.keys():
            print(i)

    def read(
        self,
        var_file="",
        datadir="data",
        proc=-1,
        ivar=-1,
        quiet=True,
        trimall=False,
        magic=None,
        sim=None,
        precision="d",
        lpersist=False,
        dtype=np.float64,
        flist=None,
        timing=True,
        fbloc=True,
        lvec=True,
        lonlyvec=False,
        range_x=None,
        range_y=None,
        range_z=None,
        irange_x=None,
        irange_y=None,
        irange_z=None
    ):
        """
        read(var_file='', datadir='data', proc=-1, ivar=-1, quiet=True,
             trimall=False, magic=None, sim=None, precision='d',
             lpersist=False, dtype=np.float64, flist=None,
             timing=True, fbloc=True, lvec=True, lonlyvec=False,
             range_x=None, range_y=None, range_z=None,
             irange_x=None, irange_y=None, irange_z=None)

        Read VAR files from Pencil Code. If proc < 0, then load all data
        and assemble, otherwise load VAR file from specified processor.

        The file format written by output() (and used, e.g. in var.dat)
        consists of the followinig Fortran records:
        1. data(mx, my, mz, nvar)
        2. t(1), x(mx), y(my), z(mz), dx(1), dy(1), dz(1), deltay(1)
        Here nvar denotes the number of slots, i.e. 1 for one scalar field, 3
        for one vector field, 8 for var.dat in the case of MHD with entropy.
        but, deltay(1) is only there if lshear is on! need to know parameters.

        Parameters
        ----------
         var_file : string
             Name of the VAR file.
             If not specified, use var.dat (which is the latest snapshot of the fields)

         datadir : string
             Directory where the data is stored.

         proc : int
             Processor to be read. If -1 read all and assemble to one array.

         ivar : int
           Index of the VAR file, if var_file is not specified.

         quiet : bool
             Flag for switching off output.

         trimall : bool
             Trim the data cube to exclude ghost zones.

         magic : list of strings
             If present list of derived values to be computed from the data, e.g. B = curl(A).

         sim : pencil code simulation object
             Contains information about the local simulation.

         precision : string
             Float 'f', double 'd' or half 'half'.

         flist : list
             If present list of exclusive basic farrays to include

         timing : bool
             Report the time taken to create the obbject

         fbloc : bool
             If memory is restricted omit duplicate farray copy

         lvec : bool
             Combine components to form a vector

         lonlyvec : bool
             If memory is restricted omit components and provide only the vector

         lpersist : bool
             Read the persistent variables if they exist

         range_[xyz] : 2-tuple of real
             coordinate range selection for subdomain (closed interval).
             If trimall=F, ghost zones will be included on either side of the interval.

         irange_[xyz] : 2-tuple of integer
             index range selection for subdomain (closed interval).
             Note that this index is in the full array (including ghost zones).
             If trimall=F, ghost zones will be included on either side of the interval.

        Returns
        -------
        DataCube
            Instance of the pencil.read.var.DataCube class.
            All of the computed fields are imported as class members.

        Examples
        --------
        Read the latest var.dat file and print the shape of the uu array:
        >>> var = pc.read.var()
        >>> print(var.uu.shape)

        Read the VAR2 file, compute the magnetic field B = curl(A),
        the vorticity omega = curl(u) and remove the ghost zones:
        >>> var = pc.read.var(var_file='VAR2', magic=['bb', 'vort'], trimall=True)
        >>> print(var.bb.shape)
        """

        if precision=="h":
            precision = "half"
        if timing:
            start_time = time.time()

        if sim is None:
            datadir = os.path.expanduser(datadir)
            dim = read.dim(datadir, proc)
            param = read.param(datadir=datadir, quiet=quiet, conflicts_quiet=True)
            index = read.index(datadir=datadir)

            try:
                grid = read.grid(datadir=datadir, quiet=True, proc=proc)
            except FileNotFoundError:
                # KG: Handling this case because there is no grid.dat in `tests/input/serial-1/proc0` and we don't want the test to fail. Should we just drop this and add a grid.dat in the test input?
                warnings.warn("Grid.dat not found. Assuming the grid is uniform.")
                grid = None
        else:
            datadir = os.path.expanduser(sim.datadir)
            dim = sim.dim
            param = read.param(datadir=datadir, quiet=True, conflicts_quiet=True)
            index = read.index(datadir=datadir)
            grid = read.grid(datadir=datadir, quiet=True) # we can't use sim.grid because we want the untrimmed one

        if param.io_strategy[-1] == '/': param.io_strategy = param.io_strategy[:-1]    # because of Python 3.10.10 bug!

        if var_file[0:2].lower() == "og":
            dim = read.ogdim(datadir, proc)
        elif var_file[0:4] == "VARd":
            dim = read.dim(datadir, proc, down=True)
            warnings.warn(
                "Reading downsampled grid is currently not implemented. Assuming the grid is uniform."
            )
            grid = None

        # Used later on to support the case where only some of the variables were written into the snapshots.
        index_max = dim.mvar + dim.maux

        run2D = param.lwrite_2d
        if param.io_strategy == "HDF5":
            grid = self._read_hdf5(
                grid=grid,
                dim=dim,
                param=param,
                var_file=var_file,
                datadir=datadir,
                precision=precision,
                quiet=quiet,
                lpersist=lpersist,
                ivar=ivar,
                irange_x=irange_x,
                irange_y=irange_y,
                irange_z=irange_z,
                range_x=range_x,
                range_y=range_y,
                range_z=range_z,
                dtype=dtype,
                index=index,
                )
        elif param.io_strategy == "dist":
            if (
                (range_x is not None) or
                (range_y is not None) or
                (range_z is not None) or
                (irange_x is not None) or
                (irange_y is not None) or
                (irange_z is not None)
                ):
                raise NotImplementedError("subdomains when IO = io_dist")

            grid = self._read_io_dist(
                grid=grid,
                dim=dim,
                param=param,
                proc=proc,
                ivar=ivar,
                var_file=var_file,
                datadir=datadir,
                precision=precision,
                quiet=quiet,
                lpersist=lpersist,
                )
        else:
            raise NotImplementedError(
                "IO strategy {} not supported by the Python module.".format(
                    param.io_strategy
                )
            )

        aatest = []
        uutest = []
        for key in index.__dict__.keys():
            if "aatest" in key:
                aatest.append(key)
            if "uutest" in key:
                uutest.append(key)
        if magic is not None:
            """
            In the functions curl and curl2, the arguments (dx,dy,dz,x,y) are ignored when grid is not None. Nevertheless, we pass them below to take care of the case where the user is trying to read a snapshot without the corresponding grid.dat being present (such as in the test test_read_var).
            """
            if "bb" in magic:
                # Compute the magnetic field before doing trimall.
                aa = self.f[index.ax - 1 : index.az, ...]
                self.bb = curl(
                        aa,
                        dx=self.dx,
                        dy=self.dy,
                        dz=self.dz,
                        x=self.x,
                        y=self.y,
                        run2D=run2D,
                        coordinate_system=param.coord_system,
                        grid=grid,
                )
                if trimall:
                    self.bb = self._trim(self.bb, dim, run2D)
            if "bbtest" in magic:
                if param.io_strategy == "HDF5":
                    # Compute the magnetic field before doing trimall.
                    for j in range(int(len(aatest) / 3)):
                        key = aatest[j*3][:-1]
                        value = index.__dict__[aatest[j*3]]
                        aa = self.f[value - 1 : value + 2, ...]
                        bb = curl(
                                aa,
                                dx=self.dx,
                                dy=self.dy,
                                dz=self.dz,
                                x=self.x,
                                y=self.y,
                                run2D=run2D,
                                coordinate_system=param.coord_system,
                                grid=grid,
                        )
                        if trimall:
                            setattr(
                                self,
                                "bb"+key[2:],
                                self._trim(bb, dim, run2D),
                                )
                        else:
                            setattr(self,"bb"+key[2:],bb)
                else:
                    if hasattr(index, "aatest1"):
                        naatest = int(len(aatest) / 3)
                        for j in range(0, naatest):
                            key = "aatest" + str(np.mod(j + 1, naatest))
                            value = index.__dict__["aatest1"] + 3 * j
                            aa = self.f[value - 1 : value + 2, ...]
                            bb = curl(
                                    aa,
                                    dx=self.dx,
                                    dy=self.dy,
                                    dz=self.dz,
                                    x=self.x,
                                    y=self.y,
                                    run2D=run2D,
                                    coordinate_system=param.coord_system,
                                    grid=grid,
                            )
                            if trimall:
                                setattr(
                                    self,
                                    "bb"+key[2:],
                                    self._trim(bb, dim, run2D),
                                    )
                            else:
                                setattr(self,"bb"+key[2:],bb)
            if "jj" in magic:
                # Compute the electric current field before doing trimall.
                aa = self.f[index.ax - 1 : index.az, ...]
                self.jj = curl2(
                        aa,
                        dx=self.dx,
                        dy=self.dy,
                        dz=self.dz,
                        x=self.x,
                        y=self.y,
                        coordinate_system=param.coord_system,
                        grid=grid,
                )
                if trimall:
                    self.jj = self._trim(self.jj, dim, run2D)
            if "vort" in magic:
                # Compute the vorticity field before doing trimall.
                uu = self.f[index.ux - 1 : index.uz, ...]
                self.vort = curl(
                        uu,
                        dx=self.dx,
                        dy=self.dy,
                        dz=self.dz,
                        x=self.x,
                        y=self.y,
                        run2D=run2D,
                        coordinate_system=param.coord_system,
                        grid=grid,
                )
                if trimall:
                    self.vort = self._trim(self.vort, dim, run2D)

        # Trim the ghost zones of the global f-array if asked.
        if trimall:
            self.x = self.x[dim.nghostx:-dim.nghostx]
            self.y = self.y[dim.nghosty:-dim.nghosty]
            self.z = self.z[dim.nghostz:-dim.nghostz]
            self.f = self._trim(self.f, dim, run2D)
        else:
            self.l1 = dim.l1
            self.l2 = dim.l2 + 1
            self.m1 = dim.m1
            self.m2 = dim.m2 + 1
            self.n1 = dim.n1
            self.n2 = dim.n2 + 1

        # Assign an attribute to self for each variable defined in
        # 'data/index.pro' so that e.g. self.ux is the x-velocity
        for key in index.__dict__.keys():
            if (
                key != "global_gg"
                and key != "keys"
                and "aatest" not in key
                and "uutest" not in key
            ):
                value = index.__dict__[key]
                if value <= index_max:
                    setattr(self, key, self.f[value - 1, ...])
        # Special treatment for vector quantities.
        if hasattr(index, "ux") and index.uz <= index_max:
            setattr(self, "uu", self.f[index.ux - 1 : index.uz, ...])
        if hasattr(index, "ax") and index.az <= index_max:
            setattr(self, "aa", self.f[index.ax - 1 : index.az, ...])
        if hasattr(index, "uu_sph") and index.uu_sphz <= index_max:
            self.uu_sph = self.f[index.uu_sphx - 1 : index.uu_sphz, ...]
        if hasattr(index, "bb_sph") and index.bb_sphz <= index_max:
            self.bb_sph = self.f[index.bb_sphx - 1 : index.bb_sphz, ...]
        # Special treatment for test method vector quantities.
        # Note index 1,2,3,...,0 last vector may be the zero field/flow
        if param.io_strategy != "HDF5":
            if hasattr(index, "aatest1"):
                naatest = int(len(aatest) / 3)
                for j in range(0, naatest):
                    key = "aatest" + str(np.mod(j + 1, naatest))
                    value = index.__dict__["aatest1"] + 3 * j
                    setattr(self, key, self.f[value - 1 : value + 2, ...])
            if hasattr(index, "uutest1"):
                nuutest = int(len(uutest) / 3)
                for j in range(0, nuutest):
                    key = "uutest" + str(np.mod(j + 1, nuutest))
                    value = index.__dict__["uutest"] + 3 * j
                    setattr(self, key, self.f[value - 1 : value + 2, ...])
        else:
            #Dummy operation to be corrected
            for j in range(int(len(aatest) / 3)):
                key = aatest[j*3][:-1]
                value = index.__dict__[aatest[j*3]]
                setattr(self, key, self.f[value - 1 : value + 2, ...])
            for j in range(int(len(uutest) / 3)):
                key = uutest[j*3][:-1]
                value = index.__dict__[uutest[j*3]]
                setattr(self, key, self.f[value - 1 : value + 2, ...])

        # Do the rest of magic after the trimall (i.e. no additional curl.)
        self.magic = magic
        if self.magic is not None:
            self.magic_attributes(param, dtype=dtype)
        if timing:
            print("object completed in {:.2f} seconds.".format(time.time()-start_time))

    def __natural_sort(self, procs_list):
        """
        Sort array in a more natural way, e.g. 9VAR < 10VAR
        """
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", key)]
        return sorted(procs_list, key=alphanum_key)

    def magic_attributes(self, param, dtype=np.float64):
        """
        Compute some additional 'magic' quantities.
        """
        for field in self.magic:
            if field == "rho" and not hasattr(self, "rho"):
                if hasattr(self, "lnrho"):
                    setattr(self, "rho", np.exp(self.lnrho))
                else:
                    raise AttributeError("Problem in magic: lnrho is missing")

            if field == "tt" and not hasattr(self, "tt"):
                if hasattr(self, "lnTT"):
                    tt = np.exp(self.lnTT)
                    setattr(self, "tt", tt)
                else:
                    if hasattr(self, "ss"):
                        if hasattr(self, "lnrho"):
                            lnrho = self.lnrho
                        elif hasattr(self, "rho"):
                            lnrho = np.log(self.rho)
                        else:
                            raise AttributeError(
                                "Problem in magic: missing rho or" + " lnrho variable"
                            )
                        cp = param.cp
                        gamma = param.gamma
                        cs20 = param.cs0 ** 2
                        lnrho0 = np.log(param.rho0)
                        lnTT0 = np.log(cs20 / (cp * (gamma - 1.0)))
                        lnTT = (
                            lnTT0
                            + gamma / cp * self.ss
                            + (gamma - 1.0) * (lnrho - lnrho0)
                        )
                        setattr(self, "tt", np.exp(lnTT))
                    else:
                        raise AttributeError("Problem in magic: ss is missing ")

            if field == "ss" and not hasattr(self, "ss"):
                cp = param.cp
                gamma = param.gamma
                cs20 = param.cs0 ** 2
                lnrho0 = np.log(param.rho0)
                lnTT0 = np.log(cs20 / (cp * (gamma - 1.0)))
                if hasattr(self, "lnTT"):
                    setattr(
                        self,
                        "ss",
                            cp
                            / gamma
                            * (
                                self.lnTT
                                - lnTT0
                                - (gamma - 1.0) * (self.lnrho - lnrho0)
                            )
                    )
                elif hasattr(self, "tt"):
                    setattr(
                        self,
                        "ss",
                            cp
                            / gamma
                            * (
                                np.log(self.tt)
                                - lnTT0
                                - (gamma - 1.0) * (self.lnrho - lnrho0)
                            )
                    )
                else:
                    raise AttributeError("Problem in magic: missing lnTT or tt")

            if field == "pp" and not hasattr(self, "pp"):
                cp = param.cp
                gamma = param.gamma
                cv = cp / gamma
                lnrho0 = np.log(param.rho0)
                if hasattr(self, "lnrho"):
                    lnrho = self.lnrho
                elif hasattr(self, "rho"):
                    lnrho = np.log(self.rho)
                else:
                    raise AttributeError("Problem in magic: missing rho or lnrho variable")
                if hasattr(self, "ss"):
                    cs20 = param.cs0 ** 2
                    lnTT0 = np.log(cs20 / (cp * (gamma - 1.0)))
                    setattr(self, "pp",
                           (cp - cv) * np.exp(lnTT0 + gamma / cp * self.ss
                                    + gamma * lnrho - (gamma - 1.0) * lnrho0))
                elif hasattr(self, "lntt"):
                    setattr(self, "pp", (cp - cv) * np.exp(self.lnTT + lnrho))
                elif hasattr(self, "tt"):
                    setattr(self, "pp", (cp - cv) * self.TT * np.exp(lnrho))
                else:
                    raise AttributeError("Problem in magic: missing ss or lntt or tt")

    def _get_persist_iodist(self, infile, precision, quiet):
        """An open Fortran file potentially containing persistent variables appended
        to the f array and grid data are read from the first proc data

        Record types provide the labels and id record for the peristent
        variables in the depricated fortran binary format
        """
        record_types = {}
        for key in read.record_types.keys():
            if read.record_types[key][1] == "d":
                record_types[key] = (read.record_types[key][0], precision)
            else:
                record_types[key] = read.record_types[key]

        try:
            tmp_id = infile.read_record("h")
        except:
            return -1
        block_id = 0
        pers_obj = _Persist()
        for i in range(2000):
            i += 1
            tmp_id = infile.read_record("h")
            block_id = tmp_id[0]
            if block_id == 2000:
                break
            for key in record_types.keys():
                #Kishore: DANGER: there is a wrong assumption here that persistent variables must be scalars. A counter-example is forcing_location.
                if record_types[key][0] == block_id:
                    tmp_val = infile.read_record(record_types[key][1])
                    pers_obj.__setattr__(key, tmp_val[0])
                    if not quiet:
                        print(key, record_types[key][0], record_types[key][1], tmp_val)

        self.persist = pers_obj

    def _read_hdf5(self, grid, dim, param, var_file, datadir, precision, quiet, lpersist, ivar, irange_x, irange_y, irange_z, range_x, range_y, range_z, dtype, index):
        import h5py

        if param.lwrite_aux:
            total_vars = dim.mvar + dim.maux
        else:
            total_vars = dim.mvar

        run2D = param.lwrite_2d

        if not var_file:
            if ivar < 0:
                var_file = "var.h5"
            else:
                var_file = "VAR" + str(ivar) + ".h5"

        file_name = os.path.join(datadir, "allprocs", var_file)

        with h5py.File(file_name, "r") as tmp:
            x = (tmp["grid/x"][:]).astype(precision)
            irange_x, mx, x = self._handle_range(range_x, irange_x, x, dim.nghostx)

            y = (tmp["grid/y"][:]).astype(precision)
            irange_y, my, y = self._handle_range(range_y, irange_y, y, dim.nghosty)

            z = (tmp["grid/z"][:]).astype(precision)
            irange_z, mz, z = self._handle_range(range_z, irange_z, z, dim.nghostz)

            if grid != None:
                grid.restrict(irange_x,irange_y,irange_z)

            # Set up the global array.
            if run2D:
                if dim.ny == 1:
                    self.f = np.zeros((total_vars, mz, mx), dtype=dtype)
                elif dim.nz == 1:
                    self.f = np.zeros((total_vars, my, mx), dtype=dtype)
                else:
                    self.f = np.zeros((total_vars, mz, my), dtype=dtype)
            else:
                self.f = np.zeros((total_vars, mz, my, mx), dtype=dtype)

            for key in tmp["data"].keys():
                if key in index.__dict__.keys():
                    self.f[index.__getattribute__(key) - 1, :, :, :] = dtype(
                            tmp["data/" + key][irange_z[0]:irange_z[1],
                                                irange_y[0]:irange_y[1],
                                                irange_x[0]:irange_x[1]]
                    )
            t = (tmp["time"][()]).astype(precision)
            dx = (tmp["grid/dx"][()]).astype(precision)
            dy = (tmp["grid/dy"][()]).astype(precision)
            dz = (tmp["grid/dz"][()]).astype(precision)
            if param.lshear:
                deltay = (tmp["persist/shear_delta_y"][(0)]).astype(precision)
            if lpersist:
                pers_obj = _Persist()
                nprocs = dim.nprocx * dim.nprocy * dim.nprocz
                for key in tmp["persist"].keys():
                    val = tmp["persist"][key][()]
                    #Note that persistent variables need not be scalars (e.g. forcing_location)
                    val_local = np.split(val, nprocs)[0].astype(precision)
                    if len(val_local) == 1:
                        val_local = val_local[0]
                    setattr(pers_obj, key, val_local)
                self.persist = pers_obj

        self.x = x
        self.y = y
        self.z = z
        self.t = t
        self.dx = dx
        self.dy = dy
        self.dz = dz
        if param.lshear:
            self.deltay = deltay

        return grid

    def _read_io_dist(self, grid, dim, param, proc, ivar, var_file, datadir, precision, quiet, lpersist):
        if param.lwrite_aux:
            total_vars = dim.mvar + dim.maux
        else:
            total_vars = dim.mvar

        run2D = param.lwrite_2d

        if dim.precision == "D":
            read_precision = "d"
        else:
            read_precision = "f"

        if not var_file:
            if ivar < 0:
                var_file = "var.dat"
            else:
                var_file = "VAR" + str(ivar)

        if proc < 0:
            proc_dirs = self.__natural_sort(
                filter(lambda s: s.startswith("proc"), os.listdir(datadir))
            )
            if proc_dirs.count("proc_bounds.dat") > 0:
                proc_dirs.remove("proc_bounds.dat")
            if param.lcollective_io:
                # A collective IO strategy is being used
                proc_dirs = ["allprocs"]
        #                else:
        #                    proc_dirs = proc_dirs[::dim.nprocx*dim.nprocy]
        else:
            proc_dirs = ["proc" + str(proc)]

        # Set up the global array.
        if not run2D:
            self.f = np.zeros((total_vars, dim.mz, dim.my, dim.mx), dtype=precision)
        else:
            if dim.ny == 1:
                self.f = np.zeros((total_vars, dim.mz, dim.mx), dtype=precision)
            else:
                self.f = np.zeros((total_vars, dim.my, dim.mx), dtype=precision)

        x = np.zeros(dim.mx, dtype=precision)
        y = np.zeros(dim.my, dtype=precision)
        z = np.zeros(dim.mz, dtype=precision)

        for directory in proc_dirs:
            if not param.lcollective_io:
                proc = int(directory[4:])
                if var_file[0:2].lower() == "og":
                    procdim = read.ogdim(datadir, proc)
                else:
                    if var_file[0:4] == "VARd":
                        procdim = read.dim(datadir, proc, down=True)
                    else:
                        procdim = read.dim(datadir, proc)
                if not quiet:
                    print(  "Reading data from processor"
                            + " {0} of {1} ...".format(proc, len(proc_dirs)))

            else:
                # A collective IO strategy is being used
                procdim = dim
            #                else:
            #                    procdim.mx = dim.mx
            #                    procdim.my = dim.my
            #                    procdim.nx = dim.nx
            #                    procdim.ny = dim.ny
            #                    procdim.ipx = dim.ipx
            #                    procdim.ipy = dim.ipy

            mxloc = procdim.mx
            myloc = procdim.my
            mzloc = procdim.mz

            # Read the data: f-array
            file_name = os.path.join(datadir, directory, var_file)
            infile = FortranFileExt(file_name,header_dtype=np.int32)
            if not run2D:
                f_loc = (infile.read_record(dtype=read_precision)).astype(precision)
                f_loc = f_loc.reshape((-1, mzloc, myloc, mxloc))
            #    print(proc,f_loc.shape,f_loc.dtype)
            else:
                if dim.ny == 1:
                    f_loc = (infile.read_record(dtype=read_precision)).astype(precision)
                    f_loc = f_loc.reshape((-1, mzloc, mxloc))
                else:
                    f_loc = (infile.read_record(dtype=read_precision)).astype(precision)
                    f_loc = f_loc.reshape((-1, myloc, mxloc))

            # Read the data: time, coordinates, etc.
            raw_etc = (infile.read_record(dtype=read_precision)).astype(precision)

            # Read the data: persistent variables
            if lpersist and directory==proc_dirs[0]:
                self._get_persist_iodist(infile=infile, precision=read_precision, quiet=quiet)
            infile.close()

            t = raw_etc[0]
            x_loc = raw_etc[1 : mxloc + 1]
            y_loc = raw_etc[mxloc + 1 : mxloc + myloc + 1]
            z_loc = raw_etc[mxloc + myloc + 1 : mxloc + myloc + mzloc + 1]
            if param.lshear:
                shear_offset = 1
                deltay = raw_etc[-1]
            else:
                shear_offset = 0

            dx = raw_etc[-3 - shear_offset]
            dy = raw_etc[-2 - shear_offset]
            dz = raw_etc[-1 - shear_offset]

            if len(proc_dirs) > 1:
                # Calculate where the local processor will go in
                # the global array.
                #
                # Don't overwrite ghost zones of processor to the
                # left (and accordingly in y and z direction -- makes
                # a difference on the diagonals)
                #
                # Recall that in NumPy, slicing is NON-INCLUSIVE on
                # the right end, ie, x[0:4] will slice all of a
                # 4-digit array, not produce an error like in idl.

                if procdim.ipx == 0:
                    i0x = 0
                    i1x = i0x + procdim.mx
                    i0xloc = 0
                    i1xloc = procdim.mx
                else:
                    i0x = procdim.ipx * procdim.nx + procdim.nghostx
                    i1x = i0x + procdim.mx - procdim.nghostx
                    i0xloc = procdim.nghostx
                    i1xloc = procdim.mx

                if procdim.ipy == 0:
                    i0y = 0
                    i1y = i0y + procdim.my
                    i0yloc = 0
                    i1yloc = procdim.my
                else:
                    i0y = procdim.ipy * procdim.ny + procdim.nghosty
                    i1y = i0y + procdim.my - procdim.nghosty
                    i0yloc = procdim.nghosty
                    i1yloc = procdim.my

                if procdim.ipz == 0:
                    i0z = 0
                    i1z = i0z + procdim.mz
                    i0zloc = 0
                    i1zloc = procdim.mz
                else:
                    i0z = procdim.ipz * procdim.nz + procdim.nghostz
                    i1z = i0z + procdim.mz - procdim.nghostz
                    i0zloc = procdim.nghostz
                    i1zloc = procdim.mz

                x[i0x:i1x] = x_loc[i0xloc:i1xloc]
                y[i0y:i1y] = y_loc[i0yloc:i1yloc]
                z[i0z:i1z] = z_loc[i0zloc:i1zloc]

                if not run2D:
                    self.f[:, i0z:i1z, i0y:i1y, i0x:i1x] = f_loc[
                        :, i0zloc:i1zloc, i0yloc:i1yloc, i0xloc:i1xloc
                    ]
                else:
                    if dim.ny == 1:
                        self.f[:, i0z:i1z, i0x:i1x] = f_loc[
                            :, i0zloc:i1zloc, i0xloc:i1xloc
                        ]
                    else:
                        self.f[i0z:i1z, i0y:i1y, i0x:i1x] = f_loc[
                            i0zloc:i1zloc, i0yloc:i1yloc, i0xloc:i1xloc
                        ]
            else:                    # reading from a single processor
                self.f = f_loc
                x = x_loc
                y = y_loc
                z = z_loc

        self.x = x
        self.y = y
        self.z = z
        self.t = t
        self.dx = dx
        self.dy = dy
        self.dz = dz
        if param.lshear:
            self.deltay = deltay

    def _parse_range(self, rang, irang, coords, nghost):
        if (rang is not None) and (len(rang) != 2):
            raise ValueError
        if (irang is not None) and (len(irang) != 2):
            raise ValueError

        ind_min = nghost
        ind_maxp1 = len(coords)-nghost
        if rang is not None:
            [inds] = np.nonzero( (coords>=rang[0]) & (coords<=rang[-1]) )
            irang = (max(inds[0],ind_min), min(inds[-1]+1,ind_maxp1))
        elif irang is not None:
            irang = (max(irang[0],ind_min), min(irang[1]+1,ind_maxp1))
        else:
            irang = (ind_min,ind_maxp1)

        #Include ghost zones as well for proper computation of magic variables
        irang = (irang[0]-nghost, irang[1]+nghost)

        return irang

    def _handle_range(self, rang, irang, coords, nghost):
        irang = self._parse_range(rang, irang, coords, nghost)
        m = irang[1] - irang[0]
        coords = coords[irang[0]:irang[1]]
        return irang, m, coords

    def _trim(self, arr, dim, run2D):
        sl_tr_x = slice(dim.nghostx, -dim.nghostx)
        sl_tr_y = slice(dim.nghosty, -dim.nghosty)
        sl_tr_z = slice(dim.nghostz, -dim.nghostz)

        if (arr.ndim == 3) or (arr.ndim == 4):
            if run2D:
                if dim.nz == 1:
                    return arr[..., sl_tr_y, sl_tr_x]
                else:
                    return arr[..., sl_tr_z, sl_tr_x]
            else:
                return arr[..., sl_tr_z, sl_tr_y, sl_tr_x]
        else:
            raise NotImplementedError

class _Persist():
    """
    Used to store the persistent variables
    """
    def keys(self):
        for i in self.__dict__.keys():
            if not i == "keys":
               print(i)
