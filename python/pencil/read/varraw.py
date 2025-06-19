# varraw.py
#
# Read the Pencil Code varfiles.
# IMPORTANT: the f array returned is Fortran-ordered: f[nx, ny, nz, nvar]
#           NOT C-ordered. Please transpose if required.
# Implementation of the IDL routine pc_read_var_raw.pro,
# with many features adopted from varfile.py
#
# Author:
# Sankalp Srivastava (sankalp.srivastava@iiap.res.in)
# Added 14-Jun-2025
"""
Contains classes and methods to read the varfile data. Fortran binary varfiles for all Pencil code I/O schemes supported.
Selective reading of variables supported. The section with vector_groups is incomplete. Please add variable names as necessary.
Reading reduced varfiles from 'data/reduced' not implemented, since that would require the same to be implemented in 
at least dims.py as well (and in grids.py also for completeness).
swap_endian option has not been tested.
Persistent variables not yet implemented.
Downsampled snapshots not implemented. OGVAR files not implemented.
IMPORTANT: the f array returned is Fortran-ordered: f[nx, ny, nz, nvar]. Please transpose if required.
"""


def varraw(*args, **kwargs):
    """
    varraw(var_file="", datadir="data", proc=-1, ivar=-1, var_list=None, quiet=False, 
    precision="d", trimall=False, swap_endian=False)

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
         Processor to be read. If -1, read all and assemble to one array.

     ivar : int
       Index of the VAR file, if var_file is not specified.

     var_list : list of str 
       List of variables to read. If not provided or empty, 
       all available variables will be read by default. 

     quiet : bool
         Flag for switching off output.

     precision : string
         Float 'f', double 'd' or half 'half'.

     trimall : bool
         Trim the data cube to exclude ghost zones.

     swap_endian : bool
         Flag for swapping the endianness.


    Returns
    -------
    Varraw
        Class instance containing the varfile information.

    Examples
    --------
    Read the latest var.dat file and print the shape of the uu array:
    >>> var = pc.read.varraw()
    >>> print(var.uu.shape)

    Read the VAR2 file, only the logarithmic temperature & density, and remove the ghost zones:
    >>> var = pc.read.varraw(var_file='VAR2', var_list=['lnTT', 'lnrho'], trimall=True)
    >>> print(var.lnTT.shape)
    or  
    >>> print(var.f[:,:,:,var.tags['lnTT']].shape)  

    Read the VAR2 file, only the vector potential & z-velocity, in single precision:
    >>> var = pc.read.varraw(ivar=2, var_list=['aa', 'uz'], precision='f')
    >>> ay = var.ay[:]
    or  
    >>> ay = var.aa[:,:,:,1]
    or
    >>> ay = var.f[:,:,:,var.tags['ay']]

    Read the VAR2 file, only the vector potential & z-velocity, in single precision:
    >>> var = pc.read.varraw(ivar=2, var_list=['ax', 'uz','ay','az'], precision='f')
    >>> ay = var.ay[:]
    or  
    >>> ay = var.aa[:,:,:,1]
    or
    >>> ay = var.f[:,:,:,var.tags['ay']]     
    """

    varraw_tmp = Varraw()
    varraw_tmp.read(*args, **kwargs)
    return varraw_tmp


class Varraw(object):
    """
    Varraw -- holds Pencil Code VAR file data.
    """

    def __init__(self):
        """
        Fill members with default values.
        """

        self.t = 0.0
        self.dx = self.dy = self.dz = 0.0
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

    def keys(self):
        for i in self.__dict__.keys():
            print(i)

    def read(self, var_file="", datadir="data", proc=-1, ivar=-1, var_list=None, quiet=False, precision="d", trimall=False, swap_endian=False):
        """
        read(var_file="", datadir="data", proc=-1, ivar=-1, var_list=None, quiet=False, 
        precision="d", trimall=False, swap_endian=False)

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
            Processor to be read. If -1, read all and assemble to one array.

        ivar : int
            Index of the VAR file, if var_file is not specified.

        var_list : list of str 
            List of variables to read. If not provided or empty, 
            all available variables will be read by default. 

        quiet : bool
            Flag for switching off output.

        precision : string
            Float 'f', double 'd' or half 'half'.

        trimall : bool
            Trim the data cube to exclude ghost zones.

        swap_endian : bool
            Flag for swapping the endianness.


        Returns
        -------
        Varraw
            Class instance containing the varfile information.

        Examples
        --------
        Read the latest var.dat file and print the shape of the uu array:
        >>> var = pc.read.varraw()
        >>> print(var.uu.shape)

        Read the VAR2 file, only the temperature & density, and remove the ghost zones:
        >>> var = pc.read.varraw(var_file='VAR2', var_list=['lnTT', 'lnrho'], trimall=True)
        >>> print(var.lnTT.shape)
        or  
        >>> print(var.f[:,:,:,var.tags['lnTT']].shape)  

        Read the VAR2 file, only the vector potential & z-velocity, in single precision:
        >>> var = pc.read.varraw(ivar=2, var_list=['aa', 'uz'], precision='f')
        >>> ay = var.ay[:]
        or  
        >>> ay = var.aa[:,:,:,1]
        or
        >>> ay = var.f[:,:,:,var.tags['ay']]

        Read the VAR2 file, only the vector potential & z-velocity, in single precision:
        >>> var = pc.read.varraw(ivar=2, var_list=['ax', 'uz','ay','az'], precision='f')
        >>> ay = var.ay[:]
        or  
        >>> ay = var.aa[:,:,:,1]
        or
        >>> ay = var.f[:,:,:,var.tags['ay']]         
        """

        import numpy as np
        import os
        from scipy.io import FortranFile
        from pencil import read
        from collections import defaultdict

        if var_list is None:
            var_list = []

        if precision == "f":
            dtype = np.float32
        elif precision == "d":
            dtype = np.float64
        elif precision == "half":
            dtype = np.float16

        datadir = os.path.expanduser(datadir)

        index = read.index(datadir=datadir)

        param = read.param(datadir=datadir, quiet=True, conflicts_quiet=True)

        if param.io_strategy == "HDF5":
            raise NotImplementedError(
                "IO strategy {} not supported by this Python module. Use read.var() instead.".format(
                    param.io_strategy
                )
            )
        else:

            allprocs = 0

            # Name of varfile to read.
            if not var_file:
                if ivar < 0:
                    var_file = "var.dat"
                else:
                    var_file = "VAR" + str(ivar)

            if os.path.isfile(os.path.join(datadir, 'proc0', var_file)) and os.path.isdir(os.path.join(datadir, 'proc1')) and not os.path.isfile(os.path.join(datadir, 'proc1', var_file)):
                allprocs = 2
            if os.path.isfile(os.path.join(datadir, 'allprocs', var_file)):
                allprocs = 1

            if allprocs != 0 and proc >= 0:
                raise ValueError(
                    " 'proc' cannot be specified for collectively written files. Set it to -1 (default value)")

            # Set f77 value according to allprocs.
            f77 = 0 if allprocs == 1 else 1

            # Global dimensions.
            dim = read.dim(datadir, proc=-1)

            # We know from param whether we have to read 2-D or 3-D data.
            run2D = param.lwrite_2d

            # coord_system = param.coord_system

            # grid = read.grid(datadir=datadir, proc=proc, precision=precision, quiet=True, trim=trimall)

            # Read local dimensions.
            nprocs = dim.nprocx * dim.nprocy * dim.nprocz
            ipx_start = ipy_start = ipz_start = 0
            if allprocs == 1:
                procdim = dim
                ipx_end = ipy_end = ipz_end = 0
            else:
                ipz_end = dim.nprocz - 1
                if allprocs == 2:
                    procdim = read.dim(datadir, proc=0)
                    ipx_end = ipy_end = 0
                    procdim.nx = dim.nxgrid
                    procdim.ny = dim.nygrid
                    procdim.mx = dim.mxgrid
                    procdim.my = dim.mygrid
                    procdim.mw = procdim.mx * procdim.my * procdim.mz
                else:
                    if proc < 0:
                        procdim = read.dim(datadir, proc=0)
                        ipx_end = dim.nprocx-1
                        ipy_end = dim.nprocy-1
                    else:
                        procdim = read.dim(datadir, proc)
                        ipx_start = procdim.ipx
                        ipy_start = procdim.ipy
                        ipz_start = procdim.ipz
                        ipx_end = ipx_start
                        ipy_end = ipy_start
                        ipz_end = ipz_start

            mvar_io = dim.mvar
            if param.lwrite_aux:
                mvar_io += dim.maux

            if dim.precision == "D":
                read_dtype = np.float64
                data_bytes = 8
            else:
                if precision == "d":
                    raise ValueError(
                        " Reading a single precision run in double precision doesn't make sense! Please set precision to 'f' or 'half' ")
                read_dtype = np.float32
                data_bytes = 4

            if swap_endian:  # not sure whether this would really work for swapping endianness: need to test!
                read_dtype = np.dtype(read_dtype).newbyteorder()

            # Define the global arrays.
            x = np.zeros(dim.mx, dtype=precision)
            y = np.zeros(dim.my, dtype=precision)
            z = np.zeros(dim.mz, dtype=precision)

            # Read information about the file's contents and set up variable lists.

            varcontent = dict(
                sorted(index.__dict__.items(), key=lambda item: item[1]))

            # Filter duplicate indices read in from index.pro, if any.
            # In case of duplicate indices, keep only the variable name ending with 'x'.
            # If the duplicate indices have multiple variable names ending with 'x', keep the last one ending with 'x'.
            # If the duplicate indices have no variable name ending with 'x', keep the last one.
            # The criteria might need more refining.
            index_to_varname = defaultdict(list)

            for k, v in varcontent.items():
                index_to_varname[v].append(k)

            filtered_varcontent = {}

            for v, keys in index_to_varname.items():
                if len(keys) == 1:
                    filtered_varcontent[keys[0]] = v
                else:
                    varnames_with_x = [k for k in keys if k.endswith('x')]
                    if varnames_with_x:
                        keptname = varnames_with_x[-1]
                        print(
                            f"Warning: From index.pro, names {keys} share index {v}. Keeping name: '{keptname}' for this index.")
                    else:
                        keptname = keys[-1]
                        print(
                            f"Warning: From index.pro, names {keys} share index {v}. Keeping name: '{keptname}' for this index.")
                    filtered_varcontent[keptname] = v

            varcontent = filtered_varcontent

            # Safety check before proceeding.
            if len(varcontent) != mvar_io:
                raise ValueError(
                    "Inconsistency: no. of variables read in from index.pro doesn't match total no. of variables from dimensions file")

            content = ', '.join(varcontent.keys())

            # Expand known vector variables (if present). Can also accomodate quantities with different number of components,i.e. not necessarily 3.
            # (Provided they are added with appropriate component names here. So please add as required. Only the most common ones are here.)
            # Please check that the names are correct.
            # It is still possible to read the grouped variables not present here by specifying all the components separately in var_list
            # (or leaving var_list empty/unspecified, so that all variables are read).
            vector_groups = {
                'uu': ['ux', 'uy', 'uz'],
                'aa': ['ax', 'ay', 'az'],
                'qq': ['qx', 'qy', 'qz'],
                'bb': ['bx', 'by', 'bz'],
                'jj': ['jx', 'jy', 'jz'],
                'ee': ['ex', 'ey', 'ez'],
                'uu_sph': ['uu_sphx', 'uu_sphy', 'uu_sphz'],
                'bb_sph': ['bb_sphx', 'bb_sphy', 'bb_sphz'],
                'adv_der_uu': ['adv_der_ux', 'adv_der_uy', 'adv_der_uz']
            }

            if var_list:
                expanded = []
                for key in var_list:
                    if key in vector_groups:
                        expanded.extend(vector_groups[key])
                    else:
                        expanded.append(key)
                var_list = expanded

            # default to all variables if var_list is empty
            var_list = var_list or list(varcontent.keys())

            # validate that all requested variables exist
            missing_vars = set(var_list) - varcontent.keys()
            if missing_vars:
                raise KeyError(
                    f"The following variables were not found in varfile content: {', '.join(missing_vars)}")

            indices = [(key, varcontent[key]) for key in var_list]
            num_read = len(var_list)
            read_content = ', '.join(var_list)

            proc_mx = procdim.mx
            proc_my = procdim.my
            proc_mz = procdim.mz

            if run2D:
                if dim.nxgrid == 1:
                    proc_mx = 1
                if dim.nygrid == 1:
                    proc_my = 1
                if dim.nzgrid == 1:
                    proc_mz = 1

            # Display information about the file's contents.

            if not quiet:
                print()
                print(f"The file {var_file} contains: {content}")
                if len(var_list) < len(varcontent):
                    print(f"Will read only: {read_content}")
                print()
                print(f"The grid dimension is {dim.mx}, {dim.my}, {dim.mz}")
                print()

            # Initialise f-array and set marker values.
            self.f = np.zeros(
                (dim.mx, dim.my, dim.mz, num_read), dtype=precision)

            markers = 0 if f77 == 0 else 1

            # Iterate over processors.
            t = None

            for ipz in range(ipz_start, ipz_end + 1):
                for ipy in range(ipy_start, ipy_end + 1):
                    for ipx in range(ipx_start, ipx_end + 1):

                        iproc = ipx + ipy * dim.nprocx + ipz * dim.nprocx * dim.nprocy

                        x_off = (ipx - ipx_start) * procdim.nx
                        y_off = (ipy - ipy_start) * procdim.ny
                        z_off = (ipz - ipz_start) * procdim.nz

                        # Setup the coordinate mappings from the processor to the full domain.
                        # (Don't overwrite ghost zones of the lower processor.)
                        x_add_glob = dim.nghostx * \
                            (ipx != ipx_start or proc_mx == 1)
                        y_add_glob = dim.nghosty * \
                            (ipy != ipy_start or proc_my == 1)
                        z_add_glob = dim.nghostz * \
                            (ipz != ipz_start or proc_mz == 1)

                        x_add_proc = 0 if proc_mx == 1 else x_add_glob
                        y_add_proc = 0 if proc_my == 1 else y_add_glob
                        z_add_proc = 0 if proc_mz == 1 else z_add_glob

                        x_end = x_off + proc_mx - 1 + x_add_glob - x_add_proc
                        y_end = y_off + proc_my - 1 + y_add_glob - y_add_proc
                        z_end = z_off + proc_mz - 1 + z_add_glob - z_add_proc

                        # Build the full path and filename.
                        if allprocs == 1:
                            procdir = "allprocs"
                        else:
                            procdir = f"proc{iproc}"
                            if allprocs == 0 and not quiet:
                                print(f"Loading chunk {iproc + 1} of {nprocs}")

                        filename = os.path.join(datadir, procdir, var_file)

                        # Check for existence.
                        if not os.path.isfile(filename):
                            raise FileNotFoundError(
                                f'ERROR: File not found "{filename}"')

                        # Open a varfile and read some data!
                        mxyz = np.int64(proc_mx) * \
                            np.int64(proc_my) * np.int64(proc_mz)

                        with open(filename, "rb") as fil:
                            f_loc = np.zeros(
                                (proc_mx, proc_my, proc_mz, num_read), dtype=precision)
                            for pos in range(num_read):
                                pa = (indices[pos][1]) - 1
                                offset = data_bytes * pa * \
                                    mxyz + np.int64(markers * 4)
                                fil.seek(offset)
                                read_buffer = np.fromfile(
                                    fil, dtype=read_dtype, count=mxyz)
                                read_buffer = read_buffer.reshape(
                                    (proc_mx, proc_my, proc_mz), order='F')
                                f_loc[:, :, :, pos] = read_buffer.astype(dtype)
                                self.f[
                                    x_off + x_add_glob: x_end + 1,
                                    y_off + y_add_glob: y_end + 1,
                                    z_off + z_add_glob: z_end + 1,
                                    pos
                                ] = (read_buffer[
                                    x_add_proc:,
                                    y_add_proc:,
                                    z_add_proc:
                                ]).astype(dtype)

                            x_end = x_off + procdim.mx - 1
                            y_end = y_off + procdim.my - 1
                            z_end = z_off + procdim.mz - 1

                            additional = data_bytes * mvar_io * \
                                mxyz + np.int64(2 * markers * 4)
                            # skip past direct-access data, move to the sequentially written part
                            fil.seek(additional)

                            # Use FortranFile to correctly handle sequential record markers.
                            f_fortran = FortranFile(fil, "r")

                            t_test = dtype(0.0)
                            # Read the sequentially written data towards the end.
                            if allprocs == 1:  # collectively written files
                                record = f_fortran.read_record(read_dtype, np.dtype((read_dtype, (dim.mx))), np.dtype(
                                    (read_dtype, (dim.my))), np.dtype((read_dtype, (dim.mz))), read_dtype, read_dtype, read_dtype)
                                t_test = dtype(record[0].item())
                                x = (record[1]).astype(dtype)
                                y = (record[2]).astype(dtype)
                                z = (record[3]).astype(dtype)
                                self.dx = dtype(record[4].item())
                                self.dy = dtype(record[5].item())
                                self.dz = dtype(record[6].item())
                            elif allprocs == 2:  # xy-collectively written files for each ipz-layer
                                record = f_fortran.read_record(read_dtype)
                                t_test = dtype(record[0].item())
                                if iproc == 0:
                                    record = f_fortran.read_record(np.dtype((read_dtype, (dim.mx))), np.dtype(
                                        (read_dtype, (dim.my))), np.dtype((read_dtype, (dim.mz))), read_dtype, read_dtype, read_dtype)
                                    x = (record[0]).astype(dtype)
                                    y = (record[1]).astype(dtype)
                                    z = (record[2]).astype(dtype)
                                    self.dx = dtype(record[3].item())
                                    self.dy = dtype(record[4].item())
                                    self.dz = dtype(record[5].item())
                            else:  # distributed files
                                if param.lshear:
                                    record = f_fortran.read_record(read_dtype, np.dtype((read_dtype, (procdim.mx))), np.dtype(
                                        (read_dtype, (procdim.my))), np.dtype((read_dtype, (procdim.mz))), read_dtype, read_dtype, read_dtype, read_dtype)
                                    t_test = dtype(record[0].item())
                                    x_loc = (record[1]).astype(dtype)
                                    y_loc = (record[2]).astype(dtype)
                                    z_loc = (record[3]).astype(dtype)
                                    self.dx = dtype(record[4].item())
                                    self.dy = dtype(record[5].item())
                                    self.dz = dtype(record[6].item())
                                    self.deltay = dtype(record[7].item())
                                else:
                                    record = f_fortran.read_record(read_dtype, np.dtype((read_dtype, (procdim.mx))), np.dtype(
                                        (read_dtype, (procdim.my))), np.dtype((read_dtype, (procdim.mz))), read_dtype, read_dtype, read_dtype)
                                    t_test = dtype(record[0].item())
                                    x_loc = (record[1]).astype(dtype)
                                    y_loc = (record[2]).astype(dtype)
                                    z_loc = (record[3]).astype(dtype)
                                    self.dx = dtype(record[4].item())
                                    self.dy = dtype(record[5].item())
                                    self.dz = dtype(record[6].item())

                                x[x_off:x_end + 1] = x_loc
                                y[y_off:y_end + 1] = y_loc
                                z[z_off:z_end + 1] = z_loc

                            if t is None:
                                t = t_test
                            if t != t_test:
                                print(
                                    f"ERROR: TIMESTAMP IS INCONSISTENT: {filename}")
                                print(f"t != t_test: {t} != {t_test}")
                                raise ValueError(
                                    "Timestamp mismatch. Please check consistency.")

            if not quiet:
                print(f" t = {np.format_float_positional(t)}\n")

            self.t = t

            if proc < 0:
                self.x = x
                self.y = y
                self.z = z
                l1 = dim.l1
                l2 = dim.l2
                m1 = dim.m1
                m2 = dim.m2
                n1 = dim.n1
                n2 = dim.n2
            else:
                self.x = x_loc
                self.y = y_loc
                self.z = z_loc
                self.f = f_loc
                l1 = procdim.l1
                l2 = procdim.l2
                m1 = procdim.m1
                m2 = procdim.m2
                n1 = procdim.n1
                n2 = procdim.n2
                if run2D:
                    if dim.nxgrid == 1:
                        self.x = x_loc[l1: l2 + 1]
                        l1 = l2 = 0
                    if dim.nygrid == 1:
                        self.y = y_loc[m1: m2 + 1]
                        m1 = m2 = 0
                    if dim.nzgrid == 1:
                        self.z = z_loc[n1: n2 + 1]
                        n1 = n2 = 0

            # Remove ghost zones if requested.
            if trimall:
                self.x = self.x[l1: l2 + 1]
                self.y = self.y[m1: m2 + 1]
                self.z = self.z[n1: n2 + 1]

                sliced_f = self.f[l1: l2 + 1, m1: m2 + 1, n1: n2 + 1, :]
                # leave out the last dimension before squeezing f
                squeeze_axes = tuple(i for i in range(
                    3) if sliced_f.shape[i] == 1)
                self.f = np.squeeze(sliced_f, axis=squeeze_axes)

            else:
                self.l1 = l1
                self.l2 = l2  # deliberately different from varfile.py where there is an extra +1 at the end
                self.m1 = m1
                self.m2 = m2
                self.n1 = n1
                self.n2 = n2

            tags = {name: pos for pos, (name, _) in enumerate(indices)}

            # Assign an attribute to self for each variable defined in
            # var_list so that e.g. self.ux is the x-velocity
            for key, pos in tags.items():
                setattr(self, key, self.f[..., pos])

            # Special treatment for vector quantities (if all components present in f).
            # Can also accomodate quantities with different number of components,i.e. not necessarily 3.

            for vec_name, components in vector_groups.items():
                if all(comp in tags for comp in components):
                    pos_in_f = [tags[comp] for comp in components]
                    if pos_in_f == list(range(pos_in_f[0], pos_in_f[0] + len(pos_in_f))):
                        setattr(self, vec_name,
                                self.f[..., pos_in_f[0]: pos_in_f[-1] + 1])
                    else:
                        # gather and concatenate explicitly
                        parts = [self.f[..., p][..., np.newaxis]
                                 for p in pos_in_f]
                        setattr(self, vec_name, np.concatenate(parts, axis=-1))

            self.tags = tags
