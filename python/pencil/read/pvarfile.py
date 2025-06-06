import numpy as np
import os
from scipy.io import FortranFile
from os.path import expanduser, isdir, join
from pencil import read
from pencil.math import is_number

def pvar(*args, **kwargs):
    """
    pvar(pvarfile='', datadir='data', proc=-1, ipvar=-1, quiet=True,
        ID=False, pflist=None, sim=None, precision='f', dtype=np.float64)

    Read PVAR files from Pencil Code. If proc < 0, then load all data
    and assemble, otherwise load VAR file from specified processor.

    The file format written by output() (and used, e.g. in pvar.dat)
    consists of the followinig Fortran records:
    1. [npar]
    2. indices(npar)
    3. pdata(npvar, npar)
    Here npvar denotes the number of slots, i.e. 1 for one scalar field, 3 for
    one vector field, 6 for pvar.dat in the case of npar particles with 3
    coordinates and 3 velocity components.

    Parameters
    ----------
     pvarfile : string
         Name of the VAR file.
         If not specified, use var.dat (which is the latest snapshot of the fields)

     datadir : string
         Directory where the data is stored.

     proc : int
         Processor to be read. If -1 read all and assemble to one array.

     ipvar : int
       Index of the VAR file, if var_file is not specified.

     quiet : bool
         Flag for switching off output.

     ID : bool
         Flag for including the particle IDs in the object.

     pflist : bool
         If present list of exclusive basic pfarrays to include

     sim : pencil code simulation object
         Contains information about the local simulation.

     precision : string
         Float 'f', double 'd' or half 'half'.

     lpersist : bool
         Read the persistent variables if they exist

    Returns
    -------
    DataCube
        Instance of the pencil.read.var.DataCube class.
        All of the computed fields are imported as class members.

    Examples
    --------
    Read the latest var.dat file and print the shape of the uu array:
    >>> pvar = pc.read.pvar()
    >>> print(pvar.px.shape)

    Read the PVAR2 file, and include only the x coordinates and velocity
    e.g., for instance to reduce memory load for large arrays.
    >>> pvar = pc.read.pvar(pvar_file='PVAR2', pflist=['px','pvx'])
    >>> print(pvar.pvx.shape)
    """

    from pencil.sim import __Simulation__

    started = None

    for a in args:
        if isinstance(a, __Simulation__):
            started = a.started()
            break

    if "sim" in kwargs.keys():
        # started = kwargs['sim'].started()

        started = True
    elif "datadir" in kwargs.keys():
        from os.path import join, exists

        if exists(join(kwargs["datadir"], "time_series.dat")):
            started = True
    else:
        from os.path import join, exists

        if exists(join("data", "time_series.dat")):
            started = True

    if not started:
        if "ipvar" in kwargs:
            if kwargs["ipvar"] != 0:
                print("ERROR: Simulation has not yet started. There are no pvar files.")
                return False

    pvar_tmp = ParticleData()
    pvar_tmp.read(*args, **kwargs)
    return pvar_tmp

class ParticleData(object):
    """
    ParticleData -- holds Pencil Code PVAR file data.
    """


    def __init__(self):
        """
        Fill members with default values.
        """

    def keys(self):
        for i in self.__dict__.keys():
            print(i)

    def read(
        self,
        pvarfile="",
        datadir="data",
        proc=-1,
        proclist=None,
        ipvar=-1,
        quiet=True,
        pflist=None,
        ID=False,
        sim=None,
        precision="f",
        dtype=np.float64,
    ):
        """
        pvar(pvar_file='', datadir='data', proc=-1, ipvar=-1, quiet=True,
            pflist=None, sim=None, precision='f', dtype=np.float64)

        Read PVAR files from Pencil Code. If proc < 0, then load all data
        and assemble, otherwise load VAR file from specified processor.

        The file format written by output() (and used, e.g. in pvar.dat)
        consists of the followinig Fortran records:
        1. [npar]
        2. indices(npar)
        3. pdata(npvar, npar)
        Here npvar denotes the number of slots, i.e. 1 for one scalar field, 3 for
        one vector field, 6 for pvar.dat in the case of npar particles with 3
        coordinates and 3 velocity components.

        Parameters
        ----------
         pvarfile : string
             Name of the VAR file.
             If not specified, use var.dat (which is the latest snapshot of the fields)

         datadir : string
             Directory where the data is stored.

         proc : int
             Processor to be read. If -1 read all and assemble to one array.

         ipvar : int
           Index of the VAR file, if var_file is not specified.

         quiet : bool
             Flag for switching off output.

         ID : bool
             Flag for including the particle IDs in the object.

         pflist : bool
             If present list of exclusive basic pfarrays to include

         sim : pencil code simulation object
             Contains information about the local simulation.

         precision : string
             Float 'f', double 'd' or half 'half'.

         lpersist : bool
             Read the persistent variables if they exist

        Returns
        -------
        DataCube
            Instance of the pencil.read.var.DataCube class.
            All of the computed fields are imported as class members.

        Examples
        --------
        Read the latest var.dat file and print the shape of the uu array:
        >>> pvar = pc.read.pvar()
        >>> print(pvar.px.shape)

        Read the PVAR2 file, and include only the x coordinates and velocity
        e.g., for instance to reduce memory load for large arrays.
        >>> pvar = pc.read.pvar(pvar_file='PVAR2', pflist=['px','pvx'])
        >>> print(pvar.pvx.shape)
        """

        if sim is None:
            datadir = expanduser(datadir)
            dim = read.dim(datadir, proc=proc)
            pdim = read.pdim(datadir)
            param = read.param(datadir=datadir, quiet=quiet, conflicts_quiet=True)
            pindex = read.index(datadir=datadir,filename="particle_index.pro")
        else:
            datadir = expanduser(sim.datadir)
            dim = read.dim(datadir, proc=proc)
            pdim = read.pdim(datadir)
            param = read.param(datadir=sim.datadir, quiet=True, conflicts_quiet=True)
            pindex = read.index(datadir=datadir,filename="particle_index.pro")

        #constrain selection of particle variables
        pfkeys = dict()
        if pflist:
            if not isinstance(pflist, list):
                pflist = list(pflist)
        else:
            pflist=list(pindex.__dict__.keys())
        for key in pflist:
            if key in pindex.__dict__.keys():
                pfkeys[key]=pindex.__getattribute__(key)
        npvar = len(pfkeys)
        if ID:
            pfkeys["ID"] = -1
        if not quiet:
            print(pfkeys.keys())

        #cleanup of varfile string
        if not ipvar<0:
            pvarfile=ipvar
        if is_number(pvarfile):
            pvarfile = "PVAR" + str(pvarfile)
        pvarfile = str(pvarfile)
        if pvarfile == "var.dat":
            pvarfile = "pvar.dat"
        if pvarfile[:3] == "VAR":
            pvarfile = "P" + varfile
        if len(pvarfile)==0:
            pvarfile = "pvar.dat"

        if param.io_strategy == "HDF5":
            import h5py
            pvarfile = str.strip(pvarfile, ".dat") + ".h5"
            with h5py.File(join(datadir, "allprocs", pvarfile), "r") as hf:
                for key in hf["part"].keys():
                    if key in pfkeys.keys():
                        setattr(self, key.lower(), hf["part"][key][()])
        else:
            self._read_pvar_nonhdf5(dim, pdim, param, proclist, proc, datadir, npvar, dtype, ID, pvarfile, pfkeys)

        try:
            #Position vector
            setattr(self, "xxp", np.array([self.xp, self.yp, self.zp]))
        except:
            pass
        try:
            #Velocity vector
            setattr(self, "vvp", np.array([self.vpx,self.vpy,self.vpz]))
        except:
            pass
        try:
            #Magnetic vector
            setattr(self, "bbp", np.array([self.bpx,self.bpy,self.bpz]))
        except:
            pass
        try:
            #Spin vector
            setattr(self, "ssp", np.array([self.psx,self.psy,self.psz]))
        except:
            pass
        try:
            #Rate if strain tensor
            setattr(self, "wij", np.array([self.up11,self.up12,self.up13],
                                          [self.up21,self.up22,self.up23],
                                          [self.up31,self.up32,self.up33]))
        except:
            pass

    def __natural_sort(self, procs_list):
        """
        Sort array in a more natural way, e.g. 9VAR < 10VAR
        """

        import re

        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", key)]
        return sorted(procs_list, key=alphanum_key)

    def _read_singleproc_dat(self, file_name, output_dtype, read_precision, mpvar, mpaux):
        """
        Note that at least with io_dist, all the particle variables (including
        aux) are written out into the snapshots (see
        io_dist.f90/output_part_snap).

        Kishore: I am not sure if the auxiliary variables should be kept in the
        Kishore: snapshot; Fred's earlier implementation seemed to assume the
        Kishore: auxiliary variables are not in the snapshot.
        """
        with FortranFile(file_name) as infile:
            ind1 = infile.read_record(dtype='i')[0]
            ids = infile.read_record(dtype='i')
            data = output_dtype(infile.read_record(dtype=read_precision))
            data = data.reshape((mpvar+mpaux,ind1))
            return ids, data, ind1

    def _read_pvar_nonhdf5(self, dim, pdim, param, proclist, proc, datadir, npvar, dtype, ID, pvarfile, pfkeys):
        if dim.precision == "D":
            read_precision = "d"
        else:
            read_precision = "f"

        if isinstance(proclist, list):
            ind1 = 0
            proc_dirs = list()
            for idir in proclist:
                if isdir(join(datadir, "proc"+str(idir))):
                    proc_dirs.append("proc" + str(idir))
                    file_name = join(datadir,"proc"+str(idir), pvarfile)
                    infile = FortranFile(file_name)
                    ind1 += infile.read_record(dtype='i')[0]
                    infile.close()
                else:
                    raise ValueError(f"{idir} is not a valid proc directory")
            npar = ind1
        elif proc < 0:
            proc_dirs = self.__natural_sort(
                filter(lambda s: s.startswith("proc"), os.listdir(datadir))
            )
            if proc_dirs.count("proc_bounds.dat") > 0:
                proc_dirs.remove("proc_bounds.dat")
            if param.lcollective_io:
                proc_dirs = ["allprocs"]
            npar = pdim.npar

        if proc_dirs is not None:
            ptmp=np.zeros((npvar,npar), dtype=dtype)
            if ID:
                idtmp=np.zeros((npar), dtype=dtype)

            ind0 = 0
            for directory in proc_dirs:
                file_name = join(datadir, directory, pvarfile)
                ids, data, ind1 = self._read_singleproc_dat(file_name, dtype, read_precision, pdim.mpvar, pdim.mpaux)
                if ID:
                    idtmp[ind0:ind0+ind1] = ids
                for idx, key in zip(range(npvar),pfkeys.keys()):
                    ptmp[idx, ind0:ind0+ind1] = data[pfkeys[key]-1]
                ind0 += ind1
        else:
            file_name = join(datadir, "proc" + str(proc), pvarfile)
            ids, data, _ = self._read_singleproc_dat(file_name, dtype, read_precision, pdim.mpvar, pdim.mpaux)
            if ID:
                idtmp = ids
            for idx, key in zip(range(npvar),pfkeys.keys()):
                ptmp[idx] = data[pfkeys[key]-1]

        for idx, key in zip(range(npvar),pfkeys.keys()):
            if "ID" in key:
                setattr(self, key.lower(), idtmp)
            else:
                setattr(self, key.lower(), ptmp[idx])
