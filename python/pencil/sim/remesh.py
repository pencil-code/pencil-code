"""
Remesh mature simulation snapshot with [nx,ny,nz] dimensions onto new
simulation with new grid dimensions and optionally alternate cpu layout
copying the base simulation files, existing start output files.

uses:
local_remesh to apply the interpolation onto a variable array
get_dstgrid to derive the new grid layout
src2dst_remesh to create the new simulation object and files
"""

from fileinput import input
import numpy as np
from scipy.interpolate import interp1d
import subprocess as sub
import sys
from sys import stdout
import os
from pencil.math.derivatives import grad
from pencil.io import open_h5, group_h5, dataset_h5
from os.path import exists
try:
    import f90nml
    lnml=True
except:
    lnml=False

def local_remesh(var, xsrc, ysrc, zsrc, xdst, ydst, zdst, quiet=True, kind="linear"):
    """
    local_remesh(var, xsrc, ysrc, zsrc, xdst, ydst, zdst, quiet=True, kind="linear")

    Parameters
    ----------
    var : np.array
        Snapshot scalar numpy array of shape [mz, my, mx].

    xsrc, ysrc, zsrc : ndarrays
        Grid x, y, z arrays from source simulation.

    xdst, ydst, zdst : ndarrays
      Grid x, y, z arrays for destination simulation.

    kind : string
      interpolation method

    quiet : bool
      Flag for switching of output.

    Useage: interpolate a 3D data array of arbitrary shape onto an equivalent
            grid of alternate shape
    """

    tmp = var.copy()
    if not quiet:
        print("x", tmp.shape, xsrc.shape, xdst.shape)
        print("x", tmp.shape, xsrc.min(), xsrc.max(), xdst.min(), xdst.max())
    if not xsrc.size == xdst.size:
        interp = interp1d(xsrc, tmp, axis=-1, kind=kind, fill_value="extrapolate")
        tmp = interp(xdst)
    if not quiet:
        print("y", tmp.shape, ysrc.shape, ydst.shape)
        print("y", tmp.shape, ysrc.min(), ysrc.max(), ydst.min(), ydst.max())
    if not ysrc.size == ydst.size:
        interp = interp1d(ysrc, tmp, axis=-2, kind=kind, fill_value="extrapolate")
        tmp = interp(ydst)
    if not quiet:
        print("z", tmp.shape, zsrc.shape, zdst.shape)
        print("z", tmp.shape, zsrc.min(), zsrc.max(), zdst.min(), zdst.max())
    if not zsrc.size == zdst.size:
        interp = interp1d(zsrc, tmp, axis=-3, kind=kind, fill_value="extrapolate")
        tmp = interp(zdst)

    return tmp


def get_dstgrid(
    srch5,
    srcsim,
    dsth5,
    dstsim,
    ncpus=[1, 1, 1],
    nxyz=None,
    multxyz=[2, 2, 2],
    fracxyz=[1, 1, 1],
    srcghost=3,
    dstghost=3,
    dtype=np.float64,
    lsymmetric=True,
    quiet=True,
    dstprecision=[b"D"],
    srcchunks=None,
):
    """
    get_dstgrid(srch5, srcsim, dsth5, dstsim, ncpus=[1,1,1], multxyz=[2,2,2],
               nxyz=None, fracxyz=[1,1,1], srcghost=3, dstghost=3, dtype=np.float64,
               lsymmetric=True, quiet=True, dstprecision=[b"D"],
               srchunks=srchunks)

    Parameters
    ----------
    srch5 : obj
        hdf5 object from source simulation.

    srcsim : simulation object
        src simulation object

    dsth5 : obj
        hdf5 object for destination simulation data.

    dstsim : simulation object
        dst simulation object

    ncpus : int
        Array of nprocx, nprocy, and nprocz to apply for new simulation.

    nxyz : bool
        integer list of lenght 3 with new size of grid excluding ghosts.

    multxyz : list
        Factors by which to multiply old sim dimensions yxz order.

    fracxyz : list
        Factors by which to divide old sim dimensions yxz order.

    srcghost : int
      Number of ghost zones from the source order of accuracy (mx-nx)/2

    dstghost : int
        Number of ghost zones for the destination order of accuracy (mx-nx)/2

    dtype : 'string'
      Precision used in destination simulation. Default double.

    lsymmetric : bool
        Option to make non-periodic grid symmetric about old sim centre.
        Otherwise the lower boundary is retained from old sim grid.

    quiet : bool
        Flag for switching of output.

    dstprecision :
        floating point precision of new simulation grid.

    srcchunks : bool
        list of index limits [l[0]:l[1],m[0]:m[1],n[0]:n[1]] for subdomain remesh

    """
    from os.path import join, abspath
    import h5py
    # TBA
    # check prime factorization of the result and display for proc options
    # if using fft check options for grid and cpu layout
    # handling non-equidistant grids tba

    # copy settings from srcsim and revise with changes to dstsim var.h5
    if dsth5.__contains__("settings"):
        dsth5.__delitem__("settings")
    if srch5:
        srcsets = srch5["settings"]
        srch5.copy("settings",dsth5)
        sets = dsth5["settings"]
    #else:
    #    sets = dsth5.require_group("settings")
    #    srcsets = dict()
    #    for key in srcsim.dim.__dict__.keys():
    #        htype=type(srcsim.dim.__getattribute__(key))
    #        if htype=="str":
    #            srcsets[key]=np.array([srcsim.dim.__getattribute__(key).encode("utf-8")])
    #        else:
    #            srcsets[key]=np.array([srcsim.dim.__getattribute__(key)])
    #        sets.create_dataset(key, data=srcsets[key][()])
    lvalid = True
    if srcchunks:
        lvalid == (len(srcchunks) == 3)
        for lls, lld in zip(srcchunks,[["l1","l2"],["m1","m2"],["n1","n2"]]):
            lvalid ==(len(lls) == len(lld))
            lvalid == (lls[0] >= srcsets[lld[0]][0])
            lvalid == (lls[1] <= srcsets[lld[1]][0])
    if not lvalid:
        print("srcchunks {} must be list of three pairs of\n".format(srcchunks)+
              "indices that fit within the source grid shape.\n"+
               "Ignoring srcchunks and remeshing full grid.")
        srcchunks=None
    # update grid dimensions
    if not nxyz:
        sets["nx"][()] = int(srcsets["nx"][0] * multxyz[0] / fracxyz[0])
        sets["mx"][()] = sets["nx"][()] + 2 * dstghost
        sets["ny"][()] = int(srcsets["ny"][0] * multxyz[1] / fracxyz[1])
        sets["my"][()] = sets["ny"][()] + 2 * dstghost
        sets["nz"][()] = int(srcsets["nz"][0] * multxyz[2] / fracxyz[2])
        sets["mz"][()] = sets["nz"][()] + 2 * dstghost
        sets["l1"][()] = dstghost
        sets["l2"][()] = sets["mx"][()] - 1 - dstghost
        sets["m1"][()] = dstghost
        sets["m2"][()] = sets["my"][()] - 1 - dstghost
        sets["n1"][()] = dstghost
        sets["n2"][()] = sets["mz"][()] - 1 - dstghost
    else:
        sets["nx"][()] = nxyz[0]
        sets["mx"][()] = sets["nx"][()] + 2 * dstghost
        sets["ny"][()] = nxyz[1]
        sets["my"][()] = sets["ny"][()] + 2 * dstghost
        sets["nz"][()] = nxyz[2]
        sets["mz"][()] = sets["nz"][()] + 2 * dstghost
        sets["l1"][()] = dstghost
        sets["l2"][()] = sets["mx"][()] - 1 - dstghost
        sets["m1"][()] = dstghost
        sets["m2"][()] = sets["my"][()] - 1 - dstghost
        sets["n1"][()] = dstghost
        sets["n2"][()] = sets["mz"][()] - 1 - dstghost
    if not ncpus == [1, 1, 1]:
        sets["nprocx"][()] = ncpus[0]
        sets["nprocy"][()] = ncpus[1]
        sets["nprocz"][()] = ncpus[2]

    # copy the grid from the srcsim to dstsim var.h5 and grid.h5
    if dsth5.__contains__("grid"):
        dsth5.__delitem__("grid")
    if srch5:
        srcgrid = srch5["grid"]
        srch5.copy("grid",dsth5)
        grid=dsth5["grid"]
    #else:
    #    grid=dsth5.require_group("grid")
    #    srcgrid = dict()
    #    for key in srcsim.ghost_grid.__dict__.keys():
    #        if not key == "t":
    #            srcgrid[key]=srcsim.ghost_grid.__getattribute__(key)
    #            grid.create_dataset(key, data=srcgrid[key])
    #    for key, i in zip(["Ox", "Oy", "Oz"],[0,1,2]):
    #        srcgrid[key]=srcsim.param["xyz0"][i]
    #        grid.create_dataset(key, data=srcgrid[key])
    # replace grid data changed for dstsim
    if srcchunks:
        for x, sch in zip(["x","y","z"],srcchunks):
            grid["O"+x][()] = srcgrid[x][sch[0]+srcghost]
            grid["L"+x][()] = srcgrid[x][sch[1]+srcghost] - grid["O"+x][()]
            x0,x1 = grid["O"+x][()], grid["O"+x][()] + grid["L"+x][()]
            dx = dtype(grid["L"+x][()]/sets["n"+x][0])
            grid["d"+x][()] = dx
            grid.__delitem__(x)
            grid.create_dataset(x, data=np.linspace(x0-dstghost*dx,
                                                    x1+(dstghost-1)*dx,
                                                    sets["m"+x][0]), dtype=dtype)
            grid.__delitem__("d"+x+"_1")
            grid.create_dataset("d"+x+"_1", (sets["m"+x][0],), dtype=dtype)
            grid["d"+x+"_1"][:] = dtype(1.0 / np.gradient(grid[x][()]))
            grid.__delitem__("d"+x+"_tilde")
            grid.create_dataset("d"+x+"_tilde", data=np.zeros(sets["m"+x][0]), dtype=dtype)
            grid["d"+x+"_tilde"][:] = dtype(1.0 / np.gradient(grid["d"+x+"_1"][()]))
    else:
        for x, ii in zip(["x","y","z"],[0,1,2]):
            if not srcsim.param["lequidist"][ii]:
                if rank == 0:
                    print(
                        "get_dstgrid WARNING: non-equidistant grid not implemented\n",
                        "continuing with equidistant grid.\n",
                        "Please implement non-equidistant grid options.",
                    )
            if not sets["m"+x][()] == srcsets["m"+x][()]:
                dx = dtype(grid["L"+x][()]/sets["n"+x][0])
                grid["d"+x][()] = dx
                if not srcsim.param["lperi"][ii]:
                    grid["O"+x][()] = srcgrid["O"+x][()] - 0.5*dx
                x0,x1 = grid["O"+x][()], grid["O"+x][()] + grid["L"+x][()]
                print("dx {}, xyz0 {}, Lxyz {}".format(dx, x0, x1-x0))
                grid.__delitem__(x)
                if srcsim.param["lperi"][ii]:
                    grid.create_dataset(x, data=np.linspace(x0-dstghost*dx,
                                                    x1+(dstghost-1)*dx,
                                                    sets["m"+x][0]), dtype=dtype)
                else:
                    grid.create_dataset(x, data=np.linspace(x0-dstghost*dx,
                                                    x1+(dstghost-1)*dx,
                                                    sets["m"+x][0]), dtype=dtype)
                if srcsim.param["lshift_origin"][ii] or lsymmetric:
                    grid[x][()] += dtype(0.5 * grid["d"+x][()])
                elif srcsim.param["lshift_origin_lower"][ii]:
                    grid[x][()] -= dtype(0.5 * grid["d"+x][()])
                grid.__delitem__("d" + x + "_1")
                grid.create_dataset(
                    "d" + x + "_1", shape = grid[x][()].shape, dtype=dtype
                )
                grid["d" + x + "_1"][()] = dtype(1.0 / np.gradient(grid[x][()]))
                grid.__delitem__("d" + x + "_tilde")
                grid.create_dataset(
                    "d" + x + "_tilde",
                    shape=grid["d" + x + "_1"][()].shape,
                    dtype=dtype,
                )
                grid["d" + x + "_tilde"][()] = dtype(np.gradient(grid["d" + x + "_1"][()]))
    if dsth5.__contains__("unit"):
        dsth5.__delitem__("unit")
    if srch5:
        srch5.copy("unit", dsth5)
    #else:
    #    dsth5.create_group("unit")
    #    for key in srcsim.param.keys():
    #        if "unit_" in key and not "_unit_" in key:
    #            dkey=key.split("_")[-1]
    #            if dkey == "system":
    #                dsth5["unit"].create_dataset(dkey, data=srcsim.param[key].encode("utf-8"))
    #            else:
    #                dsth5["unit"].create_dataset(dkey, data=srcsim.param[key])
    gridh5 = h5py.File(join(dstsim.datadir, "grid.h5"), 'w')
    dsth5.copy("settings", gridh5)
    dsth5.copy("grid", gridh5)
    dsth5.copy("unit", gridh5)
    gridh5.close()
        #else:
        #    var=pc.read.var(proc=rank, lpersist=True)
        #    pers=dsth5.require_group("persist")
        #    for key in var.persist.keys():
        #        tmp = np.zeros(nprocs)
        #        tmp[:] = var.persist.__getattribute__(key)
        #        htype=type(var.persist.__getattribute__(key))
        #        if type(srch5["persist"][key][0].item())==float:
        #            htype=dtype
        #        if htype=="str":
        #            htype=type(var.persist.__getattribute__(key).encode("utf-8"))
        #        try:
        #            pers.require_dataset(key, (nprocs,), dtype=htype)
        #        except:
        #            pers.__delitem__(key)
        #            pers.require_dataset(key, (nprocs,), dtype=htype)
        #        pers[key][()]=tmp

def src2dst_remesh(
    src=None,
    dst=None,
    h5in="var.h5",
    h5out="var.h5",
    nxyz=None,
    multxyz=[2, 2, 2],
    fracxyz=[1, 1, 1],
    srcchunks=None,
    srcghost=3,
    dstghost=3,
    srcdatadir="data/allprocs",
    dstdatadir="data/allprocs",
    dstprecision=[b"D"],
    lsymmetric=True,
    quiet=True,
    kind="linear",
    check_grid=True,
    optionals=True,
    nmin=32,
    rename_submit_script=False,
    MBmin=64.0,
    ncpus=[1, 1, 1],
    start_optionals=False,
    hostfile=None,
    submit_new=False,
    chunksize=5000.0,
    lfs=False,
    MB=32,
    count=1,
    size=1,
    rank=0,
    comm=None,
    farray=None,
    index_farray=None,
    datasets='all',
    remesh=True,
    newtime=None,
):
    """
    src2dst_remesh(src=None, dst=None, h5in='var.h5', h5out='var.h5', nxyz=None,
                   multxyz=[2, 2, 2], fracxyz=[1, 1, 1], srcchunks=None, srcghost=3,
                   dstghost=3, srcdatadir='data/allprocs', dstdatadir='data/allprocs',
                   dstprecision=[b'D'], lsymmetric=True, quiet=True, kind="linear",
                   check_grid=True, OVERWRITE=False, optionals=True, nmin=32,
                   rename_submit_script=False, MBmin=5.0, ncpus=[1, 1, 1],
                   start_optionals=False, hostfile=None, submit_new=False,
                   chunksize=1000.0, lfs=False,  MB=1, count=1, size=1, rank=0,
                   comm=None, farray=None, index_farray=None, datasets='all', remesh=True)

    Parameters
    ----------
    src : string
        Source relative or absolute path to source simulation.

    dst : string
        Destination relative or absolute path to destination simulation.

    h5in : string
        Source simulation data file to be copied and remeshed.

    h5out : string
        Destination simulation file to be written.

    nxyz : list
        If not None a list of 3 integers [nx, ny, nz] defining the size of
        the destination domain

    multxyz : list
        Factors by which to multiply old sim dimensions xyz order.

    fracxyz : list
        Factors by which to divide old sim dimensions xyz order.

    srcghost : int
        Number of ghost zones from the source order of accuracy (mx-nx)/2.

    dstghost : int
        Number of ghost zones for the destination order of accuracy (mx-nx)/2.

    srcdatadir : string
        Path from source simulation directory to data.

    dstdatadir :
        Path from destination simulation directory to data.

    dstprecision : string
        Floating point precision settings [b'S'] or [b'D'].

    lsymmetric : bool
        Option to make non-periodic grid symmetric about old sim centre.
        Otherwise the lower boundary is retained from old sim grid.

    quiet : bool
        Flag for switching of output.

    check_grid : bool
        Flag to run check on grid and cpu layout before executing remesh.

    OVERWRITE : bool
        Flag to overwrite existing simulation directory and filesin dst.

    optionals : bool
        Copy simulation files with True or specify list of names (string) for
        additional files from src sim directory.

    nmin : int
        Minimum length along coordinate after splitting by proc.

    rename_submit_script : bool
        Edit lines in submission files vcopied from src to dst.
        Not yet operational.

    MBmin : float
        Minimum size in MB of data on a sinlge proc pf ncpus total processes.

    ncpus : ndarray
        Array of nprocx, nprocy, and nprocz to apply for new simulation.

    start_optionals : bool
        Copy simulation files output by start.x with True or specify list of
        names (string) for additional files from src sim data directory.

    hostfile : string
        Specify name of host config file argument in pc_build.
        Not yet operational.

    submit_new : bool
        Execute changes to submission files, compile and run simulation.
        Not yet operational.

    chunksize : float
      Size in megabytes of snapshot variable before chunked remesh is used.

    lfs : bool
      Flag to set the striping for large file sizes to improve IO efficiency.

    MB : float
      Size of data to write contiguously before moving to new OST on lustre.

    count : int
        Number of OSTs across which the data will be shared for IO operations.

    size : int
        Number of MPI processes

    rank : int
        ID of processor

    comm :
        MPI library calls

    farray : string list
        Select subset of the farray from src.

    index_farray : integer list
        set of index at which interpolation shall start

    datasets : string
        "all" default to interpolate full farray, or specify selection

    remesh :  bool
        remesh or just copy
    """

    import h5py
    import os
    from os.path import join, abspath
    import time

    from pencil import read
    from pencil.io import mkdir
    from pencil.sim import simulation
    from pencil.math import cpu_optimal
    from pencil import is_sim_dir

    start_time = time.time()

    #------------------------------------------------------------------
    #identify source and destination paths and copy simulation files as
    #required to destination path

    if rank == 0 or rank == size - 1:
        print("started at {}".format(time.ctime(start_time)))
    # set dtype from precision
    if dstprecision[0] == b"D":
        dtype = np.float64
        bytesize = 8
    elif dstprecision[0] == b"S":
        dtype = np.float32
        bytesize = 4
    else:
        if rank == 0 or rank == size - 1:
            print("precision " + dstprecision + " not valid")
        return 1

    if not src or src==".":
        src=os.getcwd()
    if not dst or dst==".":
        dst=src+"copy"
    if is_sim_dir(src):
        srcsim = simulation(src, quiet=quiet)
    else:
        if rank == 0 or rank == size - 1:
            print('src2dst_remesh ERROR: src"' + src + '" is not a valid simulation path')
        return 1
    print("dst is sim",is_sim_dir(dst))
    lsim=False
    if rank == 0:
        if is_sim_dir(dst):
            dstsim = simulation(dst, quiet=quiet)
            dstname = str.split(dst, "/")[-1]
            dstpath = str.strip(dst, dstname)
            srcsim.copy(
                    name=dstname,
                    quiet=quiet,
                    OVERWRITE=False,
                    optionals=optionals,
                    start_optionals=start_optionals,
                    rename_submit_script=rename_submit_script,
                )
            if not os.path.isfile(join(dst,"data/allprocs")):
                mkdir(join(dst,"data/allprocs"),lfs=lfs,MB=MB,count=count)
                mode = "w"
            if os.path.isfile(join(dst,"data/allprocs",h5out)):
                lsim=True
                mode = "r+"
        else:
            mode = "a"
            print("setting up simulation")
            if rank == 0:
                dstname = str.split(dst, "/")[-1]
                dstpath = str.strip(dst, dstname)
                if len(dstpath) == 0:
                    dstpath = str.strip(srcsim.path, srcsim.name)
                dstsim = srcsim.copy(
                    path_root=dstpath,
                    name=dstname,
                    quiet=quiet,
                    OVERWRITE=True,
                    optionals=optionals,
                    start_optionals=start_optionals,
                    rename_submit_script=rename_submit_script,
                )
                mkdir(join(dst,"data","allprocs"),lfs=lfs,MB=MB,count=count)
        print("dst is sim already?",lsim)
    if rank == 0:
        if not lsim:
            # If no h5 started determine grid and settings for dst
            srch5 = h5py.File(join(srcsim.path, srcdatadir, h5in),"r")
            print("opening {} file on rank {}".format(join(srcsim.path, srcdatadir, h5in),rank))
            with h5py.File(join(dstsim.path, dstdatadir, h5out),mode) as dsth5:
                get_dstgrid(
                    srch5,
                    srcsim,
                    dsth5,
                    dstsim,
                    ncpus=ncpus,
                    nxyz=nxyz,
                    multxyz=multxyz,
                    fracxyz=fracxyz,
                    srcghost=srcghost,
                    dstghost=dstghost,
                    dtype=dtype,
                    lsymmetric=lsymmetric,
                    quiet=quiet,
                    dstprecision=dstprecision,
                    srcchunks=srcchunks,
                    )

            print("get_dstgrid completed on rank {}".format(rank))
            with h5py.File(join(dstsim.path, dstdatadir, h5out),"r+") as dsth5:
                # use settings to determine available proc dist then set ncpus
                factors = cpu_optimal(
                    dsth5["settings/nx"][0],
                    dsth5["settings/ny"][0],
                    dsth5["settings/nz"][0],
                    mvar=dsth5["settings/mvar"][0],
                    maux=dsth5["settings/maux"][0],
                    par=srcsim.param,
                    nmin=nmin,
                    MBmin=MBmin,
                    remesh=remesh
                )
                print(
                    "remesh check grid: optional cpus up to min grid of "
                    + "nmin={}\n".format(nmin)
                    + "cpu options {}\n".format(factors)
                    + "new mesh: {}, {}, {}\n".format(
                        dsth5["settings/nx"][0],
                        dsth5["settings/ny"][0],
                        dsth5["settings/nz"][0],
                    )
                    + 'To execute remesh set "check_grid=False".'
                )
                if ncpus == [1, 1, 1]:
                    ncpus = [factors[1][0], factors[1][1], factors[1][2]]
                dsth5["settings/nprocx"][0] = ncpus[0]
                dsth5["settings/nprocy"][0] = ncpus[1]
                dsth5["settings/nprocz"][0] = ncpus[2]
                dsth5["settings/precision"][0] = dstprecision
                nprocs = ncpus[0] * ncpus[1] * ncpus[2]
                srcprocs = (
                      srcsim.dim.nprocx
                    * srcsim.dim.nprocy
                    * srcsim.dim.nprocz
                )
                if srcprocs > nprocs:
                    print(
                         "\n**********************************************************\n"
                         + "remesh WARNING: {} procs reduced from {}.\n".format(
                             nprocs, srcprocs
                         )
                         + "Review multxyz {} and fracxyz {} for more\n".format(
                             multxyz, fracxyz
                         )
                         + "efficient parallel processing options."
                         + "\n**********************************************************\n"
                        )
                    if check_grid:
                        return 1
                # add persistent variables of correct size and set time
                if "persist" in srch5.keys():
                    pers=dsth5.require_group("persist")
                    for key in srch5["persist"].keys():
                        htype=type(srch5["persist"][key][0])
                        if type(srch5["persist"][key][0].item())==float:
                            htype=dtype
                        tmp = np.empty(nprocs,dtype=htype)
                        tmp[:] = srch5["persist"][key][0]
                        try:
                            pers.create_dataset(key, data=tmp)
                        except:
                            pers.__delitem__(key)
                            pers.create_dataset(key, data=tmp)
                dsth5.require_dataset("time", (), dtype=dtype)
                if newtime:
                    dsth5["time"][()] = newtime
                else:
                    dsth5["time"][()] = srch5["time"][()]
                # update correct values to param.nml files if possible
                if lnml:
                    nmldata=f90nml.read(join(srcsim.datadir,"param.nml"))
                    nmldata["io_pars"]["io_strategy"]="HDF5"
                    nmldata["init_pars"]["xyz0"]=[dsth5["grid/Ox"][()],
                                                  dsth5["grid/Oy"][()],
                                                  dsth5["grid/Oz"][()]]
                    nmldata["init_pars"]["lxyz"]=[dsth5["grid/Lx"][()],
                                                  dsth5["grid/Ly"][()],
                                                  dsth5["grid/Lz"][()]]
                    nmldata["init_pars"]["xyz1"]=[dsth5["grid/Ox"][()] + dsth5["grid/Lx"][()],
                                                  dsth5["grid/Oy"][()] + dsth5["grid/Ly"][()],
                                                  dsth5["grid/Oz"][()] + dsth5["grid/Lz"][()]]
                    nmldata.write(join(dstsim.datadir,"param.nml"),"f")
                    nmldata=f90nml.read(join(srcsim.datadir,"param2.nml"))
                    nmldata["io_pars"]["io_strategy"]="HDF5"
                    nmldata["run_pars"]["xyz0"]=[dsth5["grid/Ox"][()],
                                                  dsth5["grid/Oy"][()],
                                                  dsth5["grid/Oz"][()]]
                    nmldata["run_pars"]["lxyz"]=[dsth5["grid/Lx"][()],
                                                  dsth5["grid/Ly"][()],
                                                  dsth5["grid/Lz"][()]]
                    nmldata["run_pars"]["xyz1"]=[dsth5["grid/Ox"][()] + dsth5["grid/Lx"][()],
                                                  dsth5["grid/Oy"][()] + dsth5["grid/Ly"][()],
                                                  dsth5["grid/Oz"][()] + dsth5["grid/Lz"][()]]
                    nmldata.write(join(dstsim.datadir,"param2.nml"),"f")
                    dstsim.update(hard=True)
                else:
                    dstsim.update()
                    for key in ["nprocx","nprocy","nprocz","nxgrid","nygrid","nzgrid"]:
                        dstsim.dim.__setattr__(key,dsth5[join("settings",key.strip("grid"))][0])
                dstsim.dim.__setattr__("ncpus",nprocs)
                #The subsequest tools need to be improved to complete revision of *.local and
                #compilation if required -- see pipelines
                # update correct values to src/cparam.local
                cpar = open(join(dstsim.path,"src/cparam.local"),"a")
                for key in ["ncpus","nprocx","nprocy","nprocz","nxgrid","nygrid","nzgrid"]:
                    cpar.write("integer, parameter :: {}={}\n".format(key,dstsim.dim.__getattribute__(key)))
                cpar.close()
            srch5.close()

    #-----------------------------------------------------------------
    #update destination simulation object on all ranks

    if comm:
        comm.Barrier()
    if not rank == 0:
        dstsim = simulation(dst, quiet=quiet)
        if not lnml:
            with h5py.File(join(dstsim.path, dstdatadir, h5out),"r+") as dsth5:
                for key in ["nprocx","nprocy","nprocz","nxgrid","nygrid","nzgrid"]:
                    dstsim.dim.__setattr__(key,dsth5[join("settings",key.strip("grid"))][0])
    if srcchunks:
        nxin = srcchunks[0][1]-srcchunks[0][0]
        nyin = srcchunks[1][1]-srcchunks[1][0]
        nzin = srcchunks[2][1]-srcchunks[2][0]
        l1in,l2in=srcchunks[0][0], srcchunks[0][1]+2*srcghost
        m1in,m2in=srcchunks[1][0], srcchunks[1][1]+2*srcghost
        n1in,n2in=srcchunks[2][0], srcchunks[2][1]+2*srcghost
        xin, yin, zin = (
          srcsim.ghost_grid.x[l1in:l2in],
          srcsim.ghost_grid.y[m1in:m2in],
          srcsim.ghost_grid.z[n1in:n2in]
          )
        print("start and end values {},{}".format(l1in,l2in))
    else:
        nxin, nyin, nzin = srcsim.dim.nxgrid, srcsim.dim.nygrid, srcsim.dim.nzgrid
        l1in,l2in=srcsim.dim.l1-srcghost, srcsim.dim.l2+srcghost+1
        m1in,m2in=srcsim.dim.m1-srcghost, srcsim.dim.m2+srcghost+1
        n1in,n2in=srcsim.dim.n1-srcghost, srcsim.dim.n2+srcghost+1
        xin, yin, zin = (
              srcsim.ghost_grid.x[()],
              srcsim.ghost_grid.y[()],
              srcsim.ghost_grid.z[()]
              )
    if rank==0:
        print("nxin, nyin, nzin", xin.size, yin.size, zin.size)

    #Determine which fields from the farray need to be remeshed
    if not farray:
        fields = list(srcsim.index.__dict__.keys())
        if not index_farray:
            index_fields = list()
            for ind in range(len(fields)):
                index_fields.append(0)
    else:
        fields = list()
        index_fields = list()
        if isinstance(farray, list):
            if not index_farray:
                index_farray = list()
                for ind in range(len(farray)):
                    index_farray.append(0)
            if isinstance(index_farray, list):
                index_farray = index_farray
            else:
                index_farray = list(index_farray)
            if not isinstance(index_farray[0], int):
                if rank == 0:
                    print("index_farray {} should be a list of start indices or else None to start writing at index 0.".format(index_farray))
                index_farray = [0]
                if len(index_farray) < len(farray):
                    for ind in range(len(farray)-len(index_farray)):
                        index_farray.append(0)
            for field, ind in zip(farray,index_farray):
                if field in srch5["data"].keys():
                    fields.append(field)
                    index_fields.append(ind)
                else:
                    if rank == 0:
                        print("{} in farray list is not in {} data".format(field, srch5.filename))
        else:
            if isinstance(farray, str):
                if farray in srch5["data"].keys():
                    fields.append(farray)
                    if not index_fields:
                        index_fields = [0]
                    elif isinstance(index_farray, list):
                        if not isinstance(index_farray[0], int):
                            if rank == 0:
                                print("index_farray {} should be a list of start indices or else None to start writing at index 0.".format(index_farray))
                            index_fields = [0]
                    else:
                        if not isinstance(index_farray, int):
                            if rank == 0:
                                print("index_farray {} should be a list of start indices or else None to start writing at index 0.".format(index_farray))
                            index_fields = [0]
                else:
                    if rank == 0:
                        print("farray {} is not a string or is not in {} data".format(farray, srch5.filename))
        if len(fields) == 0:
            print("farray list empty no data to remesh")
            return 1
    if comm:
        driver = "mpio"
        dh5 = h5py.File(join(dstsim.path, dstdatadir, h5out),"r+", driver=driver, comm=comm)
    else:
        dh5 = h5py.File(join(dstsim.path, dstdatadir, h5out),"r+" )
    print("dh5 name", dh5.filename)

    #determine if chunking is required due to memory constraints
    with dh5 as dsth5:
        if comm:
            try:
                dsth5.atomic = True
            except:
                print("atomic not supported with driver {}".format(driver))
            print(dsth5.filename,"atomic status:",dsth5.atomic)
        nx, ny, nz =(
            dsth5["settings"]["nx"][0],
            dsth5["settings"]["ny"][0],
            dsth5["settings"]["nz"][0],
        )
        nprocs = (
              dsth5["settings/nprocx"][0]
            * dsth5["settings/nprocy"][0]
            * dsth5["settings/nprocz"][0]
        )
        ncpus = [dsth5["settings/nprocx"][0],dsth5["settings/nprocy"][0],dsth5["settings/nprocz"][0]]
        mx, my, mz = (
                dsth5["settings"]["mx"][0],
                dsth5["settings"]["my"][0],
                dsth5["settings"]["mz"][0],
            )
        xout, yout, zout = (
              dsth5["grid/x"][()],
              dsth5["grid/y"][()],
              dsth5["grid/z"][()]
              )
    dstchunksize = bytesize * nx * ny * nz / 1024 / 1024 / size
    if rank == 0:
        print("rank {}, nx = {} dstchunksize {} chunksize {}, ncpus {}".format(rank, nx, dstchunksize, chunksize, ncpus))
    lchunks = False
    if dstchunksize > chunksize or size > 1:
        lchunks = True
        if size > 1:
            tchunks = size
        else:
            tchunks = int(dstchunksize/chunksize+1)
        nallchunks = cpu_optimal(
                nx,
                ny,
                nz,
                MBmin=tchunks,
                nsize=tchunks,
                size=tchunks,
                remesh=remesh
                )[1]
        if rank == 0:
            print("rank,lchunks, tchunks, nallchunks",rank,lchunks, tchunks, nallchunks)
        allindx = np.array_split(np.arange(nxin) + srcghost, nallchunks[0])
        allindy = np.array_split(np.arange(nyin) + srcghost, nallchunks[1])
        allindz = np.array_split(np.arange(nzin) + srcghost, nallchunks[2])
        alloutx = np.array_split(np.arange(nx) + dstghost, nallchunks[0])
        allouty = np.array_split(np.arange(ny) + dstghost, nallchunks[1])
        alloutz = np.array_split(np.arange(nz) + dstghost, nallchunks[2])
        chunks = list()
        for izz in range(nallchunks[2]):
            for iyy in range(nallchunks[1]):
                for ixx in range(nallchunks[0]):
                    chunks.append([ixx,iyy,izz])
        #print("rank {} chunks {} len chunks {}".format(rank, chunks, len(chunks)))
        rankchunks = np.array_split(chunks, size)
        for irank in range(size):
            if rank == irank:
                chunks = rankchunks[irank]
                #print("rank {} chunks {} len chunks {}".format(rank, chunks, len(chunks)))
                nchunks = 0
                lfirstx,lfirsty,lfirstz = list(),list(),list()
                llastx,llasty,llastz = list(),list(),list()
                for ichunk in chunks:
                    nchunks += 1
                    if ichunk[0] == 0:
                        lfirstx.append("True")
                    else:
                        lfirstx.append("False")
                    if ichunk[0] == nallchunks[0] - 1:
                        llastx.append("True")
                    else:
                        llastx.append("False")
                    if ichunk[1] == 0:
                        lfirsty.append("True")
                    else:
                        lfirsty.append("False")
                    if ichunk[1] == nallchunks[1] - 1:
                        llasty.append("True")
                    else:
                        llasty.append("False")
                    if ichunk[2] == 0:
                        lfirstz.append("True")
                    else:
                        lfirstz.append("False")
                    if ichunk[2] == nallchunks[2] - 1:
                        llastz.append("True")
                    else:
                        llastz.append("False")

    else:
        lfirstx,lfirsty,lfirstz = list("True"),list("True"),list("True")
        llastx,llasty,llastz = list("True"),list("True"),list("True")
        nchunks = 1
        chunks = [0,0,0]
    print("rank {} len chunks {}, len lfirstx".format(rank, len(chunks), len(lfirstx)))

    #Process remeshing for each field specified from the farray
    for key, ind_start in zip(fields,index_fields):
        if exists(join(dst,'PYTHONSTOP')):
            if rank == 0:
                cmd = 'rm -f '+join(dst,'PYTHONSTOP')
                print("rank {} ind {} for key {} at time {}".format(rank, ind, key, time.ctime(time.time())))
                os.system(cmd)
            if comm:
                comm.Barrier()
            sys.exit()
        if comm:
            driver = "mpio"
            sh5 = h5py.File(join(srcsim.path, srcdatadir, h5in), "r", driver=driver, comm=comm)
            dh5 = h5py.File(join(dstsim.path, dstdatadir, h5out),"r+", driver=driver, comm=comm)
        else:
            driver = None
            sh5 = h5py.File(join(srcsim.path, srcdatadir, h5in), "r", driver=driver)
            dh5 = h5py.File(join(dstsim.path, dstdatadir, h5out),"r+", driver=driver)
        with sh5 as srch5:
            if comm:
                try:
                    srch5.atomic = True
                except:
                    print("atomic not supported with driver {}".format(driver))
            with dh5 as dsth5:
                if comm:
                    try:
                        dsth5.atomic = True
                    except:
                        print("atomic not supported with driver {}".format(driver))
                if not quiet:
                    if rank == 0 or rank == size - 1:
                        print("nx {}, ny {}, nz {}".format(nx, ny, nz))
                        print("mx {}, my {}, mz {}".format(mx, my, mz))
                dsth5.require_group("data")
                if rank == 0 or rank == size - 1:
                    print("remeshing " + key)
                if not lchunks:
                    #invar = srch5["data"][key][n1:n2,m1:m2,l1:l2]
                    invar = srch5["data"][key][n1in:n2in,m1in:m2in,l1in:l2in]
                    print(key, srch5["data"][key][n1in:n2in,m1in:m2in,l1in:l2in].shape)
                    print(key, zout.size,  yout.size, xout.size)
                    var = local_remesh(
                        invar,
                        xin,
                        yin,
                        zin,
                        xout,
                        yout,
                        zout,
                        quiet=quiet,
                        kind=kind,
                    )
                    if rank == 0 or rank == size - 1:
                        print("serial rank {} writing {} shape {}".format(rank,key,var.shape))
                    dsth5["data"].require_dataset(key, shape=(mz, my, mx), dtype=dtype)
                    dsth5["data"][key][()] = var
                else:
                    dsth5["data"].require_dataset(key, shape=(mz, my, mx), dtype=dtype)
                    print("xin {}, yin {}, zin {}, xout {}, yout {}, zout {}".format(xin.shape, yin.shape, zin.shape, xout.shape, yout.shape, zout.shape))
                    for [ix, iy, iz], ixyz, firstz, lastz, firsty, lasty, firstx, lastx in zip(chunks,range(nchunks),lfirstz,llastz,lfirsty,llasty,lfirstx,llastx):
                        print("ixyz {} ix {} iy {} iz {}, firstz {}, lastz {}, firsty {}, lasty {}, firstx {}, lastx {}".format(ixyz, ix, iy, iz, firstz, lastz, firsty, lasty, firstx, lastx))
                        if ixyz >= ind_start:
                            n1, n2 = alloutz[iz][0] - dstghost, alloutz[iz][-1] + dstghost
                            #print("zout[n1] {}, n1 {}, zin[:n1] {}".format(zout[n1], n1, zin[n1]))
                            try:
                                srcn1 = np.max(
                                np.where(zin < zout[n1])
                                )
                                print(iz,"srcn1 {} zin {:.3f} n1 {} zout {:.3f}".format(
                                       srcn1,   zin[srcn1],n1,   zout[n1]))
                            except:
                                srcn1 = 0
                            try:
                                srcn2 = np.min(
                                    np.where(zin > zout[n2])
                                )
                                print(iz,"srcn2 {} zin {:.3f} n2 {} zout {:.3f}".format(
                                       srcn2,   zin[srcn2],n2,   zout[n2]))
                            except:
                                srcn2 = zin.size - 1
                            n1out = n1 + dstghost
                            n2out = n2 - dstghost + 1
                            varn1 = dstghost
                            varn2 = -dstghost
                            if firstz == True:
                                n1out = 0
                                varn1 = 0
                            if lastz == True:
                                n2out = n2 + 1
                                varn2 = n2 + 1
                            if not quiet:
                                print(
                                    "n1 {}, n2 {}, srcn1 {}, srcn2 {}".format(
                                        n1, n2, srcn1, srcn2
                                    )
                                )
                            m1, m2 = allouty[iy][0] - dstghost, allouty[iy][-1] + dstghost
                            try:
                                srcm1 = np.max(
                                np.where(yin < yout[m1])
                                )
                                print(iy,"srcm1 {} yin {:.3f} m1 {} yout {:.3f}".format(
                                       srcm1,   yin[srcm1],m1,   yout[m1]))
                            except:
                                srcm1 = 0
                            try:
                                srcm2 = np.min(
                                    np.where(yin > yout[m2])
                                )
                                print(iy,"srcm2 {} yin {:.3f} m2 {} yout {:.3f}".format(
                                       srcm2,   yin[srcm2],m2,   yout[m2]))
                            except:
                                srcm2 = yin.size - 1
                            m1out = m1 + dstghost
                            m2out = m2 - dstghost + 1
                            varm1 = dstghost
                            varm2 = -dstghost
                            if firsty == True:
                                m1out = 0
                                varm1 = 0
                            if lasty == True:
                                m2out = m2 + 1
                                varm2 = m2 + 1
                            if not quiet:
                                print(
                                    "rank {} m1 {}, m2 {}, srcm1 {}, srcm2 {}".format( rank,
                                        m1, m2, srcm1, srcm2
                                    )
                                )
                            l1, l2 = alloutx[ix][0] - dstghost, alloutx[ix][-1] + dstghost
                            try:
                                srcl1 = np.max(
                                np.where(xin < xout[l1])
                                )
                                print(ix,"srcl1 {} xin {:.3f} l1 {} xout {:.3f}".format(
                                       srcl1,   xin[srcl1],l1,   xout[l1]))
                            except:
                                srcl1 = 0
                            try:
                                srcl2 = np.min(
                                    np.where(xin > xout[l2])
                                )
                                print(ix,"srcl2 {} xin {:.3f} l2 {} xout {:.3f}".format(
                                       srcl2,   xin[srcl2],l2,   xout[l2]))
                            except:
                                srcl2 = zin.size - 1
                            l1out = l1 + dstghost
                            l2out = l2 - dstghost + 1
                            varl1 = dstghost
                            varl2 = -dstghost
                            if firstx == True:
                                l1out = 0
                                varl1 = 0
                            if lastx == True:
                                l2out = l2 + 1
                                varl2 = l2 + 1
                            if not quiet:
                                print(
                                    "l1 {}, l2 {}, srcl1 {}, srcl2 {}".format(
                                        l1, l2, srcl1, srcl2
                                    )
                                )
                            if not quiet:
                                print(
                                    "remeshing "
                                    + key
                                    + " chunk {}".format([iz, iy, ix])
                                )
                                print("ixyz {} xin  {:.3f} {:.3f}, yin  {:.3f} {:.3f}, zin  {:.3f} {:.3f}".format(
                                       ixyz,
                                       xin[srcl1], xin[srcl2],
                                       yin[srcm1], yin[srcm2],
                                       zin[srcn1], zin[srcn2]))
                                print("ixyz {} xout {:.3f} {:.3f}, yout {:.3f} {:.3f}, zout {:.3f} {:.3f}".format(
                                       ixyz,
                                       xout[l1], xout[l2],
                                       yout[m1], yout[m2],
                                       zout[n1], zout[n2]))
                            if not quiet:
                                print("rank {} ind {} writing {} shape {}".format(rank,ind,key,[n2out-n1out, m2out-m1out, l2out-l1out]))
                            invar = srch5["data"][key][
                                    srcn1 : srcn2 + 1,
                                    srcm1 : srcm2 + 1,
                                    srcl1 : srcl2 + 1,
                                ]
                            var = local_remesh(
                                invar,
                                xin[srcl1 : srcl2 + 1],
                                yin[srcm1 : srcm2 + 1],
                                zin[srcn1 : srcn2 + 1],
                                xout[l1 : l2 + 1],
                                yout[m1 : m2 + 1],
                                zout[n1 : n2 + 1],
                                quiet=quiet,
                            )
                            print("invar {} var {}".format(invar.shape,var.shape))
                            if not quiet:
                                print(
                                    "writing "
                                    + key
                                    + " shape {} chunk {}".format(
                                        var.shape, [iz, iy, ix]
                                    )
                                )
                            globals()[str(ixyz)+"time"] = time.time()
                            outvar = dtype(
                                var[varn1:varn2, varm1:varm2, varl1:varl2]
                            )
                            dsth5["data"][key][n1out:n2out, m1out:m2out, l1out:l2out] = outvar
                            print("wrote rank {} nchunk {} for {} in {} seconds".format(rank,[ix,iy,iz],key,time.time()-globals()[str(ixyz)+"time"]))
                            del(var,globals()[str(ixyz)+"time"],iz, firstz, lastz, iy, firsty, lasty, ix, firstx, lastx, ixyz)

    if rank==0 and submit_new:
        cmd = 'source '+join(srcsim.path,'src','.moduleinfo')
        os.system(cmd)
        os.chdir(dstsim.path)
        cmd = 'pc_setupsrc; make cleann'
        os.system(cmd)
        cmd = 'pc_build'
        if hostfile: cmd = cmd + ' -f '+hostfile
        process = sub.Popen(cmd.split(),stdout=sub.PIPE)
        process = sub.Popen(cmd.split(),stdout=sub.PIPE)
        output, error = process.communicate()
        print(cmd,output,error)
        cmd = "rm -f "+join(dstsim.path, dstdatadir, h5out+"copy")
        os.system(cmd)
    if comm:
        comm.Barrier()
    end_time = time.time()
    if rank == 0 or rank == size - 1:
        print(
              "end at {} after {} seconds".format(
                 time.ctime(end_time), end_time - start_time)
        )


