
# pstalk2.py
#
# Read the stalker files.
# Returns fs.ipar, fs.xp, etc.
#
# python translation of Anders' pstalker IDL script.
# Called it "2" because there is another pstalk, which 
# uses the IDL-python bridge (although that one does not 
# seem to work).
#
# W. Lyra (wlyra@nmsu.edu)
# 27 nov 2025: 
#
import os
import glob
import struct
import numpy as np
import h5py

from pencil import read

class Struct(dict):
    """Simple dict-like container with attribute (dot) access."""
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as e:
            raise AttributeError(name) from e

    def __setattr__(self, name, value):
        self[name] = value

def pstalk2(
    datadir="data",
    it0=0,
    it1=-1,
    swap_endian=False,
    quiet=False,
    noutmax=None,
    single=False,
    nstalk=None,
):
    """
    Read Pencil Code particle stalker data (PSTALK) and return a Struct with
    dot-access fields:

        obj.t      : (nout,)
        obj.ipar   : (nstalk,)
        obj.<field>: (nstalk, nout) for each field in particles_stalker_header.dat

    Parameters
    ----------
    datadir : str
        Base data directory (e.g. "data"). ~ expansion is supported.
    it0 : int
        First time index to read (inclusive).
    it1 : int
        If >0, print a log line whenever it % it1 == 0.
    swap_endian : bool
        If True, swap endianness when reading the binary PSTALK files.
    quiet : bool
        If True, suppress most logging.
    noutmax : int or None
        Maximum number of times to read. None => use all available.
    single : bool
        If True, store data in float32; otherwise float64.
    nstalk : int or None
        Max number of stalked particles. Defaults to pdim.npar_stalk.

    Returns
    -------
    obj : Struct
        Container with fields t, ipar, and one field per stalk quantity.
    """
    datadir = os.path.expanduser(datadir)
    if noutmax is None:
        noutmax = -1

    # ------------------------------------------------------------------
    # Read dim, pdim, param from pencil.read
    # ------------------------------------------------------------------
    dim = read.dim(datadir)
    pdim = read.pdim(datadir)
    param = read.param(datadir=datadir, quiet=quiet, conflicts_quiet=True)

    if nstalk is None:
        nstalk = pdim.npar_stalk

    # ------------------------------------------------------------------
    # HDF5 branch: datadir/allprocs/PSTALK*.h5
    # ------------------------------------------------------------------
    h5_files = sorted(glob.glob(os.path.join(datadir, "allprocs", "PSTALK*.h5")))
    if h5_files:
        if not quiet:
            print(
                "pstalk2: WARNING: PSTALK HDF5 detected. For efficiency consider "
                "using higher-level Pencil readers where available."
            )

        num_files = len(h5_files)
        dtype = np.float32 if single else np.float64

        stalk_obj = Struct()
        quantities = None

        for pos, fname in enumerate(h5_files):
            with h5py.File(fname, "r") as f:
                t_val = np.array(f["time"])[()]

                if pos == 0:
                    distribution = np.array(f["proc"]["distribution"])
                    num_part = int(distribution.sum())

                    stalk_group = f["stalker"]
                    quantities = list(stalk_group.keys())

                    stalk_obj.t = np.full(num_files, np.nan, dtype=dtype)
                    stalk_obj.ipar = np.array(stalk_group["ID"], dtype=np.int64)

                    for q in quantities:
                        if q.upper() == "ID":
                            continue
                        stalk_obj[q] = np.full(
                            (num_part, num_files), np.nan, dtype=dtype
                        )

                for q in quantities:
                    if q.upper() == "ID":
                        continue
                    stalk_obj[q][:, pos] = np.array(
                        f["stalker"][q], dtype=dtype, copy=False
                    )

                stalk_obj.t[pos] = t_val

        return stalk_obj

    # ------------------------------------------------------------------
    # Legacy binary PSTALK (particles_stalker.dat)
    # ------------------------------------------------------------------

    if pdim.npar_stalk == 0:
        print(
            "pstalk2: pdim.npar_stalk is zero - set it in cparam.local and rerun."
        )
        return Struct()

    # Do we stalk sink particles?
    lstalk_sink_particles = int(
        getattr(param, "lstalk_sink_particles", 0)
    )

    # ------------------------------------------------------------------
    # Read tstalk.dat -> (tout, noutmaxfile)
    # ------------------------------------------------------------------
    tstalk_file = os.path.join(datadir, "tstalk.dat")
    with open(tstalk_file, "r") as f:
        parts = f.read().split()
        tout = float(parts[0])
        noutmaxfile = int(parts[1])

    if noutmax == -1 or noutmax > noutmaxfile:
        nout = noutmaxfile
    else:
        nout = noutmax

    # ------------------------------------------------------------------
    # Read header to know which fields are stored
    # ------------------------------------------------------------------
    header_file = os.path.join(datadir, "particles_stalker_header.dat")
    with open(header_file, "r") as f:
        header = f.readline().strip()

    fields = [s.strip() for s in header.split(",") if s.strip()]
    nfields = len(fields)

    if not quiet:
        print(f"Going to read the {nfields} fields:")
        print("  ", fields)
        print("at", nout, "times")

    # ------------------------------------------------------------------
    # Allocate arrays
    # ------------------------------------------------------------------
    dtype = np.float32 if single else np.float64

    t = np.zeros(nout, dtype=dtype)                         # (nout,)
    array = np.zeros((nfields, nstalk, nout), dtype=dtype)  # (nfields, nstalk, nout)

    if lstalk_sink_particles:
        ipar_stalk = np.full((nstalk, nout), -1, dtype=np.int64)
        npar_stalk_read = np.zeros(nout, dtype=np.int64)
    else:
        ipar_stalk = np.arange(nstalk, dtype=np.int64)

    # ------------------------------------------------------------------
    # Fortran unformatted record helpers
    # ------------------------------------------------------------------
    endian = ">" if swap_endian else "<"
    float_fmt = "f" if single else "d"
    float_size = 4 if single else 8

    def read_record(fh):
        """Read one Fortran unformatted sequential record and return payload bytes."""
        header_bytes = fh.read(4)
        if not header_bytes:
            return None
        (nbytes,) = struct.unpack(endian + "i", header_bytes)
        payload = fh.read(nbytes)
        # trailing record length
        fh.read(4)
        return payload

    # ------------------------------------------------------------------
    # Loop over processors and read particles_stalker.dat
    # ------------------------------------------------------------------
    nprocs = dim.nprocx * dim.nprocy * dim.nprocz

    for iproc in range(nprocs):
        if not quiet:
            print(f"Reading data from processor {iproc}/{nprocs - 1}")
            print(
                "-------- iproc ------ it --------- t ----------- npar ------- "
            )

        it = 0
        ntread = 0

        proc_dir = os.path.join(datadir, f"proc{iproc}")
        fname = os.path.join(proc_dir, "particles_stalker.dat")
        if not os.path.exists(fname):
            continue

        with open(fname, "rb") as fh:
            while ntread < nout:
                rec1 = read_record(fh)
                if rec1 is None:
                    break  # EOF

                # Record 1: t_loc (float) + npar_stalk_loc (int32)
                t_loc = struct.unpack(endian + float_fmt, rec1[:float_size])[0]
                npar_stalk_loc = struct.unpack(
                    endian + "i", rec1[float_size : float_size + 4]
                )[0]

                ipar_loc = None
                array_loc = None

                # If npar_stalk_loc > 0, we must read records 2 and 3
                if npar_stalk_loc > 0:
                    # Record 2: ipar_loc (npar_stalk_loc ints)
                    rec2 = read_record(fh)
                    if rec2 is None:
                        break
                    ipar_loc = np.frombuffer(
                        rec2, dtype=endian + "i4", count=npar_stalk_loc
                    ).astype(np.int64)

                    # Record 3: array_loc (nfields * npar_stalk_loc floats)
                    rec3 = read_record(fh)
                    if rec3 is None:
                        break
                    array_loc = np.frombuffer(
                        rec3,
                        dtype=dtype,
                        count=nfields * npar_stalk_loc,
                    ).reshape((nfields, npar_stalk_loc), order="F")

                # Store if we are at/after it0
                if it >= it0:
                    if (it1 != -1) and (it % it1 == 0) and (not quiet):
                        print(iproc, it, t_loc, npar_stalk_loc)

                    if npar_stalk_loc > 0:
                        if not lstalk_sink_particles:
                            idx = ipar_loc - 1  # 1-based -> 0-based
                            array[:, idx, it - it0] = array_loc
                        else:
                            start = npar_stalk_read[it - it0]
                            stop = start + npar_stalk_loc

                            if stop > nstalk:
                                nstalk2 = stop
                                array2 = array
                                ipar2 = ipar_stalk

                                array = np.zeros(
                                    (nfields, nstalk2, nout), dtype=dtype
                                )
                                ipar_stalk = np.full(
                                    (nstalk2, nout), -1, dtype=np.int64
                                )

                                array[:, :nstalk, :] = array2
                                ipar_stalk[:nstalk, :] = ipar2
                                nstalk = nstalk2

                            array[:, start:stop, it - it0] = array_loc
                            ipar_stalk[start:stop, it - it0] = ipar_loc
                            npar_stalk_read[it - it0] += npar_stalk_loc

                    t[it - it0] = t_loc
                    ntread += 1

                it += 1

    # ------------------------------------------------------------------
    # Sink-particle reordering / trimming
    # ------------------------------------------------------------------
    if lstalk_sink_particles:
        array2 = array
        ipar2 = ipar_stalk

        flat = ipar2.ravel()
        flat_sorted = np.sort(flat)
        # unique but keep sorted order
        uniq_vals, uniq_idx = np.unique(flat_sorted, return_index=True)
        ipar_unique = uniq_vals

        if ipar_unique.size > 1:
            # mimic IDL's skipping of first element
            ipar_unique = ipar_unique[1:]

            if ipar_unique.size > nstalk:
                array = np.zeros(
                    (nfields, ipar_unique.size, nout), dtype=dtype
                )

            for k, ip in enumerate(ipar_unique):
                for jt in range(1, nout):
                    kk = np.where(ipar2[:, jt] == ip)[0]
                    if kk.size > 1:
                        array[:, k, :] = 0.0
                        break
                    if kk.size == 1:
                        array[:, k, jt] = array2[:, kk[0], jt]
                    else:
                        array[:, k, jt] = 0.0

            ipar_stalk = ipar_unique
            array = array[:, : ipar_stalk.size, :]
        else:
            ipar_stalk = np.array([-1], dtype=np.int64)
            array = array[:, 0:1, :]

    # ------------------------------------------------------------------
    # Build result Struct with dot-access
    # ------------------------------------------------------------------
    obj = Struct()
    obj.t = t
    obj.ipar = ipar_stalk

    # reshape each field -> (nstalk, nout)
    for i, name in enumerate(fields):
        obj[name] = array[i, :, :]

    return obj
