# pstalk2.py
#
# Read the stalker files.
# Returns fs.ipar, fs.xp, etc. with dot-access.
#
# Python translation of Anders Johansen's pc_read_pstalk IDL script,
# adapted to use the Pencil Python I/O.
#
# W. Lyra (wlyra@nmsu.edu)
# 27 Nov 2025
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
    ipar=None,
):
    """
    Read Pencil Code particle stalker data (PSTALK) and return a Struct with
    dot-access fields:

        fs.t      : (nout,)
        fs.ipar   : (nstalk,)
        fs.<field>: (nstalk, nout) for each field in particles_stalker_header.dat

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
        Maximum number of output times to read. None => use all available.
        Semantics match IDL: min(noutmax, noutmaxfile), with -1 meaning "all".
    single : bool
        If True, store output data in float32; otherwise float64.
        (Binary precision is auto-detected and then cast.)
    nstalk : int or None
        Max number of stalked particles to keep. Defaults to pdim.npar_stalk.
        In non-sink mode, particles with index >= nstalk are ignored.
        In sink-particle mode, nstalk may be increased internally if needed.

    ipar : int or list of int or None
        If given, read only these particle indices (memory-efficient).

    Returns
    -------
    fs : Struct
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

    # -----------------------------
    # Normalize ipar
    # -----------------------------
    if ipar is not None:
        if isinstance(ipar, int):
            ipar = [ipar]
        ipar_set = set(ipar)
        nkeep = len(ipar_set)
        keep_map = {gidx: k for k, gidx in enumerate(sorted(ipar_set))}
    else:
        ipar_set = None
        nkeep = nstalk
        keep_map = None

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
        out_dtype = np.float32 if single else np.float64

        fs = Struct()
        quantities = None

        for pos, fname in enumerate(h5_files):
            with h5py.File(fname, "r") as f:
                t_val = np.array(f["time"])[()]

                if pos == 0:
                    distribution = np.array(f["proc"]["distribution"])
                    num_part = int(distribution.sum())

                    stalk_group = f["stalker"]
                    quantities = list(stalk_group.keys())

                    fs.t = np.full(num_files, np.nan, dtype=out_dtype)
                    fs.ipar = np.array(stalk_group["ID"], dtype=np.int64)

                    for q in quantities:
                        if q.upper() == "ID":
                            continue
                        fs[q] = np.full(
                            (num_part, num_files), np.nan, dtype=out_dtype
                        )

                for q in quantities:
                    if q.upper() == "ID":
                        continue
                    data = np.array(f["stalker"][q])
                    fs[q][:, pos] = data.astype(out_dtype, copy=False)

                fs.t[pos] = out_dtype(t_val)

        # Apply ipar selection if requested
        if ipar is not None:
            if isinstance(ipar, int):
                ipar = [ipar]
            ipar_set = set(ipar)
            mask = np.isin(fs.ipar, list(ipar_set))
            for q in quantities:
                if q.upper() == "ID":
                    continue
                fs[q] = fs[q][mask, :]
            fs.ipar = fs.ipar[mask]

        return fs

    # ------------------------------------------------------------------
    # Legacy binary PSTALK (particles_stalker.dat)
    # ------------------------------------------------------------------

    if pdim.npar_stalk == 0:
        print(
            "pstalk2: pdim.npar_stalk is zero - set it in cparam.local and rerun."
        )
        return Struct()

    if "lstalk_sink_particles" in param.__dict__:
        lstalk_sink_particles = param.lstalk_sink_particles
    else:
        lstalk_sink_particles = 0


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
    # Allocate arrays (output dtype controlled by `single`)
    # ------------------------------------------------------------------
    out_dtype = np.float32 if single else np.float64

    t = np.zeros(nout, dtype=out_dtype)                         # (nout,)

    # --------------------------------------------------------------
    # Allocate array ONLY for kept particles
    # --------------------------------------------------------------
    array = np.zeros((nfields, nkeep, nout), dtype=out_dtype)

    if lstalk_sink_particles:
        ipar_stalk = np.full((nkeep, nout), -1, dtype=np.int64)
        npar_stalk_read = np.zeros(nout, dtype=np.int64)
    else:
        if ipar_set is not None:
            ipar_stalk = np.array(sorted(ipar_set), dtype=np.int64)
        else:
            ipar_stalk = np.arange(nstalk, dtype=np.int64)

    # ------------------------------------------------------------------
    # Fortran unformatted record helpers (auto-detect float precision)
    # ------------------------------------------------------------------
    endian = ">" if swap_endian else "<"

    def read_record(fh):
        """Read one Fortran unformatted sequential record and return (payload, nbytes)."""
        header_bytes = fh.read(4)
        if not header_bytes:
            return None, 0
        (nbytes,) = struct.unpack(endian + "i", header_bytes)
        payload = fh.read(nbytes)
        fh.read(4)  # trailing record length
        return payload, nbytes

    # These will be auto-detected from the first rec1 we see:
    read_float_fmt = None      # 'f' or 'd'
    read_float_size = None     # 4 or 8
    binary_float_dtype = None  # np.float32 or np.float64

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
                rec1, nbytes1 = read_record(fh)
                if rec1 is None:
                    break  # EOF

                # Auto-detect float precision from record length rec1:
                # rec1 contains: t_loc (float or double) + npar_stalk_loc (int32)
                if binary_float_dtype is None:
                    if nbytes1 == 8:      # 4 (float32) + 4 (int32)
                        read_float_fmt = "f"
                        read_float_size = 4
                        binary_float_dtype = np.float32
                    elif nbytes1 == 12:   # 8 (float64) + 4 (int32)
                        read_float_fmt = "d"
                        read_float_size = 8
                        binary_float_dtype = np.float64
                    else:
                        raise RuntimeError(
                            f"pstalk2: cannot infer float size from record length {nbytes1}"
                        )

                # Record 1: t_loc (float) + npar_stalk_loc (int32)
                t_loc_raw = struct.unpack(
                    endian + read_float_fmt, rec1[:read_float_size]
                )[0]
                npar_stalk_loc = struct.unpack(
                    endian + "i", rec1[read_float_size : read_float_size + 4]
                )[0]

                ipar_loc = None
                array_loc = None

                # If npar_stalk_loc > 0, read records 2 and 3
                if npar_stalk_loc > 0:
                    # Record 2: ipar_loc (npar_stalk_loc ints)
                    rec2, _ = read_record(fh)
                    if rec2 is None:
                        break
                    ipar_loc = np.frombuffer(
                        rec2, dtype=endian + "i4", count=npar_stalk_loc
                    ).astype(np.int64)

                    # Record 3: array_loc (nfields * npar_stalk_loc floats)
                    rec3, _ = read_record(fh)
                    if rec3 is None:
                        break
                    array_loc_native = np.frombuffer(
                        rec3,
                        dtype=binary_float_dtype,
                        count=nfields * npar_stalk_loc,
                    ).reshape((nfields, npar_stalk_loc), order="F")

                    # Cast to output dtype if needed
                    if array_loc_native.dtype != out_dtype:
                        array_loc = array_loc_native.astype(out_dtype)
                    else:
                        array_loc = array_loc_native

                # Store if we are at/after it0
                if it >= it0:
                    out_it = it - it0
                    if out_it >= nout:
                        # we already filled requested outputs for this proc
                        break

                    if (it1 != -1) and (it % it1 == 0) and (not quiet):
                        print(iproc, it, t_loc_raw, npar_stalk_loc)

                    if npar_stalk_loc > 0:
                        if not lstalk_sink_particles:
                            # Convert global stalk indices (1-based → 0-based)
                            idx_global = ipar_loc - 1

                            if keep_map is None:
                                # --- FULL LOAD (no ipar filtering) ---
                                mask = (idx_global >= 0) & (idx_global < nstalk)
                                if np.any(mask):
                                    array[:, idx_global[mask], out_it] = array_loc[:, mask]

                            else:
                                mask = np.isin(idx_global, list(ipar_set))
                                if np.any(mask):
                                    kompact = [keep_map[g] for g in idx_global[mask]]
                                    array[:, kompact, out_it] = array_loc[:, mask]
                        else:
                            # Sink-particle mode (currently not used in IDL, but logic kept)
                            start = npar_stalk_read[out_it]
                            stop = start + npar_stalk_loc

                            if stop > nstalk:
                                nstalk2 = stop
                                array2 = array
                                ipar2 = ipar_stalk

                                array = np.zeros(
                                    (nfields, nstalk2, nout), dtype=out_dtype
                                )
                                ipar_stalk = np.full(
                                    (nstalk2, nout), -1, dtype=np.int64
                                )

                                array[:, :nstalk, :] = array2
                                ipar_stalk[:nstalk, :] = ipar2
                                nstalk = nstalk2

                            array[:, start:stop, out_it] = array_loc
                            ipar_stalk[start:stop, out_it] = ipar_loc
                            npar_stalk_read[out_it] += npar_stalk_loc

                    # Cast t_loc to output dtype when storing
                    t[out_it] = out_dtype(t_loc_raw)
                    ntread += 1

                it += 1

    # ------------------------------------------------------------------
    # Sink-particle reordering / trimming (not active if lstalk_sink_particles==0)
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
                    (nfields, ipar_unique.size, nout), dtype=out_dtype
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
    fs = Struct()
    fs.t = t
    fs.ipar = ipar_stalk

    # reshape each field -> (nstalk, nout)
    for i, name in enumerate(fields):
        fs[name] = array[i, :, :]

    # ------------------------------------------------------------------
    # Memory-efficient particle selection
    # ------------------------------------------------------------------
    if ipar_set is not None:
        mask = np.isin(fs.ipar, list(ipar_set))
        for q in fields:
            fs[q] = fs[q][mask, :]
        fs.ipar = fs.ipar[mask]

    return fs
