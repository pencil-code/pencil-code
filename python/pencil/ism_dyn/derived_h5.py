# derived_h5.py
""" Derive auxilliary data and other diagnostics from var.h5 file and
    save to new h5 file

    uses:
      compute "data" arrays of size [nz,ny,nx] as required
      store "time" of snapshot
"""
import numpy as np
from pencil.math import dot, dot2, helmholtz_fft, cpu_optimal
from pencil.math.derivatives import curl, curl2
from pencil.calc import fluid_reynolds, magnetic_reynolds
from pencil.io import group_h5, dataset_h5
from pencil import read
from scipy.ndimage import gaussian_filter as gf
import os


def is_vector(key):
    """Check if the variable denoted by the label key is a vector."""
    vec = False
    for tag in ("aa", "uu", "bb", "jj", "upot", "urot", "vort", "uxb",
                "meanuu", "meanbb", "meanjj", "meanvort",
                "etadel2a", "curluxb", "curletadel2a", "advec_force",
                "fvisc", "grav", "gradp", "shear", "coriolis", "lorentz" ):
        if key in tag:
            vec = True
    return vec


def der_limits(n1, n2, m1, m2, l1, l2, nghost):
    if n1 >= nghost:
        n1shift = n1 - nghost
    else:
        n1shift = n1
    if m1 >= nghost:
        m1shift = m1 - nghost
    else:
        m1shift = m1
    if l1 >= nghost:
        l1shift = l1 - nghost
    else:
        l1shift = l1
    n2shift = n2 + nghost
    m2shift = m2 + nghost
    l2shift = l2 + nghost
    return n1shift, n2shift, m1shift, m2shift, l1shift, l2shift


def under_limits(n1, m1, l1, n1shift, m1shift, l1shift, nghost):
    if n1shift == n1:
        n1r = 0
    else:
        n1r = nghost
    if m1shift == m1:
        m1r = 0
    else:
        m1r = nghost
    if l1shift == l1:
        l1r = 0
    else:
        l1r = nghost
    return n1r, m1r, l1r


def derive_data(
    sim_path,
    src,
    dst,
    magic=["pp", "tt"],
    par=[],
    comm=None,
    gd=[],
    grp_overwrite=False,
    overwrite=False,
    rank=0,
    size=1,
    nghost=3,
    status="a",
    chunksize=1000.0,
    dtype=np.float64,
    quiet=True,
    nmin=32,
    Reynolds_shock=False,
    lmix=False,
):

    if comm:
        overwrite = False
    if isinstance(par, list):
        os.chdir(sim_path)
        par = read.param(quiet=True, conflicts_quiet=True)
    if isinstance(gd, list):
        os.chdir(sim_path)
        gd = read.grid(quiet=True)
    # get data dimensions
    nx, ny, nz = (
        src["settings"]["nx"][0],
        src["settings"]["ny"][0],
        src["settings"]["nz"][0],
    )
    mx, my, mz = (
        src["settings"]["mx"][0],
        src["settings"]["my"][0],
        src["settings"]["mz"][0],
    )
    # split data into manageable memory chunks
    dstchunksize = 8 * nx * ny * nz / 1024 * 1024
    if dstchunksize > chunksize:
        nchunks = cpu_optimal(
            nx,
            ny,
            nz,
            quiet=quiet,
            mvar=src["settings/mvar"][0],
            maux=src["settings/maux"][0],
            MBmin=chunksize,
            nmin=nmin,
            size=size,
        )[1]
    else:
        nchunks = [1, 1, 1]
    print("nchunks {}".format(nchunks))
    # for mpi split chunks across processes
    if size > 1:
        locindx = np.array_split(np.arange(nx) + nghost, nchunks[0])
        locindy = np.array_split(np.arange(ny) + nghost, nchunks[1])
        locindz = np.array_split(np.arange(nz) + nghost, nchunks[2])
        indx = [
            locindx[
                np.mod(
                    rank + int(rank / nchunks[2]) + int(rank / nchunks[1]), nchunks[0]
                )
            ]
        ]
        indy = [locindy[np.mod(rank + int(rank / nchunks[2]), nchunks[1])]]
        indz = [locindz[np.mod(rank, nchunks[2])]]
        allchunks = 1
    else:
        locindx = np.array_split(np.arange(nx) + nghost, nchunks[0])
        locindy = np.array_split(np.arange(ny) + nghost, nchunks[1])
        locindz = np.array_split(np.arange(nz) + nghost, nchunks[2])
        indx = np.array_split(np.arange(nx) + nghost, nchunks[0])
        indy = np.array_split(np.arange(ny) + nghost, nchunks[1])
        indz = np.array_split(np.arange(nz) + nghost, nchunks[2])
        allchunks = nchunks[0] * nchunks[1] * nchunks[2]
    # save time
    dataset_h5(
        dst,
        "time",
        status=status,
        data=src["time"][()],
        comm=comm,
        size=size,
        rank=rank,
        overwrite=overwrite,
        dtype=dtype,
    )
    # ensure derived variables are in a list
    if isinstance(magic, list):
        magic = magic
    else:
        magic = [magic]
    # initialise group
    group = group_h5(
        dst,
        "data",
        status="a",
        overwrite=grp_overwrite,
        comm=comm,
        rank=rank,
        size=size,
    )
    for key in magic:
        if is_vector(key):
            dataset_h5(
                group,
                key,
                status=status,
                shape=[3, mz, my, mx],
                comm=comm,
                size=size,
                rank=rank,
                overwrite=overwrite,
                dtype=dtype,
            )
            print("writing " + key + " shape {}".format([3, mz, my, mx]))
        else:
            dataset_h5(
                group,
                key,
                status=status,
                shape=[mz, my, mx],
                comm=comm,
                size=size,
                rank=rank,
                overwrite=overwrite,
                dtype=dtype,
            )
            print("writing " + key + " shape {}".format([mz, my, mx]))
        for ichunk in range(allchunks):
            for iz in [indz[np.mod(ichunk, nchunks[2])]]:
                n1, n2 = iz[0] - nghost, iz[-1] + nghost + 1
                n1out = n1 + nghost
                n2out = n2 - nghost
                varn1 = nghost
                varn2 = -nghost
                if iz[0] == locindz[0][0]:
                    n1out = 0
                    varn1 = 0
                if iz[-1] == locindz[-1][-1]:
                    n2out = n2
                    varn2 = n2
                if not quiet:
                    print("remeshing " + key + "z chunk {}".format([iz]))
                for iy in [indy[np.mod(ichunk + int(ichunk / nchunks[2]), nchunks[1])]]:
                    m1, m2 = iy[0] - nghost, iy[-1] + nghost + 1
                    m1out = m1 + nghost
                    m2out = m2 - nghost
                    varm1 = nghost
                    varm2 = -nghost
                    if iy[0] == locindy[0][0]:
                        m1out = 0
                        varm1 = 0
                    if iy[-1] == locindy[-1][-1]:
                        m2out = m2
                        varm2 = m2
                    if not quiet:
                        print("remeshing " + key + "y chunk {}".format([iy]))
                    for ix in [
                        indx[
                            np.mod(
                                ichunk
                                + int(ichunk / nchunks[2])
                                + int(ichunk / nchunks[1]),
                                nchunks[0],
                            )
                        ]
                    ]:
                        l1, l2 = ix[0] - nghost, ix[-1] + nghost + 1
                        l1out = l1 + nghost
                        l2out = l2 - nghost
                        varl1 = nghost
                        varl2 = -nghost
                        if ix[0] == locindx[0][0]:
                            l1out = 0
                            varl1 = 0
                        if ix[-1] == locindx[-1][-1]:
                            l2out = l2
                            varl2 = l2
                        if not quiet:
                            print("remeshing " + key + "x chunk {}".format([ix]))
                        var = calc_derived_data(
                            src["data"],
                            dst["data"],
                            key,
                            par,
                            gd,
                            l1,
                            l2,
                            m1,
                            m2,
                            n1,
                            n2,
                            nghost=nghost,
                            Reynolds_shock=Reynolds_shock,
                            lmix=lmix,
                        )
                        if is_vector(key):
                            dst["data"][key][
                                :, n1out:n2out, m1out:m2out, l1out:l2out
                            ] = dtype(var[:, varn1:varn2, varm1:varm2, varl1:varl2])
                        else:
                            dst["data"][key][
                                n1out:n2out, m1out:m2out, l1out:l2out
                            ] = dtype(var[varn1:varn2, varm1:varm2, varl1:varl2])


# ==============================================================================
def calc_derived_data(
    src,
    dst,
    key,
    par,
    gd,
    l1,
    l2,
    m1,
    m2,
    n1,
    n2,
    nghost=3,
    Reynolds_shock=False,
    lmix=False,
):
    """
    compute from src data and existing dst data derived data
    """
    # ==========================================================================
    def pressure(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == "pp":
            if "rho" in src.keys():
                rho = src["rho"][n1:n2, m1:m2, l1:l2]
            elif "lnrho" in src.keys():
                rho = np.exp(src["lnrho"][n1:n2, m1:m2, l1:l2])
            else:
                print("no density used setting rho=1 in pressure calculation")
                rho = 1
            if not par.gamma == 1:
                lnTT0 = np.log(par.cs0 ** 2 / (par.cp * (par.gamma - 1)))
                cv = par.cp / par.gamma
            else:
                lnTT0 = np.log(par.cs0 ** 2 / par.cp)
                cv = 1
            lnrho0 = np.log(par.rho0)
            cv1 = 1.0 / cv
            gamma_m1 = par.gamma - 1.0

            if "ss" in src.keys():
                ss = src["ss"][n1:n2, m1:m2, l1:l2]
                var = (par.cp - cv) * np.exp(
                    cv1 * ss + par.gamma * np.log(rho) - gamma_m1 * lnrho0
                )
            elif "tt" in dst.keys():
                tt = dst["tt"][n1:n2, m1:m2, l1:l2]
                var = (par.cp - cv) * tt * rho
            else:
                if "rho" in src.keys() or "lnrho" in src.keys():
                    print(
                        "no entropy or temperature using cs^2"
                        + " in pressure calculation"
                    )
                    var = rho * par.cs0 ** 2
                else:
                    print(
                        "no density or temperature," + " pressure cannot be calculated"
                    )
                    return 1
            return var

    # ==========================================================================
    def temperature(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == "tt":
            if "rho" in src.keys():
                rho = src["rho"][n1:n2, m1:m2, l1:l2]
            elif "lnrho" in src.keys():
                rho = np.exp(src["lnrho"][n1:n2, m1:m2, l1:l2])
            else:
                print("no density used setting rho=1 in temperature calculation")
                rho = 1
            if "ss" in src.keys():
                ss = src["ss"][n1:n2, m1:m2, l1:l2]
                lnrho0 = np.log(par.rho0)
                if not par.gamma == 1:
                    lnTT0 = np.log(par.cs0 ** 2 / (par.cp * (par.gamma - 1)))
                    lnTT = (
                        lnTT0
                        + par.gamma / par.cp * ss
                        + (par.gamma - 1) * (np.log(rho) - lnrho0)
                    )
                else:
                    lnTT0 = np.log(par.cs0 ** 2 / par.cp)
                    lnTT = lnTT0 + par.gamma / par.cp * ss
            else:
                lnTT0 = np.log(par.cs0 ** 2 / (par.cp * (par.gamma - 1)))
                lnTT = (par.gamma - 1) * (np.log(rho) - lnrho0) + lnTT0
            return np.exp(lnTT)

    # ======================================================================
    def Re_number(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == "Re":
            n1shift, n2shift, m1shift, m2shift, l1shift, l2shift = der_limits(
                n1, n2, m1, m2, l1, l2, nghost
            )
            uu = np.array(
                [
                    src["ux"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift],
                    src["uy"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift],
                    src["uz"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift],
                ]
            )
            if "rho" in src.keys():
                lnrho = np.log(
                    src["rho"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift]
                )
            elif "lnrho" in src.keys():
                lnrho = src["lnrho"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift]
            else:
                lnrho = list()
            if "shock" in src.keys() and Reynolds_shock:
                shock = src["shock"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift]
            else:
                shock = list()
            var = fluid_reynolds(uu, par, gd, lnrho=lnrho, shock=shock, lmix=lmix)
            n1r, m1r, l1r = under_limits(n1, m1, l1, n1shift, m1shift, l1shift, nghost)
            return var[n1r : n2 - n1 + n1r, m1r : m2 - m1 + m1r, l1r : l2 - l1 + l1r]

    # ======================================================================
    def Rm_number(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == "Rm":
            n1shift, n2shift, m1shift, m2shift, l1shift, l2shift = der_limits(
                n1, n2, m1, m2, l1, l2, nghost
            )
            uu = np.array(
                [
                    src["ux"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift],
                    src["uy"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift],
                    src["uz"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift],
                ]
            )
            aa = np.array(
                [
                    src["ax"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift],
                    src["ay"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift],
                    src["az"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift],
                ]
            )
            if "bb" in dst.keys():
                bb = dst["bb"][:, n1shift:n2shift, m1shift:m2shift, l1shift:l2shift]
            else:
                bb = bfield(
                    src,
                    dst,
                    "bb",
                    par,
                    gd,
                    l1shift,
                    l2shift,
                    m1shift,
                    m2shift,
                    n1shift,
                    n2shift,
                    nghost,
                )
            if "jj" in dst.keys():
                jj = dst["jj"][:, n1shift:n2shift, m1shift:m2shift, l1shift:l2shift]
            else:
                jj = current(
                    src,
                    dst,
                    "jj",
                    par,
                    gd,
                    l1shift,
                    l2shift,
                    m1shift,
                    m2shift,
                    n1shift,
                    n2shift,
                    nghost,
                )
            var = magnetic_reynolds(uu, par, gd, aa=aa, bb=bb, jj=jj, lmix=lmix)
            n1r, m1r, l1r = under_limits(n1, m1, l1, n1shift, m1shift, l1shift, nghost)
            return var[n1r : n2 - n1 + n1r, m1r : m2 - m1 + m1r, l1r : l2 - l1 + l1r]

    # ======================================================================
    def Pm_number(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == "Pm":
            if "Re" in dst.keys():
                Re = dst["Re"][n1:n2, m1:m2, l1:l2]
            else:
                Re = Re_number(src, dst, "Re", par, gd, l1, l2, m1, m2, n1, n2, nghost)
            if "Rm" in dst.keys():
                Rm = dst["Rm"][n1:n2, m1:m2, l1:l2]
            else:
                Rm = Rm_number(src, dst, "Rm", par, gd, l1, l2, m1, m2, n1, n2, nghost)
            if Re.max() > 0:
                Re[np.where(Re == 0)] = Re[np.where(Re > 0)].min()
            else:
                if Rm.max() > 0:
                    Re[:] = Rm.min()
                else:
                    Re[:] = 1
            return Rm / Re

    # ======================================================================
    def rot_flow(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == "urot":
            uu = np.array(
                [
                    src["ux"][n1:n2, m1:m2, l1:l2],
                    src["uy"][n1:n2, m1:m2, l1:l2],
                    src["uz"][n1:n2, m1:m2, l1:l2],
                ]
            )
            var = helmholtz_fft(uu, gd, par, rot=True, pot=False)
            return var

    # ======================================================================
    def pot_flow(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == "upot":
            uu = np.array(
                [
                    src["ux"][n1:n2, m1:m2, l1:l2],
                    src["uy"][n1:n2, m1:m2, l1:l2],
                    src["uz"][n1:n2, m1:m2, l1:l2],
                ]
            )
            var = helmholtz_fft(uu, gd, par, pot=True, rot=False)
            return var

    # ======================================================================
    def Mach_cs(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == "Ms":
            uu = np.array(
                [
                    src["ux"][n1:n2, m1:m2, l1:l2],
                    src["uy"][n1:n2, m1:m2, l1:l2],
                    src["uz"][n1:n2, m1:m2, l1:l2],
                ]
            )
            if "tt" in dst.keys():
                tt = dst["tt"][n1:n2, m1:m2, l1:l2]
            else:
                tt = temperature(
                    src, dst, "tt", par, gd, l1, l2, m1, m2, n1, n2, nghost
                )
            if not par.gamma == 1:
                cs2 = par.cp * (par.gamma - 1) * tt
            else:
                cs2 = par.cp * tt
            print("tt min {} max {}".format(tt.min(), tt.max()))
            var = np.sqrt(dot2(uu) / cs2)
            return var

    # ======================================================================
    def Mach_Av(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == "Ma":
            if "bb" in dst.keys():
                bb = dst["bb"][:, n1:n2, m1:m2, l1:l2]
            else:
                bb = bfield(src, dst, "bb", par, gd, l1, l2, m1, m2, n1, n2, nghost)
            if "rho" in src.keys():
                rho = src["rho"][n1:n2, m1:m2, l1:l2]
            elif "lnrho" in src.keys():
                rho = np.exp(src["lnrho"][n1:n2, m1:m2, l1:l2])
            else:
                print("no density used setting rho=1 in pressure calculation")
                rho = 1
            var = np.sqrt(dot2(bb) / (par.mu0 * rho))
            return var

    # ======================================================================
    def mag_pressure(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == "pb":
            if "bb" in dst.keys():
                bb = dst["bb"][:, n1:n2, m1:m2, l1:l2]
            else:
                bb = bfield(src, dst, "bb", par, gd, l1, l2, m1, m2, n1, n2, nghost)
            var = 0.5 * dot2(bb) / par.mu0
            return var

    # ======================================================================
    def vorticity(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == "vort":
            n1shift, n2shift, m1shift, m2shift, l1shift, l2shift = der_limits(
                n1, n2, m1, m2, l1, l2, nghost
            )
            uu = np.array(
                [
                    src["ux"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift],
                    src["uy"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift],
                    src["uz"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift],
                ]
            )
            var = curl(uu, gd.dx, gd.dy, gd.dz)
            n1r, m1r, l1r = under_limits(n1, m1, l1, n1shift, m1shift, l1shift, nghost)
            return var[:, n1r : n2 - n1 + n1r, m1r : m2 - m1 + m1r, l1r : l2 - l1 + l1r]

    # ======================================================================
    def bfield(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == "bb":
            n1shift, n2shift, m1shift, m2shift, l1shift, l2shift = der_limits(
                n1, n2, m1, m2, l1, l2, nghost
            )
            aa = np.array(
                [
                    src["ax"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift],
                    src["ay"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift],
                    src["az"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift],
                ]
            )
            var = curl(aa, gd.dx, gd.dy, gd.dz)
            n1r, m1r, l1r = under_limits(n1, m1, l1, n1shift, m1shift, l1shift, nghost)
            return var[:, n1r : n2 - n1 + n1r, m1r : m2 - m1 + m1r, l1r : l2 - l1 + l1r]

    # ======================================================================
    def current(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == "jj":
            n1shift, n2shift, m1shift, m2shift, l1shift, l2shift = der_limits(
                n1, n2, m1, m2, l1, l2, nghost
            )
            aa = np.array(
                [
                    src["ax"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift],
                    src["ay"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift],
                    src["az"][n1shift:n2shift, m1shift:m2shift, l1shift:l2shift],
                ]
            )
            var = curl2(aa, gd.dx, gd.dy, gd.dz)
            n1r, m1r, l1r = under_limits(n1, m1, l1, n1shift, m1shift, l1shift, nghost)
            return var[:, n1r : n2 - n1 + n1r, m1r : m2 - m1 + m1r, l1r : l2 - l1 + l1r]

    # ======================================================================
    def kin_helicity(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == "ou":
            uu = np.array(
                [
                    src["ux"][n1:n2, m1:m2, l1:l2],
                    src["uy"][n1:n2, m1:m2, l1:l2],
                    src["uz"][n1:n2, m1:m2, l1:l2],
                ]
            )
            if "vort" in dst.keys():
                oo = dst["vort"][:, n1:n2, m1:m2, l1:l2]
            else:
                oo = vorticity(
                    src, dst, "vort", par, gd, l1, l2, m1, m2, n1, n2, nghost
                )
            var = dot(uu, oo)
            return var

    # ======================================================================
    def mag_helicity(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == "ab":
            aa = np.array(
                [
                    src["ax"][n1:n2, m1:m2, l1:l2],
                    src["ay"][n1:n2, m1:m2, l1:l2],
                    src["az"][n1:n2, m1:m2, l1:l2],
                ]
            )
            if "bb" in dst.keys():
                bb = dst["bb"][:, n1:n2, m1:m2, l1:l2]
            else:
                bb = bfield(src, dst, "bb", par, gd, l1, l2, m1, m2, n1, n2, nghost)
            var = dot(aa, bb)
            return var

    #==========================================================================
    def urand(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost):
        if key == "u2rand":
            sigma = 0.02
            if "sigma" in par.keys:
                sigma = par.sigma
            n1shift,n2shift,m1shift,m2shift,l1shift,l2shift=der_limits(
                                                     n1,n2,m1,m2,l1,l2,nghost)
            uu = np.array([
                              src["ux"][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                              src["uy"][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift],
                              src["uz"][n1shift:n2shift,m1shift:m2shift,l1shift:l2shift]
                              ])
            umean = uu.copy()
            urand = uu.copy()
            dx = min(gd.dx,gd.dy,gd.dz)
            sigma_dx = sigma/dx
            print("sigma_dx {:.2f} sigma {}".format(sigma_dx,sigma))
            for j in range(3):
                umean[j] = gf(uu[j], sigma_dx, mode="reflect")
                urand[j] = uu[j] - umean[j]
            var = dot2(urand)
            n1r,m1r,l1r = under_limits(n1,m1,l1,n1shift,m1shift,l1shift,nghost)
            return var[n1r:n2-n1+n1r,m1r:m2-m1+m1r,l1r:l2-l1+l1r]

    # ==========================================================================
    def calc_derived_item(key):
        case = {
            "urot": rot_flow(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
            "upot": pot_flow(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
            "vort": vorticity(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
            "ou": kin_helicity(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
            "tt": temperature(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
            "pp": pressure(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
            "Re": Re_number(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
            "Ms": Mach_cs(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
            "bb": bfield(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
            "pb": mag_pressure(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
            "ab": mag_helicity(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
            "jj": current(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
            "Ma": Mach_Av(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
            "Rm": Rm_number(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
            "Pm": Pm_number(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
            "u2rand": urand(src, dst, key, par, gd, l1, l2, m1, m2, n1, n2, nghost),
        }
        func = case.get(key, lambda: "No function for " + key)
        return func
        # ======================================================================

    return calc_derived_item(key)
