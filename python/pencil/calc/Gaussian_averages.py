# derived_h5.py
#
# 05-may-20
# Author: F. Gent (fred.gent.ncl@gmail.com).
#
""" Derive auxilliary data and other diagnostics from var.h5 file and
    save to new h5 file

    uses:
      compute 'data' arrays of size [nz,ny,nx] as required
      store 'time' of snapshot
      compute 'masks' for example by temperature phase
      compute summary statistics 'stats'
      compute 'structure' functions as required
"""
import numpy as np
from scipy.interpolate import interp1d
from pencil.math import dot, dot2, natural_sort, helmholtz_fft, cpu_optimal
from pencil.math.derivatives import curl, div, curl2, grad
from pencil.calc import fluid_reynolds, magnetic_reynolds
from pencil.io import open_h5, group_h5, dataset_h5
from fileinput import input
from sys import stdout
import subprocess as sub
from pencil import read
import os

from ..math import dot, dot2, cross
from ..math.derivatives import div, curl, curl2, grad, del2, del6
import numpy as np
from scipy.ndimage.filters import gaussian_filter1d, gaussian_filter
from pencil.ism_dyn import is_vector


def kernel_smooth(
    sim_path,
    src,
    dst,
    magic=["meanuu"],
    par=[],
    comm=None,
    gd=[],
    grp_overwrite=False,
    overwrite=False,
    rank=0,
    size=1,
    nghost=3,
    kernel=1.,
    status="a",
    chunksize=1000.0,
    dtype=np.float64,
    quiet=True,
    nmin=32,
    typ='piecewise',
    mode=list(),
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
    # extend gost zones to include up to 1.5 * kernel length)
    dx = max(src['grid/dx'][()],src['grid/dy'][()],src['grid/dz'][()])
    nkernel = np.int(2.5*kernel/dx)
    sigma = kernel / dx
    if not quiet:
        print('sigma {:.2f}, kernel {:.2f}, dx {:.4f}'.format(sigma,kernel,dx))
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
    if 1 in nchunks:
        mode = ["reflect","reflect","reflect"]
        for ich in range(3):
            if nchunks[ich] == 1:
                mode[2-ich] = "wrap"
            if mode[2-ich] == "reflect":
                typ = "piecewise"
            else:
                typ = "all"
    if not quiet:
        print('mode:',mode,'typ:',typ)
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
                key+str(nkernel),
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
                key+str(nkernel),
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
                if nchunks[2] == 1:
                    zextra = nghost
                else:
                    zextra = nkernel + nghost
                n1, n2 = iz[0] - zextra, iz[-1] + zextra + 1
                lindz = np.arange(n1,n2)
                n1out = n1 + zextra
                n2out = n2 - zextra
                varn1 = zextra
                varn2 = -zextra
                if iz[0] == locindz[0][0]:
                    n1out = 0
                    varn1 = zextra - nghost
                if iz[-1] == locindz[-1][-1]:
                    n2out = n2 - zextra + nghost
                    varn2 = n2 - n1 - zextra + nghost
                if n1 < 0:
                    lindz[np.where(lindz<nghost)[0]] += nz
                if n2 > mz-1:
                    lindz[np.where(lindz>mz-1-nghost)[0]] -= nz
                if not quiet:
                    print('n1out {},n2out {},varn1 {},varn2 {},zextra {}'.format(n1out,n2out,varn1,varn2,zextra))
                for iy in [indy[np.mod(ichunk + int(ichunk / nchunks[2]), nchunks[1])]]:
                    if nchunks[1] == 1:
                        yextra = nghost
                    else:
                        yextra = nkernel + nghost
                    m1, m2 = iy[0] - yextra, iy[-1] + yextra + 1
                    lindy = np.arange(m1,m2)
                    m1out = m1 + yextra
                    m2out = m2  - yextra
                    varm1 = yextra
                    varm2 = -yextra
                    if iy[0] == locindy[0][0]:
                        m1out = 0
                        varm1 = yextra - nghost
                    if iy[-1] == locindy[-1][-1]:
                        m2out = m2 - yextra + nghost
                        varm2 = m2 - m1 - yextra + nghost
                    if m1 < 0:
                        lindy[np.where(lindy<0)[0]] += ny
                    if m2 > my-1:
                        lindy[np.where(lindy>my-1-nghost)[0]] -= ny
                    if not quiet:
                        print('m1out {},m2out {},varm1 {},varm2 {},yextra {}'.format(m1out,m2out,varm1,varm2,yextra))
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
                        if nchunks[1] == 1:
                            xextra = nghost
                        else:
                            xextra = nkernel + nghost
                        l1, l2 = ix[0] - xextra,  ix[-1] + xextra + 1
                        lindx = np.arange(l1,l2)
                        l1out = l1 + xextra
                        l2out = l2 - xextra
                        varl1 = xextra
                        varl2 = -xextra
                        if ix[0] == locindx[0][0]:
                            l1out = 0
                            varl1 = xextra - nghost
                        if ix[-1] == locindx[-1][-1]:
                            l2out = l2 - xextra + nghost
                            varl2 = l2 - l1 - xextra + nghost
                        if l1 < 0:
                            lindx[np.where(lindx<0)[0]] += nx
                        if l2 > mx-1:
                            lindx[np.where(lindx>mx-1-nghost)[0]] -= nx
                        if not quiet:
                            print('l1out {},l2out {},varl1 {},varl2 {},xextra {}'.format(l1out,l2out,varl1,varl2,xextra))
                        if not quiet:
                            print("remeshing " + key + " chunk {}".format([iz, iy, ix]))
                            print('sending ichunk {} with index ranges {}'.format(ichunk, [n1, n2, m1, m2, l1, l2]))
                        var = smoothed_data(
                            src["data"],
                            dst["data"],
                            key,
                            par,
                            gd,
                            lindx,
                            lindy,
                            lindz,
                            nghost,
                            sigma,
                            typ,
                            quiet,
                            mode
                        )
                        if not quiet:
                            print('ichunk {}, var min {:.1e}, var max {:.1e}'.format(ichunk, var.min(), var.max()))
                        # print('var shape {}'.format(var.shape))
                        # if not quiet:
                        #    print('writing '+key+
                        #                   ' shape {} chunk {}'.format(
                        #                         var.shape, [iz,iy,ix]))
                        if not quiet:
                            print('ichunk: out indices {}'.format([n1out, n2out, m1out, m2out, l1out, l2out]))
                        if is_vector(key):
                            dst["data"][key+str(nkernel)][
                                :, n1out:n2out, m1out:m2out, l1out:l2out
                            ] = dtype(var[:, varn1:varn2, varm1:varm2, varl1:varl2])
                        else:
                            dst["data"][key+str(nkernel)][
                                n1out:n2out, m1out:m2out, l1out:l2out
                            ] = dtype(var[varn1:varn2, varm1:varm2, varl1:varl2])

# ======================================================================
def load_dataset(src, key, lindx, lindy, lindz, nghost):

    if is_vector(key):
        var = np.empty([3,lindz.size,lindy.size,lindx.size])
    else:
        var = np.empty([lindz.size,lindy.size,lindx.size])

    nz, ny, nx = src[key][()].shape[-3], src[key][()].shape[-2], src[key][()].shape[-1]
    if var.shape == src[key][()].shape:
        var = src[key][()]
    else:
        n3, n4, m3, m4, l3, l4 = None, None, None, None, None, None
        # define array indices in z direction
        if var.shape[-3] >= nz:
            n1, n2 = 0, lindz.size
            nn1, nn2 = 0, lindz.size
        else:
            if lindz[0] == np.min(lindz):
                n1, n2 = lindz[0], lindz[-1]+1
                nn1, nn2 = 0, lindz.size
            else:
                n1 = lindz[0]
                n2 = nz - nghost
                n3 = nghost
                n4 = lindz[-1] + 1
                nn1 = 0
                nn2 = n2 - n1
                nn3 = n2 - n1
                if lindz[-1] > nz:
                    nn4 = lindz.size
                else:
                    nn4 = np.mod(lindz.size, nz)
        # define array indices in y direction
        if var.shape[-2] >= ny:
            m1, m2 =  0, lindy.size
            mm1, mm2 = 0, lindy.size
        else:
            if lindy[0] == np.min(lindy):
                m1, m2 = lindy[0], lindy[-1]+1
                mm1, mm2 = 0, lindy.size
            else:
                m1 = lindy[0]
                m2 = ny - nghost
                m3 = nghost
                m4 = lindy[-1] + 1
                mm1 = 0
                mm2 = m2 - m1
                mm3 = m2 - m1
                if lindy[-1] > ny:
                    mm4 = lindy.size
                else:
                    mm4 = np.mod(lindy.size, ny)
        # define array indices in x direction
        if var.shape[-1] >= nx:
            l1, l2 = 0, lindx.size
            ll1, ll2 = 0, lindx.size
        else:
            if lindx[0] == np.min(lindx):
                l1, l2 = lindx[0], lindx[-1]+1
                ll1, ll2 = 0, lindx.size
            else:
                l1 = lindx[0]
                l2 = nx - nghost
                l3 = nghost
                l4 = lindx[-1] + 1
                ll1 = 0
                ll2 = l2 - l1
                ll3 = l2 - l1
                if lindx[-1] > nx:
                    ll4 = lindx.size
                else:
                    ll4 = np.mod(lindx.size, nx)
        var[:,nn1:nn2,mm1:mm2,ll1:ll2] = src[key][:,n1:n2,m1:m2,l1:l2]
        if l3:
            var[:,nn1:nn2,mm1:mm2,ll3:ll4] = src[key][:,n1:n2,m1:m2,l3:l4]
            if m3:
                var[:,nn1:nn2,mm3:mm4,ll3:ll4] = src[key][:,n1:n2,m3:m4,l3:l4]
                var[:,nn1:nn2,mm3:mm4,ll1:ll2] = src[key][:,n1:n2,m3:m4,l1:l2]
                if n3:
                    var[:,nn3:nn4,mm3:mm4,ll3:ll4] = src[key][:,n3:n4,m3:m4,l3:l4]
                    var[:,nn3:nn4,mm1:mm2,ll3:ll4] = src[key][:,n3:n4,m1:m2,l3:l4]
                    var[:,nn3:nn4,mm3:mm4,ll1:ll2] = src[key][:,n3:n4,m3:m4,l1:l2]
                    var[:,nn3:nn4,mm1:mm2,ll1:ll2] = src[key][:,n3:n4,m1:m2,l1:l2]
        elif m3:
            var[:,nn1:nn2,mm3:mm4,ll1:ll2] = src[key][:,n1:n2,m3:m4,l1:l2]
            if n3:
                var[:,nn3:nn4,mm3:mm4,ll1:ll2] = src[key][:,n3:n4,m3:m4,l1:l2]
                var[:,nn3:nn4,mm1:mm2,ll1:ll2] = src[key][:,n3:n4,m1:m2,l1:l2]
        else:
            if n3:
                var[:,nn3:nn4,mm1:mm2,ll1:ll2] = src[key][:,n3:n4,m1:m2,l1:l2]

    return var

# ======================================================================
def gauss_3Dsmooth(arr, sigma=1., typ="all", quiet=True, mode=list()):


    # if vector alternate boundary conditions optional
    if typ == 'piecewise':
        if len(mode) == 3:
            arr0 = gaussian_filter1d(arr, axis=-3, sigma=sigma,mode=mode[0])
            if not quiet:
                print('gaussian_filter1d axis -3 complete')
            arr1 = gaussian_filter1d(arr0, axis=-2, sigma=sigma,mode=mode[1])
            if not quiet:
                print('gaussian_filter1d axis -2 complete')
            sm_arr = gaussian_filter1d(arr1, axis=-1, sigma=sigma,mode=mode[2])
            if not quiet:
                print('gaussian_filter1d axis -1 complete')
        else:
            # default mode is reflect
            if len(mode) < 1:
                mode.append("reflect")
                print("Warning: piecewise requires mode list containing 3 options for axes")
            sm_arr = gaussian_filter(arr, sigma=[sigma,sigma,sigma], mode=mode[0])
    elif typ == "all":
        # default mode is reflect
        if len(mode) < 1:
            mode.append("reflect")
        sm_arr = gaussian_filter(arr, sigma=[sigma,sigma,sigma], mode=mode[0])
        if not quiet:
            print('gaussian_filter complete')
    return sm_arr

# ======================================================================
def smoothed_data(
                  src,
                  dst,
                  key,
                  par,
                  gd,
                  lindx,
                  lindy,
                  lindz,
                  nghost,
                  sigma,
                  typ,
                  quiet,
                  mode
                 ):


    # ======================================================================
    def mean_flow(src, dst, key, par, gd, lindx, lindy, lindz, nghost, sigma, typ, quiet, mode):
        if key == "meanuu":
            ux = load_dataset(src, 'ux', lindx, lindy, lindz, nghost)
            uy = load_dataset(src, 'uy', lindx, lindy, lindz, nghost)
            uz = load_dataset(src, 'uz', lindx, lindy, lindz, nghost)
            uu = np.array([ux,uy,uz])
            var = uu.copy()
            for j in range(3):
                var[j] = gauss_3Dsmooth(uu[j], sigma=sigma, typ=typ, quiet=quiet, mode=mode)
            return var

    # ======================================================================
    def mean_density(src, dst, key, par, gd, lindx,lindy, lindz, nghost, sigma, typ, quiet, mode):
        if key == "meanrho":
            if 'rho' in src.keys():
                rho = load_dataset(src, 'rho', lindx, lindy, lindz, nghost)
            elif 'lnrho' in src.keys():
                var = load_dataset(src, 'rho', lindx, lindy, lindz, nghost)
                rho = np.exp(var)
            var = gauss_3Dsmooth(rho, sigma=sigma, typ=typ, quiet=quiet, mode=mode)
            return var

    # ======================================================================
    def mean_mag_pressure(src, dst, key, par, gd, lindx, lindy, lindz, nghost, sigma, typ, quiet, mode):
        if key == "meaneb":
            if "meanbb" in dst.keys():
                bb = load_dataset(dst, 'meanbb', lindx, lindy, lindz, nghost)
            else:
                bb = mean_bfield(src, dst, 'meanbb', par, gd, lindx, lindy, lindz, nghost, sigma, typ, quiet, mode)
            eb = 0.5 * dot2(bb) / par.mu0
            #var = gauss_3Dsmooth(eb, sigma=sigma, typ=typ, mode=mode)
            return eb

    # ======================================================================
    def mean_bfield(src, dst, key, par, gd, lindx, lindy, lindz, nghost, sigma, typ, quiet, mode):
        if key == "meanbb":
            if "bb" in dst.keys():
                bb = load_dataset(dst, 'bb', lindx, lindy, lindz, nghost)
            else:
                ax = load_dataset(src, 'ax', lindx, lindy, lindz, nghost)
                ay = load_dataset(src, 'ay', lindx, lindy, lindz, nghost)
                az = load_dataset(src, 'az', lindx, lindy, lindz, nghost)
                aa = np.array([ax,ay,az])
                bb  = curl(aa, gd.dx, gd.dy, gd.dz)
            if not quiet:
                print('bb range {:.2e} to {:.2e}'.format(bb.min(),bb.max()))
            var = bb.copy()
            for j in range(3):
                if not quiet:
                    print('gaussian smoothing:',j)
                var[j] = gauss_3Dsmooth(bb[j], sigma=sigma, typ=typ, quiet=quiet, mode=mode)
            if not quiet:
                print('meanbb range {:.2e} to {:.2e}'.format(var.min(),var.max()))

            return var

    # ======================================================================
    def mean_mag_helicity(src, dst, key, par, gd, lindx, lindy, lindz, nghost, sigma, typ, quiet, mode):
        if key == "meanab":
            if 'ab' in dst.keys():
                ab = load_dataset(dst, 'ab', lindx, lindy, lindz, nghost)
            else:
                ax = load_dataset(src, 'ax', lindx, lindy, lindz, nghost)
                ay = load_dataset(src, 'ay', lindx, lindy, lindz, nghost)
                az = load_dataset(src, 'az', lindx, lindy, lindz, nghost)
                aa = np.array([ax,ay,az])
                if "bb" in dst.keys():
                    bb = load_dataset(dst, 'bb', lindx, lindy, lindz, nghost)
                else:
                    bb  = curl(aa, gd.dx, gd.dy, gd.dz)
                ab = dot(aa,bb)
            var = gauss_3Dsmooth(ab, sigma=sigma, typ=typ, quiet=quiet, mode=mode)
            return var

    # ======================================================================
    def mean_cur_helicity(src, dst, key, par, gd, lindx, lindy, lindz, nghost, sigma, typ, quiet, mode):
        if key == "meanjb":
            if 'jb' in dst.keys():
                jb = load_dataset(dst, 'jb', lindx, lindy, lindz, nghost)
            else:
                if "bb" in dst.keys():
                    bb = load_dataset(dst, 'bb', lindx, lindy, lindz, nghost)
                else:
                    ax = load_dataset(src, 'ax', lindx, lindy, lindz, nghost)
                    ay = load_dataset(src, 'ay', lindx, lindy, lindz, nghost)
                    az = load_dataset(src, 'az', lindx, lindy, lindz, nghost)
                    aa = np.array([ax,ay,az])
                    bb  = curl(aa, gd.dx, gd.dy, gd.dz)
                if "jj" in dst.keys():
                    jj = load_dataset(dst, 'jj', lindx, lindy, lindz, nghost)
                else:
                    jj  = curl(bb, gd.dx, gd.dy, gd.dz)
                jb = dot(jj,bb)
            var = gauss_3Dsmooth(jb, sigma=sigma, typ=typ, quiet=quiet, mode=mode)
            return var

    # ======================================================================
    def mean_kin_helicity(src, dst, key, par, gd, lindx, lindy, lindz, nghost, sigma, typ, quiet, mode):
        if key == "meanou":
            if 'ou' in dst.keys():
                ou = load_dataset(dst, 'ou', lindx, lindy, lindz, nghost)
            else:
                ux = load_dataset(src, 'ux', lindx, lindy, lindz, nghost)
                uy = load_dataset(src, 'uy', lindx, lindy, lindz, nghost)
                uz = load_dataset(src, 'uz', lindx, lindy, lindz, nghost)
                uu = np.array([ux,uy,uz])
                if 'vort' in dst.keys():
                    vort = load_dataset(dst, 'vort', lindx, lindy, lindz, nghost)
                else:
                    vort  = curl(uu, gd.dx, gd.dy, gd.dz)
                ou = dot(vort,uu)
            var = gauss_3Dsmooth(ou, sigma=sigma, typ=typ, quiet=quiet, mode=mode)
            return var

    # ==========================================================================
    def calc_smoothed_item(key):
        case = {
            "meanuu":         mean_flow(src, dst, key, par, gd, lindx, lindy, lindz, nghost, sigma, typ, quiet, mode),
            "meanrho":     mean_density(src, dst, key, par, gd, lindx, lindy, lindz, nghost, sigma, typ, quiet, mode),
            "meanbb":       mean_bfield(src, dst, key, par, gd, lindx, lindy, lindz, nghost, sigma, typ, quiet, mode),
            "meaneb": mean_mag_pressure(src, dst, key, par, gd, lindx, lindy, lindz, nghost, sigma, typ, quiet, mode),
            "meanab": mean_mag_helicity(src, dst, key, par, gd, lindx, lindy, lindz, nghost, sigma, typ, quiet, mode),
            "meanjb": mean_cur_helicity(src, dst, key, par, gd, lindx, lindy, lindz, nghost, sigma, typ, quiet, mode),
            "meanou": mean_kin_helicity(src, dst, key, par, gd, lindx, lindy, lindz, nghost, sigma, typ, quiet, mode),
        }
        func = case.get(key, lambda: "No function for " + key)
        return func
        # ======================================================================

    return calc_smoothed_item(key)
