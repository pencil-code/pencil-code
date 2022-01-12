# struct.py
#
# Calculate structure functions from arrays in D1, D2 or D3 or time
#

"""
Contains the functions to extract structure functions from data.
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import time
import sys

# ------------------------------------------------------------------------------
# spatial fitting function
def fit_gaussm1(ell, sigma, L0):
    """Fit for 2nd order structure function, Eq. 6 DOI 10.3847/1538-4357/aa93e7
    correlation length = sqrt(pi/2)L0
    """
    return 2 * sigma ** 2 * (1 - np.exp(-0.5 * ell ** 2 / L0 ** 2))


# ------------------------------------------------------------------------------
# standard gaussian fitting function
def fit_gauss(ell, sigma, L0):
    """Fit for 2nd order structure function, Eq. 6 DOI 10.3847/1538-4357/aa93e7
    correlation length = sqrt(pi/2)L0
    """
    return sigma * np.exp(-0.5 * ell ** 2 / L0 ** 2)


# ------------------------------------------------------------------------------
# temporal fitting function
def fit_expm1(t, sigma, L0):
    """Fit for 2nd order structure function, Eq. 5 DOI 10.3847/1538-4357/aa93e7
    correlation length = L0
    """
    return 2 * sigma ** 2 * (1 - np.exp(-t / L0))


# ------------------------------------------------------------------------------
# standard exponential fitting function
def fit_exp(t, sigma, L0):
    """Fit for 2nd order structure function, Eq. 5 DOI 10.3847/1538-4357/aa93e7
    correlation length = L0
    """
    return sigma * np.exp(t * L0)


# ------------------------------------------------------------------------------


def space_struct(
    arr,
    dims=[[], [], []],
    Dorder=2,
    dirs="zyx",
    lperi=(True, True, True),
    lplot=True,
    InitialGuess=[1.0, 1.0],
    deltay=0.0,
    chunksize=1000.0,
    figname=None,
    dlabel="data",
    quiet=True,
    downsample=[1, 1, 1],
    maxl=[-1, -1, -1],
    loopsample=[1, 1, 1],
):
    """
    Calculate the structure function in space.

    call signature:

    space_struct(arr, dims=[[],[],[]], Dorder=2, dirs='zyx',
                 lperi=(True,True,True), lplot=True,
                 InitialGuess = [1.0,1.0],
                 deltay=0.,
                 figname = None, dlabel = 'data',
                 quiet=False)

    Keyword arguments:

    *arr*:
      numpy data array of shape, e.g., [nz,ny,nz].

    *z*:
      z dimension full array to replace default empty list.

    *y*:
      y dimension full array to replace default empty list.

    *x*:
      x dimension full array to replace default empty list.

    *Dorder*:
      Order of structure function, by default 2nd.

    *lperi*:
      Flag indicates for each dimension, whether boundary is periodic.

    *lplot*
      Flag to plot the result.

    *InitialGuess*:
      Initial parameters for curve fitting, default [1.0,1.0].

    *figname*:
      String name for plot, if saving rather than only display.

    *dlabel*:
      legend label for data.

    *quiet*
      Flag for switching off output.

    Returns
    -------
    D1 Structure function, length array, fit parameters [sigma, L0].

    Notes
    -----
    Must provide list of dimension arrays to match dirs required.

    Examples
    --------
    >>> D,ell,fit = pc.math.stats.space_struct(var.rho, dims=[gd.z,gd.y,gd.x],
                                               dirs='zyx')
    >>> D, ell, fit
    [0, ...,], [0, 0.0041, ...], [2.5,0.04]
    """
    start_time = time.time()
    ldirmatch = True
    # determine dimensions of dataset
    arrshape = arr.shape
    if not quiet:
        print("arrshape", arrshape)
    D3 = len(arrshape) == 3
    D2 = len(arrshape) == 2
    D1 = len(arrshape) == 1
    if not quiet:
        print("D3 {}, D2 {}, D1 {}".format(D3, D2, D1))
    if D3:
        arrchunksize = 8 * arrshape[0] * arrshape[1] * arrshape[2] / 1024 / 1024
    elif D2:
        arrchunksize = 8 * arrshape[0] * arrshape[1] / 1024 / 1024
    else:
        arrchunksize = 8 * arrshape[0] / 1024 / 1024

    # include which dimensions in function D?
    axis = 0
    Dshape = []
    Daxis = []
    if "z" in dirs:
        # first dimension of function D is z
        nz = int(len(dims[axis]) / 2)
        Dshape.append(nz)
        Daxis.append(axis)
        if not quiet:
            print("z in dirs {}: Dshape {}, Daxis {}".format(dirs, Dshape, Daxis))
        axis += 1
        ldirmatch = 2 * nz == arrshape[0]
        if not ldirmatch:
            print(
                "Please correct size dims(z) "
                + "{} does not equal arr shape {}".format(2 * nz, arrshape[0])
            )
            sys.exit()
    else:
        # first dimension of function D is not z
        nz = 0
    if "y" in dirs:
        # next dimension of function D is y
        if D3 and "z" not in dirs:
            axis += 1
        ny = int(len(dims[axis]) / 2)
        Dshape.append(ny)
        Daxis.append(axis)
        if not quiet:
            print("y in dirs {}: Dshape {}, Daxis {}".format(dirs, Dshape, Daxis))
        axis += 1
        ldirmatch = 2 * ny == arrshape[axis]
        if not ldirmatch:
            print(
                "Please correct size dims(y) "
                + "{} does not equal arr shape {}".format(2 * ny, arrshape[axis])
            )
            sys.exit()
    else:
        # y is not dimension of function D
        ny = 0
    if "x" in dirs:
        # next dimension of function D is x
        if D3 and "y" not in dirs:
            if "z" not in dirs:
                axis += 2
            else:
                axis += 1
        if D2 and "y" not in dirs and "z" not in dirs:
            axis += 1
        nx = int(len(dims[axis]) / 2)
        Dshape.append(nx)
        Daxis.append(axis)
        if not quiet:
            print("x in dirs {}: Dshape {}, Daxis {}".format(dirs, Dshape, Daxis))
        ldirmatch = 2 * nx == arrshape[axis]
        if not ldirmatch:
            print(
                "Please correct size dims(x) "
                + "{} does not equal arr shape {}".format(2 * nx, arrshape[axis])
            )
            sys.exit()
    else:
        # x is not dimension of function D
        nx = 0
    if len(Dshape) == 0:
        print(
            "length arrays are empty lists; at least 1 dimension from"
            " z, y or x must be provided"
        )
    if len(Dshape) > len(arrshape):
        print(
            "N dirs {} cannot exceed N dims {} of arr".format(
                len(Dshape), len(arrshape)
            )
        )
    # D container of correlations ell container of separation lengths N sample count
    D = np.zeros(Dshape)
    N = np.zeros(Dshape)
    ell = np.zeros(Dshape)
    D[:] = np.nan
    N[:] = np.nan
    ell[:] = np.nan
    # skip calculating D when ell > lmax, based on shortest domain edge
    lmax = 1e8
    if D1:
        lmax = min(
            lmax, (Dshape[0] + 1) / 2 / Dshape[0] * abs(dims[0][-1] - dims[0][0])
        )
        if maxl[0] < 0:
            maxl[0] = Dshape[0]
        D[0] = 0.0
        try:
            N[0] = arr[arr.mask == False].size / arr.size
        except:
            N[0] = 1.0
        ell[0] = 0.0
    if D2:
        lmax = np.min(
            [
                lmax,
                np.sqrt(
                    ((Dshape[0] + 1) / 2 / Dshape[0] * abs(dims[0][-1] - dims[0][0]))
                    ** 2
                    + ((Dshape[1] + 1) / 2 / Dshape[1] * abs(dims[1][-1] - dims[1][0]))
                    ** 2
                ),
            ]
        )
        for j in range(2):
            if maxl[j] < 0:
                maxl[j] = Dshape[j]
        D[0, 0] = 0.0
        try:
            N[0, 0] = arr[arr.mask == False].size / arr.size
        except:
            N[0, 0] = 1.0
        ell[0, 0] = 0.0
    if D3:
        lmax = np.min(
            [
                lmax,
                np.sqrt(
                    ((Dshape[0] + 1) / 2 / Dshape[0] * abs(dims[0][-1] - dims[0][0]))
                    ** 2
                    + ((Dshape[1] + 1) / 2 / Dshape[1] * abs(dims[1][-1] - dims[1][0]))
                    ** 2
                    + ((Dshape[2] + 1) / 2 / Dshape[2] * abs(dims[2][-1] - dims[2][0]))
                    ** 2
                ),
            ]
        )
        for j in range(3):
            if maxl[j] < 0:
                maxl[j] = Dshape[j]
        D[0, 0, 0] = 0.0
        try:
            N[0, 0, 0] = arr[arr.mask == False].size / arr.size
        except:
            N[0, 0, 0] = 1.0
        ell[0, 0, 0] = 0.0
    if not quiet:
        print(lmax)
    # compute correlations
    zskip = max(1, int(np.random.uniform() * loopsample[0]))
    for iz in range(1, Dshape[0], zskip):
        # downsample array in steps of zstep over subset starting at iz0
        zstep = max(np.mod(iz, downsample[0]), 1)
        iz0 = int(np.random.uniform() / zstep * (nz - 1))
        if not len(Dshape) > 1:
            if D1:
                marr = np.ma.array(arr[::zstep])
            elif D2:
                marr = np.ma.array(arr[::zstep, ::zstep])
            else:
                marr = np.ma.array(arr[::zstep, ::zstep, ::zstep])
            ell[iz] = np.sqrt((dims[Daxis[0]][iz] - dims[Daxis[0]][0]) ** 2)
            zshift = np.roll(marr, np.mod(zstep, iz), axis=Daxis[0])
            shift_diff = zshift[iz0 : iz0 + maxl[0]] - marr[iz0 : iz0 + maxl[0]]
            shift_power = np.power(shift_diff, Dorder)
            D[iz] = np.mean(shift_power[np.logical_not(np.isnan(shift_diff))])
            if D1:
                N[iz] = (
                    marr[iz0 : iz0 + maxl[0]][
                        marr[iz0 : iz0 + maxl[0]].mask == False
                    ].size
                    / arr[::zstep][iz0 : iz0 + maxl[0]].size
                )
            elif D2:
                N[iz] = (
                    marr[iz0 : iz0 + maxl[0]][
                        marr[iz0 : iz0 + maxl[0]].mask == False
                    ].size
                    / arr[::zstep, ::ystep][iz0 : iz0 + maxl[0]].size
                )
            elif D3:
                N[iz] = (
                    marr[iz0 : iz0 + maxl[0]][
                        marr[iz0 : iz0 + maxl[0]].mask == False
                    ].size
                    / arr[::zstep, ::ystep, ::xstep][iz0 : iz0 + maxl[0]].size
                )
        else:
            marr = np.ma.array(arr[::zstep])
            zshift = np.roll(marr, int(iz / zstep), axis=Daxis[0])[iz0 : iz0 + maxl[0]]
            yskip = max(1, int(np.random.uniform() * loopsample[1]))
            for iy in range(0, Dshape[1], yskip):
                # downsample array in steps of ystep over subset starting at iy0
                ystep = max(np.mod(iy, downsample[1]), 1)
                iy0 = int(np.random.uniform() / ystep * (ny - 1))
                if not len(Dshape) > 2:
                    ell[iz, iy] = np.sqrt(
                        (dims[Daxis[0]][iz] - dims[Daxis[0]][0]) ** 2
                        + (dims[Daxis[1]][iy] - dims[Daxis[1]][0]) ** 2
                    )
                    # calculate D if ell <= lmax
                    if ell[iz, iy] <= lmax:
                        if D2:
                            marr = np.ma.array(arr[::zstep, ::ystep])
                        elif D3:
                            marr = np.ma.array(arr[::zstep, ::ystep, ::ystep])
                        yshift = np.roll(
                            zshift[:, ::ystep, ::ystep], int(iy / ystep), axis=Daxis[1]
                        )
                        # account for shear if D2 'yx'
                        if "x" in dirs and iy > 0 and not deltay == 0.0:
                            ishear = round(
                                2 * ny * deltay / (max(dims[0]) - min(dims[0])) / ystep
                            )
                            yshift[:, : int(iy / ystep)] = np.roll(
                                yshift[:, : int(iy / istep)], ishear, axis=0
                            )
                        shift_diff = (
                            yshift[:, iy0 : iy0 + maxl[1]]
                            - marr[iz0 : iz0 + maxl[0], iy0 : iy0 + maxl[1]]
                        )
                        shift_power = np.power(shift_diff, Dorder)
                        D[iz, iy] = np.ma.mean(
                            shift_power[np.logical_not(np.isnan(shift_diff))]
                        )
                        if D2:
                            N[iz, iy] = (
                                marr[iz0 : iz0 + maxl[0], iy0 : iy0 + maxl[1]][
                                    marr[iz0 : iz0 + maxl[0], iy0 : iy0 + maxl[1]].mask
                                    == False
                                ].size
                                / arr[::zstep, ::ystep][
                                    iz0 : iz0 + maxl[0], iy0 : iy0 + maxl[1]
                                ].size
                            )
                        elif D3:
                            N[iz, iy] = (
                                marr[iz0 : iz0 + maxl[0], iy0 : iy0 + maxl[1]][
                                    marr[iz0 : iz0 + maxl[0], iy0 : iy0 + maxl[1]].mask
                                    == False
                                ].size
                                / arr[::zstep, ::ystep, ::xstep][
                                    iz0 : iz0 + maxl[0], iy0 : iy0 + maxl[1]
                                ].size
                            )
                else:
                    marr = np.ma.array(arr[::zstep, ::ystep])
                    yshift = np.roll(
                        zshift[:, ::ystep], int(iy / ystep), axis=Daxis[1]
                    )[:, iy0 : iy0 + maxl[1]]
                    xskip = max(1, int(np.random.uniform() * loopsample[2]))
                    # print('iz {}, iy {}, zskip {} yskip, xskip {}'.format(iz,iy,zskip,yskip,zskip))
                    for ix in range(0, Dshape[2], xskip):
                        # downsample array in steps of xstep over subset starting at ix0
                        xstep = max(np.mod(ix, downsample[2]), 1)
                        ix0 = int(np.random.uniform() / xstep * (nx - 1))
                        ell[iz, iy, ix] = np.sqrt(
                            (dims[Daxis[0]][iz] - dims[Daxis[0]][0]) ** 2
                            + (dims[Daxis[1]][iy] - dims[Daxis[1]][0]) ** 2
                            + (dims[Daxis[2]][ix] - dims[Daxis[2]][0]) ** 2
                        )
                        # calculate D if ell <= lmax
                        if ell[iz, iy, ix] <= lmax:
                            marr = np.ma.array(arr[::zstep, ::ystep, ::xstep])
                            xshift = np.roll(
                                yshift[:, :, ::xstep], int(ix / xstep), axis=Daxis[2]
                            )
                            # account for shear if D3 'zyx'
                            if ix > 0 and not deltay == 0.0:
                                ishear = round(
                                    2
                                    * ny
                                    * deltay
                                    / (max(dims[1]) - min(dims[1]))
                                    / ystep
                                )
                                xshift[:, :, : int(ix / xstep)] = np.roll(
                                    xshift[:, :, : int(ix / xstep)], ishear, axis=1
                                )
                            shift_diff = (
                                xshift[:, :, ix0 : ix0 + maxl[2]]
                                - marr[
                                    iz0 : iz0 + maxl[0],
                                    iy0 : iy0 + maxl[1],
                                    ix0 : ix0 + maxl[2],
                                ]
                            )
                            shift_power = np.power(shift_diff, Dorder)
                            D[iz, iy, ix] = np.mean(
                                shift_power[np.logical_not(np.isnan(shift_diff))]
                            )
                            # N[iz,iy,ix] = marr[iz0:iz0+maxl[0],iy0:iy0+maxl[1],ix0:ix0+maxl[2]][marr[iz0:iz0+maxl[0],iy0:iy0+maxl[1],ix0:ix0+maxl[2]].mask==False].size/arr[::zstep,::ystep,::xstep][iz0:iz0+maxl[0],iy0:iy0+maxl[1],ix0:ix0+maxl[2]].size
                            N[iz, iy, ix] = (
                                shift_power[shift_power.mask == False].size
                                / shift_power.size
                            )  # marr[iz0:iz0+maxl[0],iy0:iy0+maxl[1],ix0:ix0+maxl[2]][marr[iz0:iz0+maxl[0],iy0:iy0+maxl[1],ix0:ix0+maxl[2]].mask==False].size/arr[::zstep,::ystep,::xstep][iz0:iz0+maxl[0],iy0:iy0+maxl[1],ix0:ix0+maxl[2]].size
    ltmp = np.unique(ell[np.logical_not(np.isnan(ell))])
    l1dtmp = np.linspace(0, lmax, int(lmax / ltmp[1]))
    print("l1dtmp.size", l1dtmp.size)
    l1d = [0]
    D1d = [0]
    if D1:
        N1d = [N[0]]
    if D2:
        N1d = [N[0, 0]]
    if D3:
        N1d = [N[0, 0, 0]]
    D1dsd = [0]
    for il in range(1, l1dtmp.size):
        Dtmp = np.ma.array(D.copy())[np.logical_not(np.isnan(D))]
        Ntmp = np.ma.array(N.copy())[np.logical_not(np.isnan(D))]
        elltmp = np.ma.array(ell.copy())[np.logical_not(np.isnan(D))]
        Dtmp[elltmp <= l1dtmp[il - 1]] = np.ma.masked
        Dtmp[elltmp > l1dtmp[il]] = np.ma.masked
        Ntmp[elltmp <= l1dtmp[il - 1]] = np.ma.masked
        Ntmp[elltmp > l1dtmp[il]] = np.ma.masked
        if Dtmp[Dtmp.mask == False].size > 0:
            tmp = np.ma.mean(Dtmp)
            Nmp = np.ma.mean(Ntmp)
            std = np.ma.std(Dtmp)
            if not quiet:
                print(il, tmp)
            if not np.isnan(tmp):
                D1d.append(tmp)
                N1d.append(Nmp)
                D1dsd.append(std)
                l1d.append(l1dtmp[il])
                if not quiet:
                    print(l1dtmp[il], tmp)

    D1d = np.array(D1d)
    N1d = np.array(N1d)
    s1d = np.array(D1dsd)
    l1d = np.array(l1d)
    print("size D1d {} and l1d {}".format(D1d.size, l1d.size))
    # perform curve fit
    if D1d.size > 1:
        popt, pcov = curve_fit(fit_gaussm1, l1d, D1d / N1d, InitialGuess)
    else:
        popt, pcov = [0, 1], [0, 1]
    print("Fit parameters of sigma {:.2e} and L0 {:.2e}".format(popt[0], popt[1]))

    print(
        "Calculated spatial structure function of "
        + "{:.02f} MB in {:.02f} min.".format(
            arrchunksize, (time.time() - start_time) / 60
        )
    )
    if lplot:
        if D1d.size > 1:
            plt.figure()
            plt.plot(l1d, D1d / N1d, "-", label=dlabel)
            plt.xlabel(r"$\ell$")
            plt.ylabel(r"$D(\ell)$")
            plt.plot(
                l1d,
                fit_gaussm1(l1d, *popt),
                ".",
                label=r"fit $\sigma={:.2e},L_0={:.2e}$".format(*popt),
            )
            plt.legend()
            if figname:
                plt.savefig(figname)
                plt.close()
            else:
                plt.show()

    return D1d, l1d, N1d, np.array(popt), D1dsd


# ------------------------------------------------------------------------------


def space_struct_mpi(
    arr,
    dims=[[], [], []],
    Dorder=2,
    dirs="zyx",
    lperi=(True, True, True),
    lplot=True,
    InitialGuess=[1.0, 1.0],
    deltay=0.0,
    chunksize=1000.0,
    figname=None,
    dlabel="data",
    quiet=True,
):
    """
    Clone from above -- not yet implemented
    Calculate the structure function in space.

    call signature:

    space_struct(arr, dims=[[],[],[]], Dorder=2, dirs='zyx',
                 lperi=(True,True,True), lplot=True,
                 InitialGuess = [1.0,1.0],
                 deltay=0.,
                 figname = None, dlabel = 'data',
                 quiet=False)

    Keyword arguments:

    *arr*:
      numpy data array of shape, e.g., [nz,ny,nz].

    *z*:
      z dimension full array to replace default empty list.

    *y*:
      y dimension full array to replace default empty list.

    *x*:
      x dimension full array to replace default empty list.

    *Dorder*:
      Order of structure function, by default 2nd.

    *lperi*:
      Flag indicates for each dimension, whether boundary is periodic.

    *lplot*
      Flag to plot the result.

    *InitialGuess*:
      Initial parameters for curve fitting, default [1.0,1.0].

    *figname*:
      String name for plot, if saving rather than only display.

    *dlabel*:
      legend label for data.

    *quiet*
      Flag for switching off output.

    Returns
    -------
    D1 Structure function, length array, fit parameters [sigma, L0].

    Notes
    -----
    Must provide list of dimension arrays to match dirs required.

    Examples
    --------
    >>> D,ell,fit = pc.math.stats.space_struct(var.rho, dims=[gd.z,gd.y,gd.x],
                                               dirs='zyx')
    >>> D, ell, fit
    [0, ...,], [0, 0.0041, ...], [2.5,0.04]
    """
    start_time = time.time()
    ldirmatch = True
    arrshape = arr.shape
    if not quiet:
        print("arrshape", arrshape)
    D3 = len(arrshape) == 3
    D2 = len(arrshape) == 2
    D1 = len(arrshape) == 1
    if not quiet:
        print("D3 {}, D2 {}, D1 {}".format(D3, D2, D1))
    if D3:
        arrchunksize = 8 * arrshape[0] * arrshape[1] * arrshape[2] / 1024 / 1024
    elif D2:
        arrchunksize = 8 * arrshape[0] * arrshape[1] / 1024 / 1024
    else:
        arrchunksize = 8 * arrshape[0] / 1024 / 1024
    # include which dimensions in function D?
    axis = 0
    Dshape = []
    Daxis = []
    if "z" in dirs:
        nz = int(len(dims[axis]) / 2)
        Dshape.append(nz)
        Daxis.append(axis)
        if not quiet:
            print("z in dirs {}: Dshape {}, Daxis {}".format(dirs, Dshape, Daxis))
        axis += 1
        ldirmatch = 2 * nz == arrshape[0]
        if not ldirmatch:
            print(
                "Please correct size dims(z) "
                + "{} does not equal arr shape {}".format(2 * nz, arrshape[0])
            )
            sys.exit()
    else:
        nz = 0
    if "y" in dirs:
        if D3 and "z" not in dirs:
            axis += 1
        ny = int(len(dims[axis]) / 2)
        Dshape.append(ny)
        Daxis.append(axis)
        if not quiet:
            print("y in dirs {}: Dshape {}, Daxis {}".format(dirs, Dshape, Daxis))
        axis += 1
        ldirmatch = 2 * ny == arrshape[axis]
        if not ldirmatch:
            print(
                "Please correct size dims(y) "
                + "{} does not equal arr shape {}".format(2 * ny, arrshape[axis])
            )
            sys.exit()
    else:
        ny = 0
    if "x" in dirs:
        if D3 and "y" not in dirs:
            if "z" not in dirs:
                axis += 2
            else:
                axis += 1
        if D2 and "y" not in dirs and "z" not in dirs:
            axis += 1
        nx = int(len(dims[axis]) / 2)
        Dshape.append(nx)
        Daxis.append(axis)
        if not quiet:
            print("x in dirs {}: Dshape {}, Daxis {}".format(dirs, Dshape, Daxis))
        ldirmatch = 2 * nx == arrshape[axis]
        if not ldirmatch:
            print(
                "Please correct size dims(x) "
                + "{} does not equal arr shape {}".format(2 * nx, arrshape[axis])
            )
            sys.exit()
    else:
        nx = 0
    if len(Dshape) == 0:
        print(
            "length arrays are empty lists; at least 1 dimension from"
            " z, y or x must be provided"
        )
    if len(Dshape) > len(arrshape):
        print(
            "N dirs {} cannot exceed N dims {} of arr".format(
                len(Dshape), len(arrshape)
            )
        )
    lmax = 0
    if D1:
        lmax = max(
            lmax, (Dshape[0] + 1) / 2 / Dshape[0] * abs(dims[0][-1] - dims[0][0])
        )
    if D2:
        lmax = np.max(
            [
                lmax,
                (Dshape[0] + 1) / 2 / Dshape[0] * abs(dims[0][-1] - dims[0][0]),
                (Dshape[1] + 1) / 2 / Dshape[1] * abs(dims[1][-1] - dims[1][0]),
            ]
        )
    if D3:
        lmax = np.max(
            [
                lmax,
                (Dshape[0] + 1) / 2 / Dshape[0] * abs(dims[0][-1] - dims[0][0]),
                (Dshape[1] + 1) / 2 / Dshape[1] * abs(dims[1][-1] - dims[1][0]),
                (Dshape[2] + 1) / 2 / Dshape[2] * abs(dims[2][-1] - dims[2][0]),
            ]
        )
    if not quiet:
        print(lmax)
    mar = np.ma.array(arr)  # permit handling of masked arrays
    # compute correlations
    try:
        from mpi4py import MPI

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        l_mpi = True
        l_mpi = l_mpi and (size != 1)
    except ImportError:
        rank = 0
        size = 1
        comm = None
        l_mpi = False
    Dtmp = np.zeros(Dshape)
    elltmp = np.zeros(Dshape)
    if l_mpi:
        iindz = np.array_split(np.arange(Dshape[0]), size)
        indz = iindz[rank]
        D = Dtmp[indz]
        ell = elltmp[indz]
        marr = mar[indz]
    else:
        indz = np.arange(Dshape[0])
        D = Dtmp[indz]
        ell = elltmp[indz]
        marr = mar[indz]
    print("rank {} D.shape {} ell.shape {}".format(rank, D.shape, ell.shape))
    bzshift = marr.copy()
    for iz in range(nz):
        print("rank {}, iz {}".format(rank, iz))
        fzshift = np.roll(bzshift, 1, axis=Daxis[0])
        if l_mpi:
            comm.send(
                bzshift[-np.mod(iz, indz.size) :],
                np.mod(rank - int(iz / indz.size), size),
                tag=rank,
            )
            fzshift[: np.mod(iz, indz.size)] = comm.recv(
                np.mod(rank + int(iz / indz.size), size),
                tag=np.mod(rank + int(iz / indz.size)),
            )
            print(
                "rank {} sent to {} data size {}".format(
                    rank, np.mpd(rank + int(iz / indz.size)), iz
                )
            )
            print(
                "rank {} received from to {} data size {}".format(
                    rank, np.mpd(rank - int(iz / indz.size)), iz
                )
            )
        bzshift = fzshift
        if not len(Dshape) > 1:
            ell[iz] = np.sqrt((dims[Daxis[0]][iz] - dims[Daxis[0]][0]) ** 2)
            shift_diff = zshift - marr
            shift_power = np.power(shift_diff, Dorder)
            D[iz] = np.mean(shift_power[np.logical_not(np.isnan(shift_diff))])
        else:
            for iy in range(0, Dshape[1]):
                yshift = np.roll(fzshift, iy, axis=Daxis[1])
                if not len(Dshape) > 2:
                    ell[iz, iy] = np.sqrt(
                        (dims[Daxis[0]][iz] - dims[Daxis[0]][0]) ** 2
                        + (dims[Daxis[1]][iy] - dims[Daxis[1]][0]) ** 2
                    )
                    if ell[iz, iy] <= lmax:
                        # account for shear if D2 'yx'
                        if "x" in dirs and iy > 0 and not deltay == 0.0:
                            ishear = round(
                                2 * ny * deltay / (max(dims[0]) - min(dims[0]))
                            )
                            yshift[:, :iy] = np.roll(yshift[:, :iy], ishear, axis=0)
                        shift_diff = yshift - marr
                        shift_power = np.power(shift_diff, Dorder)
                        D[iz, iy] = np.ma.mean(
                            shift_power[np.logical_not(np.isnan(shift_diff))]
                        )
                else:
                    for ix in range(0, Dshape[2]):
                        ell[iz, iy, ix] = np.sqrt(
                            (dims[Daxis[0]][iz] - dims[Daxis[0]][0]) ** 2
                            + (dims[Daxis[1]][iy] - dims[Daxis[1]][0]) ** 2
                            + (dims[Daxis[2]][ix] - dims[Daxis[2]][0]) ** 2
                        )
                        if ell[iz, iy, ix] <= lmax:
                            xshift = np.roll(yshift, ix, axis=Daxis[2])
                            # account for shear if D3 'zyx'
                            if ix > 0 and not deltay == 0.0:
                                ishear = round(
                                    2 * ny * deltay / (max(dims[1]) - min(dims[1]))
                                )
                                xshift[:, :, :ix] = np.roll(
                                    xshift[:, :, :ix], ishear, axis=1
                                )
                            shift_diff = xshift - marr
                            shift_power = np.power(shift_diff, Dorder)
                            D[iz, iy, ix] = np.mean(
                                shift_power[np.logical_not(np.isnan(shift_diff))]
                            )
    if l_mpi:
        for i in range(size):
            Dtmp[indz] = comm.bcast(D, root=i)
            elltmp[indz] = comm.bcast(ell, root=i)
    else:
        Dtmp[indz] = D
        elltmp[indz] = ell
    D = Dtmp.copy()
    ell = elltmp.copy()
    print("rank {} ell.shape {}".format(rank, ell.shape))
    ltmp = np.unique(ell)
    l1dtmp = np.linspace(0, lmax, int(lmax / ltmp[1]))
    print("l1dtmp.size", l1dtmp.size)
    l1d = [0]
    D1d = [0]
    D1dsd = [0]
    for il in range(1, l1dtmp.size):
        Dtmp = np.ma.array(D.copy())
        Dtmp[ell <= l1dtmp[il - 1]] = np.ma.masked
        Dtmp[ell > l1dtmp[il]] = np.ma.masked
        tmp = np.ma.mean(Dtmp)
        std = np.ma.std(Dtmp)
        if not quiet:
            print(il, tmp)
        if not np.isnan(tmp):
            D1d.append(tmp)
            D1dsd.append(std)
            l1d.append(l1dtmp[il])
            if not quiet:
                print(l1dtmp[il], tmp)

    D1d = np.array(D1d)
    s1d = np.array(D1dsd)
    l1d = np.array(l1d)
    print("size D1d {} and l1d {}".format(D1d.size, l1d.size))
    # perform curve fit
    popt, pcov = curve_fit(fit_gaussm1, l1d, D1d, InitialGuess)
    print("Fit parameters of sigma {:.2e} and L0 {:.2e}".format(popt[0], popt[1]))

    print(
        "Calculated spatial structure function of "
        + "{:.02f} MB in {:.02f} min.".format(
            arrchunksize, (time.time() - start_time) / 60
        )
    )
    if lplot:
        plt.figure()
        plt.plot(l1d, D1d, "-", label=dlabel)
        plt.xlabel(r"$\ell$")
        plt.ylabel(r"$D(\ell)$")
        plt.plot(
            l1d,
            fit_gaussm1(l1d, *popt),
            ".",
            label=r"fit $\sigma={:.2f},L_0={:.2f}$".format(*popt),
        )
        plt.legend()
        if figname:
            plt.savefig(figname)
            plt.close()
        else:
            plt.show()

    return D1d, l1d, np.array(popt), D1dsd


# ------------------------------------------------------------------------------
def time_struct(
    arr,
    time=[],
    order=2,
    lplot=True,
    InitialGuess=[1.0, 1.0],
    figname=None,
    dlabel="data",
    quiet=False,
):
    """
    Calculate the structure function in time.

    call signature:

    time_struct(arr, time=[], order=2, lplot=True, InitialGuess = [1.0,1.0],
            figname = None, dlabel = 'data',   quiet=False)

    Keyword arguments:

    *arr*:
      numpy data array of shape [ntime,:].

    *time*:
      time series for time correlations. By default empty list.

    *Dorder*:
      Order of structure function, by default 2nd.

    *lplot*
      Flag to plot the result.

    *InitialGuess*:
      Initial parameters for curve fitting, default [1.0,1.0].

    *figname*:
      String name for plot, if saving rather than only display.

    *dlabel*:
      legend label for data.

    *quiet*
      Flag for switching off output.

    Returns
    -------
    D1 Structure function, time array, fit parameters [sigma, L0].

    Notes
    -----
    Provide time series of D1 -- D3 sample sets of fixed spatial location.

    Examples
    --------
    >>> Dt,t,fit = pc.math.stats.space_struct(var.rho, dims=[gd.z,gd.y,gd.x],
                                               dirs='zyx')
    >>> Dt, t, fit
    [0, ...,], [0, 0.0041, ...], [2.5,0.04]
    """

    arrshape = arr.shape
    if np.size(time) > 0:
        ltime_corr = arrshape[0] == np.size(time)
        if not ltime_corr:
            print(
                "arr shape {} does not match time shape {}".format(
                    arrshape, np.size(time)
                )
            )
    nt = min(arrshape[0], np.size(time))
    Dt = list()
    Dtsd = list()
    t = list()
    # compute time correlations with fixed positions
    marr = np.ma.array(arr)  # permit handling of masked arrays
    for it in range(0, nt):
        tshift = np.roll(marr, it, axis=0)
        shift_diff = tshift - marr
        shift_order = np.power(shift_diff, Dorder)
        tmp = np.ma.mean(shift_order)
        if isinstance(tmp, float):
            Dt.append(tmp)
            Dtsd.append(np.std(shift_order))
            t.append(time[it])
    # perform curve fit
    popt, pcov = curve_fit(fit_expm1, t, Dt, InitialGuess)
    print(popt)

    if lplot:
        plt.figure()
        plt.plot(t, Dt, "-", label=dlabel)
        plt.xlabel(r"$t$")
        plt.ylabel(r"$D(t)$")
        plt.plot(
            t,
            fit_expm1(t, *popt),
            ".",
            label=r"fit $\sigma={:.2f},L_0={:.2f}$".format(*popt),
        )
        plt.legend()
        if figname:
            plt.savefig(figname)
            plt.close()
        else:
            plt.show()

    return np.array(Dt), np.array(t), np.array(popt), np.array(Dtsd)


##------------------------------------------------------------------------------
# def structure_function(*args, **kwargs):
#    """
#    Compute the set of structure functions for a simulation.
#
#    Signature:
#
#    structure_function(sim, snapshots=[], magic=['uu', 'bb', 'rho'],
#                       time=[], z=[], y=[], x=[], Dorder=2, quiet=False)
#
#    Parameters
#    ----------
#    *datadir*:  Directory where the data is stored.
#
#    *file_name*:  Filename to read.
#          If a filename is given, only that power spectrum is read.
#          By default it reads all the power spectrum files.
#
#    *quiet*:    Flag for switching off output.
#
#    Returns
#    -------
#    Class containing the different power spectrum as attributes.
#
#    Notes
#    -----
#    Use the attribute keys to get a list of attributes
#
#    Examples
#    --------
#    >>> D = pc.math.stat.structure_function()
#    >>> D.keys()
#    t
#    kin
#    krms
#    hel_kin
#    """
#
#    struct_tmp = Struct()
#    struct_tmp.calc(*args, **kwargs)
#    return struct_tmp
#
#
# class Struct(object):
#    """
#    Struct -- holds structure function data.
#    """
#
#    def __init__(self):
#        """
#        Fill members with default values.
#        """
#
#        self.t = []
#
#    def keys(self):
#        for i in self.__dict__.keys():
#            print(i)
#
#    def calc(self, var, snapshots=[], keys=['uu', 'rho'],
#                       time=[], z=[], y=[], x=[], Dorder=2,
#                 lperi=(True,True,True), lplot=True,
#                 InitialGuess = [1.0,1.0],
#                 deltay=0.,
#                 figname = None, dlabel = 'data',
#                 quiet=False):
#        for key in keys:
#            if not key in var.__dict__.keys():
#                print(key+' must be included in var object, skipping')
#            else:
#                if var.__getattribute__(key).shape[0]==3:
#                    arr = dot2(var.__getattribute__(key))
#                else:
#                    arr = var.__getattribute__(key).copy()
#            self.t = var.t
#
