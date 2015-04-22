#!/usr/bin/python3
# Last Modification: $Id$
#=======================================================================
# animate.py
#
# Facilities for animating the Pencil Code data.
#
# Chao-Chin Yang, 2015-04-22
#=======================================================================
def slices(field, datadir='./data', drange='full', **kwarg):
    """Dispatches to the respective animator of video slices.

    Positional Argument:
        field
            Field name.

    Keyword Arguments:
        datadir
            Data directory.
        drange
            Type of data range to be plotted; see _get_range().
        **kwarg
            Keyword arguments passed to matplotlib.pyplot.figure().
    """
    # Chao-Chin Yang, 2015-04-23
    from . import read
    # Read the slices.
    t, s = read.slices(field, datadir=datadir)
    # Get the dimensions, parameters, and grid.
    dim = read.dimensions(datadir=datadir)
    par = read.parameters(datadir=datadir)
    # Determine the data range.
    vmin, vmax = _get_range(t, s, drange)
    # Dispatch.
    ndim = (dim.nxgrid > 1) + (dim.nygrid > 1) + (dim.nzgrid > 1)
    if ndim == 1:
        raise NotImplementedError("1D run")
    elif ndim == 2:
        if par.coord_system != 'cartesian':
            raise NotImplementedError("2D, curvilinear model")
        _slices_2D_rectangular(field, t, s, dim, par, vmin, vmax, **kwarg)
    elif ndim == 3:
        raise NotImplementedError("3D run")
#=======================================================================
###### Local Functions ######
#=======================================================================
def _get_range(t, data, drange):
    """Determines the data range to be plotted.

    Positional Arguments:
        t
            Sequence of time points.
        data
            Sequence in time of numpy arrays.
        drange
            Type of data range:
                'full'
                    Absolute minimum and maximum.
                'dynamic'
                    Dynamic minimum and maximum at each time.
            Otherwise, user-defined range and returned as is.
    """
    # Chao-Chin Yang, 2015-04-22
    from collections.abc import Sequence
    import numpy as np
    # User-defined range.
    if not isinstance(drange, str):
        return drange
    # Find the data range at each time.
    nt = len(t)
    vmin, vmax = np.empty(nt), np.empty(nt)
    for i, d in enumerate(data):
        if d.dtype.names is None:
            vmin1, vmax1 = d.min(), d.max()
        else:
            vmin1, vmax1 = np.inf, -np.inf
            for n in d.dtype.names:
                vmin1 = min(vmin1, d[n].min())
                vmax1 = max(vmax1, d[n].max())
        vmin[i], vmax[i] = vmin1, vmax1
    # Check the type of range requested.
    if drange == 'dynamic':
        pass
    elif drange == 'full':
        vmin, vmax = vmin.min(), vmax.max()
    else:
        raise ValueError("Unknown type of range '{}'".format(drange))
    return vmin, vmax
#=======================================================================
def _slices_2D_rectangular(field, t, slices, dim, par, vmin, vmax, **kwarg):
    """Animates video slices from a 2D, rectangular model.

    Positional Arguments:
        field
            Field name.
        t
            Array of time points.
        slices
            Record array of video slices supplied by read.slices().
        dim
            Dimensions supplied by read.dimensions().
        par
            Parameters supplied by read.parameters().
        vmin
            Scalar or array of the same size as t; minimum data value.
        vmax
            Scalar or array of the same size as t; maximum data value.

    Keyword Arguments:
        **kwarg
            Keyword arguments passed to matplotlib.pyplot.figure().
    """
    # Chao-Chin Yang, 2015-04-23
    from collections.abc import Sequence
    import matplotlib.pyplot as plt
    import numpy as np
    # Determine the slice plane.
    if dim.nxgrid == 1:
        plane = 'yz'
        xdir, ydir = 1, 2
    elif dim.nygrid == 1:
        plane = 'xz'
        xdir, ydir = 0, 2
    else:
        plane = 'xy'
        xdir, ydir = 0, 1
    xmin, xmax = par.xyz0[xdir], par.xyz1[xdir]
    ymin, ymax = par.xyz0[ydir], par.xyz1[ydir]
    if not (par.lequidist[xdir] and par.lequidist[ydir]):
        raise NotImplementedError("Irregular grid")
    # Check the specified range.
    seq = lambda a: isinstance(a, Sequence) or isinstance(a, np.ndarray)
    vmin_dynamic = seq(vmin)
    vmax_dynamic = seq(vmax)
    vmin0 = vmin[0] if vmin_dynamic else vmin
    vmax0 = vmax[0] if vmax_dynamic else vmax
    # Create the first image.
    fig = plt.figure(**kwarg)
    ax = fig.gca()
    ax.minorticks_on()
    ax.set_xlabel("xyz"[xdir])
    ax.set_ylabel("xyz"[ydir])
    im = ax.imshow(slices[0][plane].transpose(), origin='lower', extent=(xmin, xmax, ymin, ymax), vmin=vmin0, vmax=vmax0,
                   interpolation='nearest')
    cb = plt.colorbar(im)
    cb.set_label(field)
    plt.show(block=False)
    # Loop over each slice and update the image.
    for i in range(1,len(t)):
        ax.set_title("t = {:#.3G}".format(t[i]))
        im.set_data(slices[i][plane].transpose())
        if vmin_dynamic: im.set_clim(vmin=vmin[i])
        if vmax_dynamic: im.set_clim(vmax=vmax[i])
        fig.canvas.draw()
