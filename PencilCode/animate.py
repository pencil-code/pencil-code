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
    # Get the dimensions and the parameters.
    dim = read.dimensions(datadir=datadir)
    par = read.parameters(datadir=datadir)
    # Dispatch.
    ndim = (dim.nxgrid > 1) + (dim.nygrid > 1) + (dim.nzgrid > 1)
    if ndim == 1:
        raise NotImplementedError("1D run")
    elif ndim == 2:
        _slices_2D(field, t, s, dim, par, drange, **kwarg)
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
def _rectangular_2D(t, maps, extent, vmin, vmax, xlabel=None, ylabel=None, clabel=None, **kwarg):
    """Animates the evolution of a 2D rectangular plane.

    Positional Arguments:
        t
            Array of time points.
        maps
            Sequence in time of the data maps.
        extent
            Four-element arrays for the lower limit and the upper limit
            horizontally and then those vertically.
        vmin
            Scalar or array of same size as t; minimum data value(s).
        vmax
            Scalar or array of same size as t; maximum data value(s).

    Keyword Arguments:
        xlabel
            Label for the x axis.
        ylabel
            Label for the y axis.
        clabel
            Label for the color bar.
        **kwarg
            Keyword arguments passed to matplotlib.pyplot.figure().
    """
    # Chao-Chin Yang, 2015-04-25
    from collections.abc import Sequence
    import matplotlib.pyplot as plt
    import numpy as np
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
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("t = {:#.3G}".format(t[0]))
    im = ax.imshow(maps[0].transpose(), origin='lower', extent=extent, vmin=vmin0, vmax=vmax0, interpolation='nearest')
    cb = plt.colorbar(im)
    cb.set_label(clabel)
    plt.show(block=False)
    # Loop over each time and update the image.
    for i in range(1,len(t)):
        ax.set_title("t = {:#.3G}".format(t[i]))
        im.set_data(maps[i].transpose())
        if vmin_dynamic: im.set_clim(vmin=vmin[i])
        if vmax_dynamic: im.set_clim(vmax=vmax[i])
        fig.canvas.draw()
#=======================================================================
def _slices_2D(field, t, slices, dim, par, drange, **kwarg):
    """Animates video slices from a 2D model.

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
        drange
            Type of data range to be plotted; see _get_range().

    Keyword Arguments:
        **kwarg
            Keyword arguments passed to matplotlib.pyplot.figure().
    """
    # Chao-Chin Yang, 2015-04-23
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
    # Check the coordinate system.
    if par.coord_system != 'cartesian':
        raise NotImplementedError("2D, curvilinear model")
    xlabel, ylabel = "xyz"[xdir], "xyz"[ydir]
    # Irregular grid not implemented.
    if not (par.lequidist[xdir] and par.lequidist[ydir]):
        raise NotImplementedError("Irregular grid")
    maps = slices[:][plane]
    # Get the dimensions.
    xmin, xmax = par.xyz0[xdir], par.xyz1[xdir]
    ymin, ymax = par.xyz0[ydir], par.xyz1[ydir]
    extent = (xmin, xmax, ymin, ymax)
    # Get the data range.
    vmin, vmax = _get_range(t, maps, drange)
    # Send to the animator.
    _rectangular_2D(t, maps, extent, vmin, vmax, xlabel=xlabel, ylabel=ylabel, clabel=field, **kwarg)
#=======================================================================
