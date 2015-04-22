#!/usr/bin/python3
# Last Modification: $Id$
#=======================================================================
# animate.py
#
# Facilities for animating the Pencil Code data.
#
# Chao-Chin Yang, 2015-04-22
#=======================================================================
def slices(field, datadir='./data', **kwarg):
    """Dispatches to the respective animator of video slices.

    Positional Argument:
        field
            Field name.

    Keyword Arguments:
        datadir
            Data directory.
        **kwarg
            Keyword arguments passed to matplotlib.pyplot.figure().
    """
    # Chao-Chin Yang, 2015-04-22
    from . import read
    # Read the slices.
    t, s = read.slices(field, datadir=datadir)
    # Get the dimensions, parameters, and grid.
    dim = read.dimensions(datadir=datadir)
    par = read.parameters(datadir=datadir)
    # Dispatch.
    ndim = (dim.nxgrid > 1) + (dim.nygrid > 1) + (dim.nzgrid > 1)
    if ndim == 1:
        raise NotImplementedError("1D run")
    elif ndim == 2:
        if par.coord_system != 'cartesian':
            raise NotImplementedError("2D, curvilinear model")
        _slices_2D_rectangular(field, t, s, dim, par, **kwarg)
    elif ndim == 3:
        raise NotImplementedError("3D run")
#=======================================================================
###### Local Functions ######
#=======================================================================
def _slices_2D_rectangular(field, t, slices, dim, par, **kwarg):
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

    Keyword Arguments:
        **kwarg
            Keyword arguments passed to matplotlib.pyplot.figure().
    """
    # Chao-Chin Yang, 2015-04-22
    import matplotlib.pyplot as plt
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
    # Check the data range.
    vmin, vmax = slices[plane].min(), slices[plane].max()
    # Create the first image.
    fig = plt.figure(**kwarg)
    ax = fig.gca()
    ax.minorticks_on()
    ax.set_xlabel("xyz"[xdir])
    ax.set_ylabel("xyz"[ydir])
    im = ax.imshow(slices[0][plane].transpose(), origin='bottom', extent=(xmin, xmax, ymin, ymax), vmin=vmin, vmax=vmax,
                   interpolation='nearest')
    cb = plt.colorbar(im)
    cb.set_label(field)
    plt.show(block=False)
    # Loop over each slice and update the image.
    for i in range(1,len(t)):
        ax.set_title("t = {:#.3G}".format(t[i]))
        im.set_data(slices[i][plane].transpose())
        fig.canvas.draw()
