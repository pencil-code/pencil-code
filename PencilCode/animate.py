#!/usr/bin/python3
# Last Modification: $Id$
#=======================================================================
# animate.py
#
# Facilities for animating the Pencil Code data.
#
# Chao-Chin Yang, 2015-04-22
#=======================================================================
def avg2d(name, direction, datadir='./data', drange='full', **kwarg):
    """Animates the time sequence of a line average.

    Positional Argument:
        name
            Name of the average.
        direction
            Direction of the average: 'x', 'y', or 'z'.

    Keyword Arguments:
        datadir
            Path to the data directory.
        drange
            Type of data range to be plotted; see _get_range().
        **kwarg
            Other keyword arguments passed on.
    """
    # Chao-Chin Yang, 2015-04-30
    from . import read
    # Check the direction of the line average.
    if direction == 'x':
        xdir, ydir = 1, 2
    elif direction == 'y':
        xdir, ydir = 0, 2
    elif direction == 'z':
        xdir, ydir = 0, 1
    else:
        raise ValueError("Unknown direction '{}'. ".format(direction))
    # Get the dimensions, parameters, and grid.
    dim = read.dimensions(datadir=datadir)
    par = read.parameters(datadir=datadir)
    grid = read.grid(datadir=datadir, interface=True, par=par)
    # Read the averages.
    t, avg = read.avg2d(datadir=datadir, direction=direction)
    # Check the coordinate system.
    if par.coord_system != 'cartesian':
        raise NotImplementedError("curvilinear model")
    xlabel, ylabel = "xyz"[xdir], "xyz"[ydir]
    x, y = getattr(grid, xlabel), getattr(grid, ylabel)
    # Animate.
    _frame_rectangle(t, x, y, avg[name], drange, xlabel=xlabel, ylabel=ylabel, clabel=name, **kwarg)
#=======================================================================
def slices(field, datadir='./data', drange='full', **kwarg):
    """Dispatches to the respective animator of video slices.

    Positional Argument:
        field
            Field name.

    Keyword Arguments:
        datadir
            Path to the data directory.
        drange
            Style of the data range; see _get_range().
        **kwarg
            Keyword arguments passed to matplotlib.pyplot.figure().
    """
    # Chao-Chin Yang, 2015-04-29
    from . import read
    # Read the slices.
    t, s = read.slices(field, datadir=datadir)
    # Get the dimensions, parameters, and grid.
    dim = read.dimensions(datadir=datadir)
    par = read.parameters(datadir=datadir)
    grid = read.grid(datadir=datadir, interface=True, par=par)
    # Dispatch.
    ndim = (dim.nxgrid > 1) + (dim.nygrid > 1) + (dim.nzgrid > 1)
    if ndim == 1:
        raise NotImplementedError("1D run")
    elif ndim == 2:
        _slices2d(field, t, s, dim, par, grid, drange, **kwarg)
    elif ndim == 3:
        _slices3d(field, t, s, dim, par, grid, drange, **kwarg)
#=======================================================================
###### Local Functions ######
#=======================================================================
def _frame_rectangle(t, x, y, c, drange, xlabel=None, ylabel=None, clabel=None, **kwarg):
    """Animates a sequence of two-dimensional rectangular data.

    Positional Arguments:
        t
            1D array of the time points.
        x
            1D array of the x-coordinates of the cell interfaces.
        y
            1D array of the y-coordinates of the cell interfaces.
        c
            3D array for the sequence in time of the data.  c[i,j,k] is
            the value at time t[i] and cell with corners at
            (x[j-1],y[k]), (x[j+1],y[k]), (x[j],y[k-1]), and
            (x[j],y[k+1]).
        drange
            Style of the data range; see _get_range().

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
    # Chao-Chin Yang, 2015-04-30
    from collections.abc import Sequence
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    # Get the data range.
    vmin, vmax = _get_range(t, c, drange)
    seq = lambda a: isinstance(a, Sequence) or isinstance(a, np.ndarray)
    vmin_dynamic = seq(vmin)
    vmax_dynamic = seq(vmax)
    vmin0 = vmin[0] if vmin_dynamic else vmin
    vmax0 = vmax[0] if vmax_dynamic else vmax
    # Create the first plot.
    fig = plt.figure(**kwarg)
    ax = fig.gca()
    pc = ax.pcolorfast(x, y, c[0].transpose(), vmin=vmin0, vmax=vmax0)
    ax.minorticks_on()
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(y[0], y[-1])
    ax.set_aspect('equal')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("t = {:#.4G}".format(t[0]))
    cb = plt.colorbar(pc)
    cb.set_label(clabel)
    plt.show(block=False)
    # Loop over each time and update the plot.
    for i in range(1,len(t)):
        ax.set_title("t = {:#.4G}".format(t[i]))
        if mpl.__version__ == '1.3.1':
            pc.set_data(x, y, c[i].transpose())
        else:
            pc.set_data(c[i].transpose())
        if vmin_dynamic: pc.set_clim(vmin=vmin[i])
        if vmax_dynamic: pc.set_clim(vmax=vmax[i])
        fig.canvas.draw()
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
def _slices2d(field, t, slices, dim, par, grid, drange, **kwarg):
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
        grid
            Grid coordinates supplied by read.grid(interface=True).
        drange
            Style of the data range; see _get_range().

    Keyword Arguments:
        **kwarg
            Keyword arguments passed to matplotlib.pyplot.figure().
    """
    # Chao-Chin Yang, 2015-04-30
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
    # Check the coordinate system.
    if par.coord_system != 'cartesian':
        raise NotImplementedError("2D, curvilinear model")
    xlabel, ylabel = "xyz"[xdir], "xyz"[ydir]
    x, y = getattr(grid, xlabel), getattr(grid, ylabel)
    # Send to the animator.
    _frame_rectangle(t, x, y, slices[:][plane], drange, xlabel=xlabel, ylabel=ylabel, clabel=field, **kwarg)
#=======================================================================
def _slices3d(field, t, slices, dim, par, grid, drange, **kwarg):
    """Animates video slices from a 3D model.

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
        grid
            Grid coordinates supplied by read.grid(trim=True).
        drange
            Type of data range to be plotted; see _get_range().

    Keyword Arguments:
        **kwarg
            Keyword arguments passed to matplotlib.pyplot.figure().
    """
    # Chao-Chin Yang, 2015-04-27
    from collections.abc import Sequence
    from matplotlib.animation import FuncAnimation
    from matplotlib import cm
    from matplotlib.colors import Normalize
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np
    # Check the coordinate system.
    if par.coord_system != 'cartesian':
        raise NotImplementedError("3D, curvilinear model")
    print("Processing the data...")
    # Set the grid.
    g = lambda x, dx_1, x0, x1: np.concatenate(((x0,), (x[:-1] * dx_1[:-1] + x[1:] * dx_1[1:]) / (dx_1[:-1] + dx_1[1:]), (x1,)))
    x = g(grid.x, grid.dx_1, par.xyz0[0], par.xyz1[0])
    y = g(grid.y, grid.dy_1, par.xyz0[1], par.xyz1[1])
    z = g(grid.z, grid.dz_1, par.xyz0[2], par.xyz1[2])
    # Set the surfaces.
    nt = len(t)
    dtype = lambda nx, ny: [('value', np.float, (nt,nx,ny)),
                            ('x', np.float, (nx+1,ny+1)), ('y', np.float, (nx+1,ny+1)), ('z', np.float, (nx+1,ny+1))]
    planes = slices.dtype.names
    surfaces = []
    if 'xy' in planes or 'xy2' in planes:
        xmesh, ymesh = np.meshgrid(x, y, indexing='ij')
        if 'xy' in planes:
            zmesh = np.full(xmesh.shape, par.xyz0[2] - 0.7 * par.lxyz[2])
            surfaces.append(np.rec.array((slices.xy, xmesh, ymesh, zmesh), dtype=dtype(dim.nxgrid,dim.nygrid)))
        if 'xy2' in planes:
            zmesh = np.full(xmesh.shape, par.xyz1[2])
            surfaces.append(np.rec.array((slices.xy2, xmesh, ymesh, zmesh), dtype=dtype(dim.nxgrid,dim.nygrid)))
    if 'xz' in planes:
        xmesh, zmesh = np.meshgrid(x, z, indexing='ij')
        ymesh = np.full(xmesh.shape, par.xyz0[1])
        surfaces.append(np.rec.array((slices.xz, xmesh, ymesh, zmesh), dtype=dtype(dim.nxgrid,dim.nzgrid)))
    if 'yz' in planes:
        ymesh, zmesh = np.meshgrid(y, z, indexing='ij')
        xmesh = np.full(ymesh.shape, par.xyz0[0])
        surfaces.append(np.rec.array((slices.yz, xmesh, ymesh, zmesh), dtype=dtype(dim.nygrid,dim.nzgrid)))
    # Check the data range.
    vmin, vmax = _get_range(t, slices, drange)
    seq = lambda a: isinstance(a, Sequence) or isinstance(a, np.ndarray)
    vmin_dynamic = seq(vmin)
    vmax_dynamic = seq(vmax)
    vmin0 = vmin[0] if vmin_dynamic else vmin
    vmax0 = vmax[0] if vmax_dynamic else vmax
    # Set up the 3D view.
    print("Initializing...")
    fig = plt.figure(**kwarg)
    ax = Axes3D(fig)
    ax.view_init(15, -135)
    ax.dist = 15
    ax.set_xlim(par.xyz0[0], par.xyz1[0])
    ax.set_ylim(par.xyz0[1], par.xyz1[1])
    ax.set_zlim(par.xyz0[2], par.xyz1[2])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    # Plot the first surfaces.
    cmap = cm.get_cmap()
    norm = Normalize(vmin=vmin0, vmax=vmax0, clip=True)
    cols = [ax.plot_surface(s.x, s.y, s.z, facecolors=cmap(norm(s.value[0])), linewidth=0, shade=False)
            for s in surfaces]
    mappable = cm.ScalarMappable(cmap=cmap, norm=norm)
    mappable.set_array([])
    cb = fig.colorbar(mappable, shrink=0.8, aspect=20)
    cb.set_label(field)
    timestamp = ax.set_title("t = {:#.4G}".format(t[0]))
    # Animate the sequence.
    print("Animating...")
    def update(num, surfaces):
        nonlocal cols
        oldcols = cols
        cols = [ax.plot_surface(s.x, s.y, s.z, facecolors=cmap(norm(s.value[num])), linewidth=0, shade=False)
                for s in surfaces]
        for c in oldcols:
            ax.collections.remove(c)
        if vmin_dynamic or vmax_dynamic:
            if vmin_dynamic:
                norm.vmin = vmin[num]
            if vmax_dynamic:
                norm.vmax = vmax[num]
            mappable.set_norm(norm)
        timestamp.set_text("t = {:#.4G}".format(t[num]))
        return cols
    anim = FuncAnimation(fig, update, nt, fargs=(surfaces,), interval=1, blit=False, repeat=False)
    #anim.save('slices3d.mp4', writer='mencoder')
    plt.show()
