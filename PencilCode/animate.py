#!/usr/bin/python3
# Last Modification: $Id$
#=======================================================================
# animate.py
#
# Facilities for animating the Pencil Code data.
#
# Chao-Chin Yang, 2015-04-22
#=======================================================================
def avg2d(name, direction, datadir='./data', **kwarg):
    """Animates the time sequence of a line average.

    Positional Argument:
        name
            Name of the average.
        direction
            Direction of the average: 'x', 'y', or 'z'.

    Keyword Arguments:
        datadir
            Path to the data directory.
        **kwarg
            Keyword arguments passed to _get_range().
    """
    # Chao-Chin Yang, 2015-05-04
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
    _frame_rectangle(t, x, y, avg[name], xlabel=xlabel, ylabel=ylabel, clabel=name, **kwarg)
#=======================================================================
def slices(field, datadir='./data', **kwarg):
    """Dispatches to the respective animator of video slices.

    Positional Argument:
        field
            Field name.

    Keyword Arguments:
        datadir
            Path to the data directory.
        **kwarg
            Keyword arguments passed to _get_range().
    """
    # Chao-Chin Yang, 2015-05-04
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
        _slices2d(field, t, s, dim, par, grid, **kwarg)
    elif ndim == 3:
        _slices3d(field, t, s, dim, par, grid, **kwarg)
#=======================================================================
###### Local Functions ######
#=======================================================================
def _frame_rectangle(t, x, y, c, xlabel=None, ylabel=None, clabel=None, **kwarg):
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

    Keyword Arguments:
        xlabel
            Label for the x axis.
        ylabel
            Label for the y axis.
        clabel
            Label for the color bar.
        **kwarg
            Keyword arguments passed to _get_range().
    """
    # Chao-Chin Yang, 2015-05-26
    from collections.abc import Sequence
    from matplotlib.colors import LogNorm, Normalize
    import matplotlib.pyplot as plt
    import numpy as np
    # Get the data range.
    vmin, vmax = _get_range(t, c, **kwarg)
    logscale = kwarg.pop("logscale", False)
    seq = lambda a: isinstance(a, Sequence) or isinstance(a, np.ndarray)
    vmin_dynamic = seq(vmin)
    vmax_dynamic = seq(vmax)
    vmin0 = vmin[0] if vmin_dynamic else vmin
    vmax0 = vmax[0] if vmax_dynamic else vmax
    norm = LogNorm(vmin=vmin0, vmax=vmax0) if logscale else Normalize(vmin=vmin0, vmax=vmax0)
    if logscale:
        c = c.clip(min(vmin) if vmin_dynamic else vmin, np.inf)
    # Create the first plot.
    fig = plt.figure()
    ax = fig.gca()
    pc = ax.pcolormesh(x, y, c[0].transpose(), norm=norm)
    ax.minorticks_on()
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(y[0], y[-1])
    ax.set_aspect('equal')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("$t = {:#.4G}$".format(t[0]))
    cb = plt.colorbar(pc)
    cb.set_label(clabel)
    plt.show(block=False)
    # Loop over each time and update the plot.
    for i in range(1,len(t)):
        ax.set_title("$t = {:#.4G}$".format(t[i]))
        pc.set_array(c[i].ravel(order='F'))
        if vmin_dynamic and vmax_dynamic:
            pc.set_clim(vmin[i], vmax[i])
        elif vmin_dynamic:
            pc.set_clim(vmin=vmin[i])
        elif vmax_dynamic:
            pc.set_clim(vmax=vmax[i])
        fig.canvas.draw()
#=======================================================================
def _get_range(t, data, center=False, drange='full', logscale=False, tmin=None):
    """Determines the data range to be plotted.

    Positional Arguments:
        t
            Sequence of time points.
        data
            Sequence in time of numpy arrays.
        center
            If center is True, the range is centered at 0 for
            logscale = False, or 1 for logscale = True.  If center is
            a number, the range is centered at this number.
        drange
            Type of data range:
                'dynamic'
                    Dynamic minimum and maximum at each time.
                'full'
                    Absolute minimum and maximum.
                'mean'
                    Time-averaged minimum and maximum.
            Otherwise, user-defined range and returned as is.
        tmin
            If not None, the range determination is restricted to
            t >= tmin.  No effect if drange is 'dynamic'.
        logscale
            Whether or not the color map is in logarithmic scale.
    """
    # Chao-Chin Yang, 2015-05-05
    from collections.abc import Sequence
    import numpy as np
    from scipy.integrate import simps
    from sys import float_info
    # User-defined range.
    if not isinstance(drange, str):
        return drange
    # Find the data range at each time.
    nt = len(t)
    vmin, vmax, vposmin = np.empty(nt), np.empty(nt), np.empty(nt)
    tiny = float_info.min
    for i, d in enumerate(data):
        if d.dtype.names is None:
            vmin1, vmax1 = d.min(), d.max()
            indices = d > 0
            vposmin1 = d[indices].min() if indices.any() else np.inf
        else:
            vmin1, vmax1, vposmin1 = np.inf, -np.inf, np.inf
            for n in d.dtype.names:
                d1 = d[n]
                vmin1 = min(vmin1, d1.min())
                vmax1 = max(vmax1, d1.max())
                indices = d1 > 0
                vposmin1 = min(vposmin1, (d1[indices].min() if indices.any() else np.inf))
        vmin[i], vmax[i] = vmin1, vmax1
        vposmin[i] = vposmin1 if vposmin1 is not np.inf else tiny
    # Check the type of range requested.
    if drange == 'dynamic':
        pass
    else:
        if tmin is not None:
            indices = t >= tmin
            t, vmin, vmax, vposmin = t[indices], vmin[indices], vmax[indices], vposmin[indices]
            nt = len(t)
        if drange == 'full':
            vmin, vmax = vmin.min(), vmax.max()
            vposmin = vposmin.min()
        elif drange == 'mean':
            dt = t[-1] - t[0]
            vmin, vmax = simps(vmin, t) / dt, simps(vmax, t) / dt
            vposmin = vposmin.min()
        else:
            raise ValueError("Unknown type of range '{}'".format(drange))
    # Guard against non-positive number if logscale is True.
    if logscale:
        vmin = max(vmin, vposmin)
    # Center the range if requested.
    if isinstance(center, bool) and center or not isinstance(center, bool) and isinstance(int, float):
        # Function to center the range.
        c = (1 if logscale else 0) if isinstance(center, bool) else center
        def get_centered_range(vmin, vmax):
            if logscale:
                a, b = c / vmin, vmax / c
                if a > 1 and b > 1:
                    if a > b:
                        vmin = c / b
                    elif a < b:
                        vmax = c * a
            else:
                a, b = c - vmin, vmax - c
                if a > 0 and b > 0:
                    if a > b:
                        vmin = c - b
                    elif a < b:
                        vmax = c + a
            return vmin, vmax
        # Update the range.
        if type(vmin) is np.ndarray:
            for i in range(nt):
                vmin[i], vmax[i] = get_centered_range(vmin[i], vmax[i])
        else:
            vmin, vmax = get_centered_range(vmin, vmax)
    return vmin, vmax
#=======================================================================
def _slices2d(field, t, slices, dim, par, grid, **kwarg):
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

    Keyword Arguments:
        **kwarg
            Keyword arguments passed to matplotlib.pyplot.figure().
    """
    # Chao-Chin Yang, 2015-05-27
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
    # Replace the labels if in **kwarg.
    xlabel = kwarg.pop("xlabel", xlabel)
    ylabel = kwarg.pop("ylabel", ylabel)
    clabel = kwarg.pop("clabel", field)
    # Send to the animator.
    _frame_rectangle(t, x, y, slices[:][plane], xlabel=xlabel, ylabel=ylabel, clabel=clabel, **kwarg)
#=======================================================================
def _slices3d(field, t, slices, dim, par, grid, **kwarg):
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

    Keyword Arguments:
        **kwarg
            Keyword arguments passed to matplotlib.pyplot.figure().
    """
    # Chao-Chin Yang, 2015-05-05
    from collections.abc import Sequence
    from matplotlib.animation import FuncAnimation
    from matplotlib import cm
    from matplotlib.colors import LogNorm, Normalize
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np
    # Check the coordinate system.
    if par.coord_system != 'cartesian':
        raise NotImplementedError("3D, curvilinear model")
    print("Processing the data...")
    # Get the data range.
    planes = slices.dtype.names
    vmin, vmax = _get_range(t, slices, **kwarg)
    logscale = kwarg.pop("logscale", False)
    seq = lambda a: isinstance(a, Sequence) or isinstance(a, np.ndarray)
    vmin_dynamic = seq(vmin)
    vmax_dynamic = seq(vmax)
    vmin0 = vmin[0] if vmin_dynamic else vmin
    vmax0 = vmax[0] if vmax_dynamic else vmax
    norm = LogNorm(vmin=vmin0, vmax=vmax0) if logscale else Normalize(vmin=vmin0, vmax=vmax0)
    if logscale:
        a = min(vmin) if vmin_dynamic else vmin
        for p in planes:
            slices[p] = slices[p].clip(a, np.inf)
    # Set the surfaces.
    nt = len(t)
    dtype = lambda nx, ny: [('value', np.float, (nt,nx,ny)),
                            ('x', np.float, (nx+1,ny+1)), ('y', np.float, (nx+1,ny+1)), ('z', np.float, (nx+1,ny+1))]
    surfaces = []
    if 'xy' in planes or 'xy2' in planes:
        xmesh, ymesh = np.meshgrid(grid.x, grid.y, indexing='ij')
        if 'xy' in planes:
            zmesh = np.full(xmesh.shape, par.xyz0[2] - 0.7 * par.lxyz[2])
            surfaces.append(np.rec.array((slices.xy, xmesh, ymesh, zmesh), dtype=dtype(dim.nxgrid,dim.nygrid)))
        if 'xy2' in planes:
            zmesh = np.full(xmesh.shape, par.xyz1[2])
            surfaces.append(np.rec.array((slices.xy2, xmesh, ymesh, zmesh), dtype=dtype(dim.nxgrid,dim.nygrid)))
    if 'xz' in planes:
        xmesh, zmesh = np.meshgrid(grid.x, grid.z, indexing='ij')
        ymesh = np.full(xmesh.shape, par.xyz0[1])
        surfaces.append(np.rec.array((slices.xz, xmesh, ymesh, zmesh), dtype=dtype(dim.nxgrid,dim.nzgrid)))
    if 'yz' in planes:
        ymesh, zmesh = np.meshgrid(grid.y, grid.z, indexing='ij')
        xmesh = np.full(ymesh.shape, par.xyz0[0])
        surfaces.append(np.rec.array((slices.yz, xmesh, ymesh, zmesh), dtype=dtype(dim.nygrid,dim.nzgrid)))
    # Set up the 3D view.
    print("Initializing...")
    fig = plt.figure()
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
    cols = [ax.plot_surface(s.x, s.y, s.z, facecolors=cmap(norm(s.value[0])), linewidth=0, shade=False)
            for s in surfaces]
    mappable = cm.ScalarMappable(cmap=cmap, norm=norm)
    mappable.set_array([])
    cb = fig.colorbar(mappable, shrink=0.8, aspect=20)
    cb.set_label(field)
    timestamp = ax.set_title("$t = {:#.4G}$".format(t[0]))
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
        timestamp.set_text("$t = {:#.4G}$".format(t[num]))
        return cols
    anim = FuncAnimation(fig, update, nt, fargs=(surfaces,), interval=1, blit=False, repeat=False)
    #anim.save('slices3d.mp4', writer='mencoder')
    plt.show()
