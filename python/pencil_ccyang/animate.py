#!/usr/bin/python3
# Last Modification: $Id$
#=======================================================================
# animate.py
#
# Facilities for animating the Pencil Code data.
#
# Chao-Chin Yang, 2015-04-22
#=======================================================================
def avg2d(name, direction, datadir='./data', save=False, **kwarg):
    """Animates the time sequence of a line average.

    Positional Argument:
        name
            Name of the average.
        direction
            Direction of the average: 'x', 'y', or 'z'.

    Keyword Arguments:
        datadir
            Path to the data directory.
        save
            If True, save the animation.
        **kwarg
            Keyword arguments passed to _get_range().
    """
    # Chao-Chin Yang, 2015-08-09
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
    _frame_rectangle(t, x, y, avg[name], xlabel=xlabel, ylabel=ylabel, clabel=name, save=save, **kwarg)
#=======================================================================
def slices(field, datadir='./data', tmin=None, **kwarg):
    """Dispatches to the respective animator of video slices.

    Positional Argument:
        field
            Field name.

    Keyword Arguments:
        datadir
            Path to the data directory.
        tmin
            If not None, the slices with time earlier than tmin are
            truncated.
        **kwarg
            Keyword arguments passed to _get_range().
    """
    # Author: Chao-Chin Yang
    # Created: 2015-04-22
    # Last Modified: 2020-11-24
    from . import read
    import numpy as np

    # Read the slices.
    t, s = read.slices(field, datadir=datadir)

    if tmin is not None:
        # Truncate the beginning of the slices.
        indices = np.where(t >= tmin)
        for v in s._fields:
            s = s._replace(**{v: getattr(s, v)[indices]})
        t = t[indices]

    # Get the dimensions.
    dim = read.dimensions(datadir=datadir)
    ndim = (dim.nxgrid > 1) + (dim.nygrid > 1) + (dim.nzgrid > 1)

    # Get the parameters.
    par = read.parameters(datadir=datadir)

    # Get the grid.
    grid = read.grid(datadir=datadir, interface = ndim > 1, par=par)

    # Dispatch.
    if ndim == 1:
        _slices1d(field, t, s, dim, par, grid, **kwarg)
    elif ndim == 2:
        _slices2d(field, t, s, dim, par, grid, **kwarg)
    elif ndim == 3:
        _slices3d(field, t, s, dim, par, grid, **kwarg)
#=======================================================================
###### Local Functions ######
#=======================================================================
def _frame_rectangle(t, x, y, c, xlabel=None, ylabel=None, clabel=None, save=False, **kwarg):
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
        save
            If True, the animation is saved to a video file.
        **kwarg
            Keyword arguments passed to _get_range().
    """
    # Author: Chao-Chin Yang
    # Created: 2015-04-22
    # Last Modified: 2020-06-01
    from collections.abc import Sequence
    from matplotlib.animation import FuncAnimation, writers
    from matplotlib.colors import LogNorm, Normalize
    import matplotlib.pyplot as plt
    from matplotlib.ticker import LogLocator
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
    # Initialize the plot.
    fig = plt.figure()
    ax = fig.gca()
    time_template = r"$t = {:#.4G}$"
    pc = ax.pcolormesh(x, y, c[0].transpose(), norm=norm)
    ax.minorticks_on()
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(y[0], y[-1])
    ax.set_aspect('equal')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(time_template.format(t[0]))
    cb = plt.colorbar(pc)
    if logscale:
        cb.set_ticks(LogLocator(subs=range(10)).tick_values(vmin0,vmax0))
    else:
        cb.ax.minorticks_on()
    cb.set_label(clabel)
    # Loop over each time and update the plot.
    def update(i):
        print("\rAnimating ({:6.1%})......".format((i+1)/len(t)),
              end='', flush=True)
        ax.set_title(time_template.format(t[i]))
        pc.set_array(c[i].ravel(order='F'))
        if vmin_dynamic and vmax_dynamic:
            pc.set_clim(vmin[i], vmax[i])
        elif vmin_dynamic:
            pc.set_clim(vmin=vmin[i])
        elif vmax_dynamic:
            pc.set_clim(vmax=vmax[i])
    # Save or show the animation.
    anim = FuncAnimation(fig, update, len(t), interval=40, repeat=False)
    if save:
        fname = "pc_anim.mp4"
        FFMpegWriter = writers["ffmpeg"]
        metadata = dict(title='Pencil Code Animation', artist='ccyang')
        writer = FFMpegWriter(fps=60, metadata=metadata, bitrate=16000)
        anim.save(fname, writer=writer)
        print("Saved the animation in {}. ".format(fname))
    else:
        plt.show(block=True)
#=======================================================================
def _get_range(t, data, center=False, drange='full', logscale=False, tmin=None):
    """Determines the data range to be plotted.

    Positional Arguments:
        t
            Sequence of time points.
        data
            Sequence in time of numpy arrays.

    Keyword Arguments:
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
        logscale
            Whether or not the color map is in logarithmic scale.
        tmin
            If not None, the range determination is restricted to
            t >= tmin.  No effect if drange is 'dynamic'.
    """
    # Chao-Chin Yang, 2015-07-30
    from collections.abc import Sequence
    import numpy as np
    from scipy.integrate import simps
    from sys import float_info
    # User-defined range.
    if not isinstance(drange, str):
        return drange
    # Find the data range at each time.
    nt = len(t)
    vmin, vmax = np.empty(nt), np.empty(nt)
    tiny = float_info.min
    for i, d in enumerate(data):
        if d.dtype.names is None:
            vmin1, vmax1 = d.min(), d.max()
            indices = d > 0
            vposmin = d[indices].min() if indices.any() else np.inf
        else:
            vmin1, vmax1, vposmin = np.inf, -np.inf, np.inf
            for n in d.dtype.names:
                d1 = d[n]
                vmin1 = min(vmin1, d1.min())
                vmax1 = max(vmax1, d1.max())
                indices = d1 > 0
                vposmin = min(vposmin, (d1[indices].min() if indices.any() else np.inf))
        vmin[i], vmax[i] = vmin1, vmax1
        # Guard against non-positive number if logscale is True.
        if logscale and vmin[i] <= 0:
            vmin[i] = vposmin if vposmin is not np.inf else tiny
    # Check the type of range requested.
    if drange == 'dynamic':
        pass
    else:
        if tmin is not None:
            indices = t >= tmin
            t, vmin, vmax = t[indices], vmin[indices], vmax[indices]
            nt = len(t)
        if drange == 'full':
            vmin, vmax = vmin.min(), vmax.max()
        elif drange == 'mean':
            dt = t[-1] - t[0]
            vmin, vmax = simps(vmin, t) / dt, simps(vmax, t) / dt
        else:
            raise ValueError("Unknown type of range '{}'".format(drange))
    # Function to center the range.
    if not isinstance(center, bool) and not isinstance(center, (int, float)):
        raise ValueError("Invalid value of the keyword center = {}".format(center))
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
    # Center the range if requested.
    if center:
        if type(vmin) is np.ndarray:
            for i in range(nt):
                vmin[i], vmax[i] = get_centered_range(vmin[i], vmax[i])
        else:
            vmin, vmax = get_centered_range(vmin, vmax)
    return vmin, vmax
#=======================================================================
def _slices1d(field, t, slices, dim, par, grid, ylog=False, **kwarg):
    """Animates video slices from a 1D model.

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
            Grid coordinates supplied by read.grid(interface=False).
        ylog
            If True, y-axis is plotted logarithmically.

    Keyword Arguments:
        **kwarg
            Keyword arguments passed to matplotlib.
    """
    # Author: Chao-Chin Yang
    # Created: 2018-06-05
    # Last Modified: 2019-10-03
    from matplotlib.animation import FuncAnimation
    import matplotlib.pyplot as plt

    # Determine the slice line.
    if dim.nxgrid > 1:
        s = slices.xy[:,:,0]
        idir = 0
    elif dim.nygrid > 1:
        s = slices.xy[:,0,:]
        idir = 1
    else:
        s = slices.xz[:,0,:]
        idir = 2

    # Name the direction.
    direction = "xyz"[idir]

    # Replace the labels if in **kwarg.
    xlabel = kwarg.pop("xlabel", direction)
    ylabel = kwarg.pop("ylabel", field)

    # Replace the y limits if in **kwarg.
    ylim = kwarg.pop("ylim", (s.min(), s.max()))

    # Get the x coordinates.
    x = getattr(grid, direction)

    # Generate the first frame.
    fig = plt.figure(**kwarg)
    ax = fig.gca()
    if ylog:
        line, = ax.semilogy(x, s[0], "o:")
    else:
        line, = ax.plot(x, s[0], "o:")
    ax.set_xlim(par.xyz0[idir], par.xyz1[idir])
    if ylim[0] < ylim[1]: ax.set_ylim(ylim)
    ax.minorticks_on()
    ax.grid()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    time_template = r"$t = {:#.4G}$"
    time = ax.set_title(time_template.format(t[0]))

    # Define the update function for animation.
    def update(i):
        time.set_text(time_template.format(t[i]))
        line.set_ydata(s[i])

    # Define the animator.
    anim = FuncAnimation(fig, update, frames=len(s), interval=40, repeat=False)

    plt.show(block=True)
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
    # Author: Chao-Chin Yang
    # Created: 2015-04-22
    # Last Modified: 2015-05-27
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
    _frame_rectangle(t, x, y, getattr(slices, plane),
            xlabel=xlabel, ylabel=ylabel, clabel=clabel, **kwarg)
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
    # Author: Chao-Chin Yang
    # Created: 2015-04-27
    # Last Modified: 2020-06-01
    from collections.abc import Sequence
    from matplotlib.animation import FuncAnimation, writers
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
    if logscale:
        norm = LogNorm(vmin=vmin0, vmax=vmax0)
    else:
        norm = Normalize(vmin=vmin0, vmax=vmax0)
    if logscale:
        a = min(vmin) if vmin_dynamic else vmin
        for p in planes:
            slices[p] = slices[p].clip(a, np.inf)
    # Set the surfaces.
    nt = len(t)
    dtype = lambda nx, ny: [('value', np.float, (nt,nx,ny)),
                            ('x', np.float, (nx+1,ny+1)),
                            ('y', np.float, (nx+1,ny+1)),
                            ('z', np.float, (nx+1,ny+1))]
    surfaces = []
    if 'xy' in planes or 'xy2' in planes:
        xmesh, ymesh = np.meshgrid(grid.x, grid.y, indexing='ij')
        if 'xy' in planes:
            zmesh = np.full(xmesh.shape, par.xyz0[2] - 0.7 * par.lxyz[2])
            surfaces.append(np.rec.array(
                    (slices.xy, xmesh, ymesh, zmesh),
                    dtype=dtype(dim.nxgrid,dim.nygrid)))
        if 'xy2' in planes:
            zmesh = np.full(xmesh.shape, par.xyz1[2])
            surfaces.append(np.rec.array(
                    (slices.xy2, xmesh, ymesh, zmesh),
                    dtype=dtype(dim.nxgrid,dim.nygrid)))
    if 'xz' in planes:
        xmesh, zmesh = np.meshgrid(grid.x, grid.z, indexing='ij')
        ymesh = np.full(xmesh.shape, par.xyz0[1])
        surfaces.append(np.rec.array(
                (slices.xz, xmesh, ymesh, zmesh),
                dtype=dtype(dim.nxgrid,dim.nzgrid)))
    if 'yz' in planes:
        ymesh, zmesh = np.meshgrid(grid.y, grid.z, indexing='ij')
        xmesh = np.full(ymesh.shape, par.xyz0[0])
        surfaces.append(np.rec.array(
                (slices.yz, xmesh, ymesh, zmesh),
                dtype=dtype(dim.nygrid,dim.nzgrid)))
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
    cols = [ax.plot_surface(s.x, s.y, s.z, facecolors=cmap(norm(s.value[0])),
                            linewidth=0, shade=False)
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
        cols = [ax.plot_surface(
                    s.x, s.y, s.z, facecolors=cmap(norm(s.value[num])),
                    linewidth=0, shade=False)
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
    anim = FuncAnimation(fig, update, nt, fargs=(surfaces,), interval=1,
                         blit=False, repeat=False)
    # Save or show the animation.
    if "ffmpeg" in writers.list():
        writer = "ffmpeg"
    elif "mencoder" in writers.list():
        writer = "mencoder"
    else:
        raise RuntimeError("no appropriate animation writer available. ")

    anim = FuncAnimation(fig, update, len(t), interval=40, repeat=False)
    if True:
        fname = "pc_anim.mp4"
        FFMpegWriter = writers["ffmpeg"]
        metadata = dict(title='Pencil Code Animation', artist='ccyang')
        writer = FFMpegWriter(fps=60, metadata=metadata, bitrate=16000)
        anim.save(fname, writer=writer)
        print("Saved the animation in {}. ".format(fname))
    else:
        plt.show()
    plt.show()
