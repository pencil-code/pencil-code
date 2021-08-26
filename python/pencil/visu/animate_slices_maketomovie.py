def make_movie(
    field="uu1",
    datadir="data/",
    proc=-1,
    extension="xz",
    format="native",
    tmin=0.0,
    tmax=1.0e38,
    amin=0.0,
    amax=1.0,
    transform="",
    oldfile=False,
    norm=None,
    save=None,
    figsize=(16, 4),
):
    """
    read 2D slice files and assemble an animation in a mpg movie.

    Quickly written from the example at http://matplotlib.sourceforge.net/faq/howto_faq.html

    Options:

     field --- which variable to slice
     datadir --- path to data directory
     proc --- an integer giving the processor to read a slice from
     extension --- which plane of xy,xz,yz,Xz. for 2D this should be overwritten.
     format --- endian. one of little, big, or native (default)
     tmin --- start time
     tmax --- end time
     amin --- minimum value for image scaling
     amax --- maximum value for image scaling
     transform --- insert arbitrary numerical code to modify the slice
     norm --- scales calar data
     save -- directory to save file
     figsize --- tuple containing the size of the figure
    """
    import os
    from pencil.io import npfile
    from pencil import read
    import numpy as np
    import pylab as plt
    from matplotlib import colors

    # Global configuration:
    # lines
    plt.rcParams["lines.linewidth"] = 2
    plt.rcParams["lines.color"] = "k"
    # font
    plt.rcParams["font.size"] = 30

    plt.rcParams["font.family"] = "serif"
    # legend
    plt.rcParams["legend.fontsize"] = 20
    plt.rcParams["legend.fancybox"] = False
    plt.rcParams["legend.numpoints"] = 2
    plt.rcParams["legend.shadow"] = False
    plt.rcParams["legend.frameon"] = False
    # latex
    plt.rc("text", usetex=True)
    plt.rcParams["text.latex.preamble"] = [r"\usepackage{amsmath}"]

    datadir = os.path.expanduser(datadir)
    if proc < 0:
        filename = os.path.join(datadir, "slice_" + field + "." + extension)
    else:
        filename = os.path.join(
            datadir, "proc" + str(proc) + "/slice_" + field + "." + extension
        )

    # global dim
    # param = read.param(datadir)

    dim = read.dim(datadir, proc)

    if dim.precision == "D":
        precision = "d"
    else:
        precision = "f"

    grid = read.grid(datadir=datadir, trim=True)
    # set up slice plane
    if extension == "xy" or extension == "Xy":
        hsize = dim.nx
        vsize = dim.ny
        xlabel = "x"
        ylabel = "y"
        x = grid.x
        y = grid.y
    if extension == "xz":
        hsize = dim.nx
        vsize = dim.nz
        xlabel = "x"
        ylabel = "z"
        x = grid.x
        y = grid.z
    if extension == "yz":
        hsize = dim.ny
        vsize = dim.nz
        xlabel = "y"
        ylabel = "z"
        x = grid.y
        y = grid.z

    plane = np.zeros((vsize, hsize), dtype=precision)

    infile = npfile(filename, endian=format)

    files = []
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(
        left=0.12, bottom=0.1, right=0.98, top=0.96, wspace=0.23, hspace=0.2
    )
    ax = fig.add_subplot(111)

    ifirst = True
    islice = 0
    while 1:
        try:
            raw_data = infile.fort_read(precision)
        except ValueError:
            break
        except TypeError:
            break

        if oldfile:
            t = raw_data[-1]
            plane = raw_data[:-1].reshape(vsize, hsize)
        else:
            slice_z2pos = raw_data[-1]
            t = raw_data[-2]
            plane = raw_data[:-2].reshape(vsize, hsize)

        if transform:
            exec("plane = plane" + transform)

        if t > tmin and t < tmax:
            ax.cla()
            title = "t = %11.3e" % t
            ax.set_title(title)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

            ax.imshow(
                plane,
                origin="lower",
                vmin=amin,
                vmax=amax,
                norm=norm,
                cmap="hot",
                extent=[x[0], x[-1], y[0], y[-1]],
                aspect=1,
            )
            fname = "_tmp%03d.png" % islice
            print("Saving frame", fname)
            fig.savefig(fname)
            files.append(fname)

            if ifirst:
                print("----islice----------t---------min-------max-------delta")
            print(
                "%10i %10.3e %10.3e %10.3e %10.3e"
                % (islice, t, plane.min(), plane.max(), plane.max() - plane.min())
            )

            ifirst = False
            islice += 1
        if t > tmax:
            break

    print("Making movie animation.mpg - this make take a while")
    os.system(
        "mencoder 'mf://_tmp*.png' -mf type=png:fps=24 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o animation.mpg"
    )
    if save:
        os.system(f"mv _tmp*.png {save}")
        print(f"Moving files to {save}")
    else:
        os.system("rm _tmp*.png")
        print("Removing all files")
    infile.close()

