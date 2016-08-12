def animate_slices_compareruns(field='uu1', datadir1='data/', datadir2='data/', proc=-1, extension='xz',
                   format='native', tmin=0., tmax=1.e38, wait=0.,
                   amin=0., amax=1., transform='', oldfile=False,
                   makemovie=False):
    """
    read 2D slice files from two different runs
    note that the runs must have same precision and sizes

    Options:

     field --- which variable to slice
     datadir1 --- path to data directory of first simulation
     datadir2 --- path to data directory of second imulation
     proc --- an integer giving the processor to read a slice from
     extension --- which plane of xy,xz,yz,Xz. for 2D this should be overwritten.
     format --- endian. one of little, big, or native (default)
     tmin --- start time
     tmax --- end time
     amin --- minimum value for image scaling
     amax --- maximum value for image scaling
     transform --- insert arbitrary numerical code to modify the slice
     wait --- pause in seconds between animation slices
     makemovie --- assemble an animation if makemovie option is set
    """

    import pylab as plt

    datadir1 = os.path.expanduser(datadir1)
    if proc < 0:
        filename1 = join(datadir1, 'slice_' + field + '.' + extension)
    else:
        filename1 = join(datadir1, 'proc' + str(proc),
                        'slice_' + field + '.' + extension)

    datadir2 = os.path.expanduser(datadir2)
    if proc < 0:
        filename2 = join(datadir2, 'slice_' + field + '.' + extension)
    else:
        filename2 = join(datadir2, 'proc' + str(proc),
                        'slice_' + field + '.' + extension)

    # Read the global dimensions.
    dim = read_dim(datadir1, proc)
    if dim.precision == 'D':
        precision = 'd'
    else:
        precision = 'f'

    # Set up slice plane.
    if extension == 'xy' or extension == 'Xy':
        hsize = dim.nx
        vsize = dim.ny
    if extension == 'xz':
        hsize = dim.nx
        vsize = dim.nz
    if extension == 'yz':
        hsize = dim.ny
        vsize = dim.nz
    plane1 = np.zeros((vsize, hsize), dtype=precision)
    plane2 = np.zeros((vsize, hsize), dtype=precision)

    infile1 = npfile(filename1, endian=format)
    infile2 = npfile(filename2, endian=format)

    #ax = plt.axes()
    #ax.set_xlabel('')
    #ax.set_ylabel('')
    #ax.set_ylim
    #ax.get_xaxis().set_visible(False)
    #ax.get_yaxis().set_visible(False)

    fig, (ax1,ax2) = plt.subplots(1,2)
    #fig.suptitle('Re = 400', fontsize=20)
    image1 = ax1.imshow(plane1, vmin=amin, vmax=amax)
    image2 = ax2.imshow(plane2, vmin=amin, vmax=amax)
    ax1.set_xlabel('')
    ax1.set_ylabel('')
    ax1.get_xaxis().set_visible(False)
    ax1.get_yaxis().set_visible(False)
    ax2.set_xlabel('')
    ax2.set_ylabel('')
    ax2.get_xaxis().set_visible(False)
    ax2.get_yaxis().set_visible(False)

    # Get the figure manager for real-time image display.
    manager = plt.get_current_fig_manager()
    manager.show()

    ifirst = True
    islice = 0
    files = []

    while True:
        try:
            raw_data1 = infile1.fort_read(precision)
            raw_data2 = infile2.fort_read(precision)
        except ValueError:
            break
        except TypeError:
            break

        if oldfile:
            t = raw_data1[-1]
            plane1 = raw_data1[:-1].reshape(vsize, hsize)
            plane2 = raw_data2[:-1].reshape(vsize, hsize)
        else:
            t = raw_data1[-2]
            plane1 = raw_data1[:-2].reshape(vsize, hsize)
            plane2 = raw_data2[:-2].reshape(vsize, hsize)

        if transform:
            exec('plane = plane' + transform)

        if t > tmin and t < tmax:
            title = 't = %11.3e' % t
            #fig.set_title(title)
            image1.set_data(plane1)
            image2.set_data(plane2)
            manager.canvas.draw()

            if ifirst:
                #print "----islice----------t---------min-------max-------delta" # Python 2
                print("----islice----------t---------min-------max-------delta")
            #print "%10i %10.3e %10.3e %10.3e %10.3e" \ # Python 2
                #% (islice, t, plane.min(), plane.max(), plane.max() - plane.min()) # Python 2
            print("{0:10} {1:10.3e} {2:10.3e} {3:10.3e} {4:10.3e}".format(islice, t, plane1.min(), plane1.max(), plane1.max() - plane1.min()))

            if(makemovie):
                fname = '_tmp%03d.png' % islice
                fig.savefig(fname)
                files.append(fname)

            ifirst = False
            islice += 1

            sleep(wait)

    infile1.close()
    infile2.close()

    if(makemovie):
        print('Making movie animation.mpg - this make take a while')
        # SC: Not all systems use mencoder. Need to change this into ffmpeg.
        os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=24 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o animation.mpg")
        os.system("rm _tmp*.png")
