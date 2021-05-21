def make_movie(field='uu1',datadir='data/',proc=-1,extension='xz',format='native',
                tmin=0.,tmax=1.e38,amin=0.,amax=1.,transform='',oldfile=False):
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
    """
    import os
    from pencil.io import npfile

    datadir = os.path.expanduser(datadir)
    if proc < 0:
        filename = datadir+'/slice_'+field+'.'+extension
    else:
        filename = datadir+'/proc'+str(proc)+'/slice_'+field+'.'+extension

    # global dim
    param = read_param(datadir)

    dim = read_dim(datadir,proc)
    if dim.precision == 'D':
        precision = 'd'
    else:
        precision = 'f'

    # set up slice plane
    if (extension == 'xy' or extension == 'Xy'):
        hsize = dim.nx
        vsize = dim.ny
    if (extension == 'xz'):
        hsize = dim.nx
        vsize = dim.nz
    if (extension == 'yz'):
        hsize = dim.ny
        vsize = dim.nz
    plane = N.zeros((vsize,hsize),dtype=precision)

    infile = npfile(filename,endian=format)

    files = []
    fig = P.figure(figsize=(5,10))
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
            plane = raw_data[:-1].reshape(vsize,hsize)
        else:
            slice_z2pos = raw_data[-1]
            t = raw_data[-2]
            plane = raw_data[:-2].reshape(vsize,hsize)

        if transform:
            exec('plane = plane'+transform)

        if (t > tmin and t < tmax):
            ax.cla()
            ax.imshow(plane)
            fname = '_tmp%03d.png'%islice
            print('Saving frame', fname)
            fig.savefig(fname)
            files.append(fname)

            if ifirst:
                print("----islice----------t---------min-------max-------delta")
            print("%10i %10.3e %10.3e %10.3e %10.3e" % (islice,t,plane.min(),plane.max(),plane.max()-plane.min()))

            ifirst = False
            islice += 1

    print('Making movie animation.mpg - this make take a while')
    os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=24 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o animation.mpg")
    os.system("rm _tmp*.png")
    infile.close()

