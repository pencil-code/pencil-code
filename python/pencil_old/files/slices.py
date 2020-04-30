# $Id$
#
# Read slice files.
#
# Author: J. Oishi (joishi@amnh.org).
#
#
import os
import numpy as np
from ..files.npfile import npfile
from ..files.dim import read_dim
from time import sleep
from os.path import join

# slice file format is either
#   plane,t (old style)
#   plane,t,slice_z2pos (new style)


def read_slices(field='uu1', datadir='data', proc=-1,
                extension='xz', format='native', oldfile=False):
    """
    Read 2D slice files and return an array of (nslices, vsize, hsize).
    """
    datadir = os.path.expanduser(datadir)
    if proc < 0:
        filename = join(datadir, 'slice_' + field + '.' + extension)
    else:
        filename = join(datadir, 'proc' + str(proc),
                        'slice_' + field + '.' + extension)

    # Read the global dimensions.
    dim = read_dim(datadir, proc)
    if dim.precision == 'D':
        precision = 'd'
    else:
        precision = 'f'

    # Set up slice plane.
    if extension.startswith('xy'):
        hsize = dim.nx
        vsize = dim.ny
    if extension.startswith('xz'):
        hsize = dim.nx
        vsize = dim.nz
    if extension.startswith('yz'):
        hsize = dim.ny
        vsize = dim.nz

    infile = npfile(filename, endian=format)

    islice = 0
    t = np.zeros(1, dtype=precision)
    slices = np.zeros(1, dtype=precision)

    while True:
        try:
            raw_data = infile.fort_read(precision)
        except ValueError:
            break
        except TypeError:
            break

        if oldfile:
            t = np.concatenate((t, raw_data[-1:]))
            slices = np.concatenate((slices, raw_data[:-1]))
        else:
            t = np.concatenate((t, raw_data[-2:-1]))
            slices = np.concatenate((slices, raw_data[:-2]))
        islice += 1

    output = slices[1:].reshape(islice, vsize, hsize)

    return output, t[1:]


def animate_slices(field='uu1', datadir='data/', proc=-1, extension='xz',
                   format='native', tmin=0., tmax=1.e38, wait=0.,
                   amin=0., amax=1., transform='', oldfile=False):
    """
    read 2D slice files and assemble an animation.

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
     wait --- pause in seconds between animation slices
    """

    import pylab as plt

    datadir = os.path.expanduser(datadir)
    if proc < 0:
        filename = join(datadir, 'slice_' + field + '.' + extension)
    else:
        filename = join(datadir, 'proc' + str(proc),
                        'slice_' + field + '.' + extension)

    # Read the global dimensions.
    dim = read_dim(datadir, proc)
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
    plane = np.zeros((vsize, hsize), dtype=precision)

    infile = npfile(filename, endian=format)

    ax = plt.axes()
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_ylim

    image = plt.imshow(plane, vmin=amin, vmax=amax)

    # Get the figure manager for real-time image display.
    manager = plt.get_current_fig_manager()
    manager.show()

    ifirst = True
    islice = 0
    while True:
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
            t = raw_data[-2]
            plane = raw_data[:-2].reshape(vsize, hsize)

        if transform:
            exec('plane = plane' + transform)

        if t > tmin and t < tmax:
            title = 't = %11.3e' % t
            ax.set_title(title)
            image.set_data(plane)
            manager.canvas.draw()

            if ifirst:
                #print "----islice----------t---------min-------max-------delta" # Python 2
                print("----islice----------t---------min-------max-------delta")
            #print "%10i %10.3e %10.3e %10.3e %10.3e" \ # Python 2
                #% (islice, t, plane.min(), plane.max(), plane.max() - plane.min()) # Python 2
            print("{0:10} {1:10.3e} {2:10.3e} {3:10.3e} {4:10.3e}".format(islice, t, plane.min(), plane.max(), plane.max() - plane.min()))

            ifirst = False
            islice += 1

            sleep(wait)

    infile.close()


def animate_multislices(field=['uu1'], datadir='data/', proc=-1,
                        extension='xz', format='native', tmin=0., tmax=1.e38,
                        amin=0., amax=1., transform='plane[0]',
                        oldfile=False, outfile=""):
    """
    Read a list of 2D slice files, combine them, and assemble an animation.

    Options:

     field --- list of variables to slice
     datadir --- path to data directory
     proc --- an integer giving the processor to read a slice from
     extension --- which plane of xy,xz,yz,Xz. for 2D this should be overwritten.
     format --- endian. one of little, big, or native (default)
     tmin --- start time
     tmax --- end time
     amin --- minimum value for image scaling
     amax --- maximum value for image scaling
     transform --- insert arbitrary numerical code to combine the slices
     outfile --- if set, write the slice values in the text file
    """

    import pylab as plt

    datadir = os.path.expanduser(datadir)
    if outfile != "":
        outslice = open(outfile, "w")
    filename = []
    if proc < 0:
        for i in field:
            filename += [datadir + '/slice_' + i + '.' + extension]
    else:
        for i in field:
            filename += [datadir + '/proc' +
                         str(proc) + '/slice_' + i + '.' + extension]

    # Read the global dimensions.
    dim = read_dim(datadir, proc)
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
    plane = []
    infile = []
    for i in filename:
        plane += [np.zeros((vsize, hsize), dtype=precision)]

        infile += [npfile(i, endian=format)]

    ax = plt.axes()
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_ylim

    exec('plotplane =' + transform)
    image = plt.imshow(plotplane, vmin=amin, vmax=amax)

    # Get the figure manager for real-time image display.
    manager = plt.get_current_fig_manager()
    manager.show()

    ifirst = True
    islice = 0
    while True:
        try:
            raw_data = []
            for i in infile:
                raw_data += [i.fort_read(precision)]
        except ValueError:
            break
        except TypeError:
            break

        if oldfile:
            t = raw_data[0][-1]
            for i in range(len(raw_data)):
                plane[i] = raw_data[i][:-1].reshape(vsize, hsize)
        else:
            t = raw_data[0][-2]
            for i in range(len(raw_data)):
                plane[i] = raw_data[i][:-2].reshape(vsize, hsize)

        exec('plotplane =' + transform)

        if t > tmin and t < tmax:
            title = 't = %11.3e' % t
            ax.set_title(title)
            image.set_data(plotplane)
            manager.canvas.draw()

            if ifirst:
                #print "----islice----------t---------min-------max-------delta" # Python 2
                print("----islice----------t---------min-------max-------delta")
            #print "%10i %10.3e %10.3e %10.3e %10.3e" % \ # Python 2
                #(islice, t, plotplane.min(), plotplane.max(), # Python 2
                 #plotplane.max() - plotplane.min()) # Python 2
            print("{0:10} {1:10.3e} {2:10.3e} {3:10.3e} {4:10.3e}".format(islice, t, plotplane.min(), plotplane.max(),
                 plotplane.max() - plotplane.min()))
            if outfile != "":
                #outslice.write("%10i %10.3e %10.3e %10.3e %10.3e" % # Python 2
                               #(islice, t, plotplane.min(), plotplane.max(), # Python 2
                                #plotplane.max() - plotplane.min())) # Python 2
                outslice.write("{0:10} {1:10.3e} {2:10.3e} {3:10.3e} {4:10.3e}".format(islice, t, plotplane.min(), plotplane.max(), plotplane.max() - plotplane.min()))
                outslice.write("\n")

            ifirst = False
            islice += 1

    for i in infile:
        i.close()
    if outfile != "":
        outslice.close()


def time_slices(field=['uu1'], datadir='data/', proc=-1, extension='xz',
                format='native', tmin=0., tmax=1.e38, amin=0., amax=1.,
                transform='plane[0]', dtstep=1, deltat=0,
                oldfile=False, outfile=""):
    """
    Read a list of 1D slice files, combine them, and plot the slice in
    one dimension, and time in the other one.

    Options:

     field --- list of variables to slice
     datadir --- path to data directory
     proc --- an integer giving the processor to read a slice from
     extension --- which plane of xy,xz,yz,Xz. for 2D this should be overwritten.
     format --- endian. one of little, big, or native (default)
     tmin --- start time
     tmax --- end time
     amin --- minimum value for image scaling
     amax --- maximum value for image scaling
     transform --- insert arbitrary numerical code to combine the slices
     dtstep --- only plot every dt step
     deltat --- if set to nonzero, plot at fixed time interval rather than step
     outfile --- if set, write the slice values in the text file
    """

    import pylab as plt

    datadir = os.path.expanduser(datadir)
    if outfile != "":
        outslice = open(outfile, "w")
    filename = []
    if proc < 0:
        for i in field:
            filename += [datadir + '/slice_' + i + '.' + extension]
    else:
        for i in field:
            filename += [datadir + '/proc' +
                         str(proc) + '/slice_' + i + '.' + extension]

    # Read the global dimensions.
    dim = read_dim(datadir, proc)
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
    plane = []
    infile = []
    for i in filename:
        plane += [np.zeros((vsize, hsize), dtype=precision)]

        infile += [npfile(i, endian=format)]

    ifirst = True
    islice = 0
    plotplane = []
    dt = 0
    nextt = tmin
    while True:
        try:
            raw_data = []
            for i in infile:
                raw_data += [i.fort_read(precision)]
        except ValueError:
            break
        except TypeError:
            break

        if oldfile:
            t = raw_data[0][-1]
            for i in range(len(raw_data)):
                plane[i] = raw_data[i][:-1].reshape(vsize, hsize)
        else:
            t = raw_data[0][-2]
            for i in range(len(raw_data)):
                plane[i] = raw_data[i][:-2].reshape(vsize, hsize)

        exec('tempplane =' + transform)

        if t > tmin and t < tmax:
            if dt == 0:
                plotplane += tempplane.tolist()

                if ifirst:
                    #print "----islice----------t---------min-------max-------delta" # Python 2
                    print("----islice----------t---------min-------max-------delta")
                #print "%10i %10.3e %10.3e %10.3e %10.3e" % \ # Python 2
                    #(islice, t, tempplane.min(), tempplane.max(), # Python 2
                     #tempplane.max() - tempplane.min()) # Python 2
                print("{0:10} {1:10.3e} {2:10.3e} {3:10.3e} {4:10.3e}".format(islice, t, tempplane.min(), tempplane.max(), tempplane.max() - tempplane.min()))
                if outfile != "":
                    outslice.write(
                        #"%10i %10.3e %10.3e %10.3e %10.3e" % # Python 2
                        #(islice, # Python 2
                         #t, # Python 2
                         #tempplane.min(), # Python 2
                            #tempplane.max(), # Python 2
                            #tempplane.max() - # Python 2
                            #tempplane.min())) # Python 2
                        "{0:10} {1:10.3e} {2:10.3e} {3:10.3e} {4:10.3e}".format(
                            islice,
                            t,
                            tempplane.min(),
                            tempplane.max(),
                            tempplane.max() -
                            tempplane.min()))                        
                    outslice.write("\n")

                ifirst = False
                islice += 1
                nextt = t + deltat
            if deltat == 0:
                dt = (dt + 1) % dtstep
            elif t >= nextt:
                dt = 0
                nextt = t + deltat
            else:
                dt = 1

    ax = plt.axes()
    ax.set_xlabel('t')
    ax.set_ylabel('y')
    ax.set_ylim
    plt.imshow(np.array(plotplane).reshape(islice, vsize).transpose(),
               vmin=amin, vmax=amax)
    manager = plt.get_current_fig_manager()
    manager.show()

    for i in infile:
        i.close()
    if outfile != "":
        outslice.close()


def make_movie(field='uu1', datadir='data/', proc=-1, extension='xz',
               format='native', tmin=0., tmax=1.e38, amin=0., amax=1.,
               transform='', oldfile=False):
    """
    read 2D slice files and assemble an animation in a mpg movie.

    Quickly written from the example at
    http://matplotlib.sourceforge.net/faq/howto_faq.html

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

    import pylab as plt

    datadir = os.path.expanduser(datadir)
    if proc < 0:
        filename = datadir + '/slice_' + field + '.' + extension
    else:
        filename = datadir + '/proc' + \
            str(proc) + '/slice_' + field + '.' + extension

    # Read the global dimensions.
    dim = read_dim(datadir, proc)
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
    plane = np.zeros((vsize, hsize), dtype=precision)

    infile = npfile(filename, endian=format)

    files = []
    fig = plt.figure(figsize=(5, 10))
    ax = fig.add_subplot(111)

    ifirst = True
    islice = 0
    while True:
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
            t = raw_data[-2]
            plane = raw_data[:-2].reshape(vsize, hsize)

        if transform:
            exec('plane = plane' + transform)

        if t > tmin and t < tmax:
            ax.cla()
            ax.imshow(plane, vmin=amin, vmax=amax)
            fname = '_tmp%03d.png' % islice
            print('Saving frame' + fname)
            fig.savefig(fname)
            files.append(fname)

            if ifirst:
                #print "----islice----------t---------min-------max-------delta" # Python 2
                print("----islice----------t---------min-------max-------delta")
            #print "%10i %10.3e %10.3e %10.3e %10.3e" % \ # Python 2
                #(islice, t, plane.min(), plane.max(), plane.max() - plane.min()) # Python 2
            print("{0:10} {1:10.3e} {2:10.3e} {3:10.3e} {4:10.3e}".format(islice, t, plane.min(), plane.max(), plane.max() - plane.min()))

            ifirst = False
            islice += 1

    #print 'Making movie animation.mpg - this make take a while'
    print('Making movie animation.mpg - this make take a while')
    # SC: Not all systems use mencoder. Need to change this into ffmpeg.
    os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=24 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o animation.mpg")
    os.system("rm _tmp*.png")
    infile.close()

def make_movie_cart(field='uu1', datadir='data/', proc=-1, extension='xz',
               format='native', tmin=0., tmax=1.e38, amin=0., amax=1.,
               transform='', oldfile=False):
    """
    read 2D slice files and assemble an animation in a mpg movie.

    Quickly written from the example at
    http://matplotlib.sourceforge.net/faq/howto_faq.html

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

    import pylab as plt
    import sys
    from pencil_old.files.var import read_var

    mkmvvar = read_var(trimall=True)    
    r2d,phi2d = np.meshgrid(mkmvvar.x,mkmvvar.y)
    x2d=r2d*np.cos(phi2d)
    y2d=r2d*np.sin(phi2d)

    datadir = os.path.expanduser(datadir)
    if proc < 0:
        filename = datadir + '/slice_' + field + '.' + extension
    else:
        filename = datadir + '/proc' + \
            str(proc) + '/slice_' + field + '.' + extension

    # Read the global dimensions.
    dim = read_dim(datadir, proc)
    if dim.precision == 'D':
        precision = 'd'
    else:
        precision = 'f'

    # Set up slice plane.
    if extension == 'xy' or extension == 'Xy':
        hsize = dim.nx
        vsize = dim.ny
    if extension == 'xz':
        print('only works for xy')
        sys.stop
    if extension == 'yz':
        print('only works for xy')
        sys.stop
    plane = np.zeros((vsize, hsize), dtype=precision)

    infile = npfile(filename, endian=format)

    files = []
    fig = plt.figure(figsize=(5, 10))
    ax = fig.add_subplot(111)

    ifirst = True
    islice = 0
    while True:
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
            t = raw_data[-2]
            plane = raw_data[:-2].reshape(vsize, hsize)

        if transform:
            exec('plane = plane' + transform)

        if t > tmin and t < tmax:
            ax.set_aspect('equal')
            ax.cla()
            ax.contourf(x2d, y2d, plane, 256)
            fname = '_tmp%03d.png' % islice
            print('Saving frame' + fname)
            fig.savefig(fname)
            files.append(fname)

            if ifirst:
                #print "----islice----------t---------min-------max-------delta" # Python 2
                print("----islice----------t---------min-------max-------delta")
            #print "%10i %10.3e %10.3e %10.3e %10.3e" % \ # Python 2
                #(islice, t, plane.min(), plane.max(), plane.max() - plane.min()) # Python 2
            print("{0:10} {1:10.3e} {2:10.3e} {3:10.3e} {4:10.3e}".format(islice, t, plane.min(), plane.max(), plane.max() - plane.min()))

            ifirst = False
            islice += 1

    #print 'Making movie animation.mpg - this make take a while'
    print('Making movie animation.mpg - this make take a while')
    # SC: Not all systems use mencoder. Need to change this into ffmpeg.
    os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=24 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o animation.mpg")
    os.system("rm _tmp*.png")
    infile.close()



def animate_slices_multi(field='uu1', datadir1='data/', datadir2='data/', proc=-1, extension='xz',
                   format='native', tmin=0., tmax=1.e38, wait=0.,
                   amin=0., amax=1., transform='', oldfile=False,
                   makemovie=False):
    """
    read 2D slice files and assemble an animation.
    version that does this for two different runs, neat for comparrison
    runs must have same precision and sizes

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

def make_movie_crossflow(field='uu1', datadir='data/', proc=-1, extension='yz',
               format='native', tmin=0., tmax=1.e38, amin=0., amax=1.,
               transform='', oldfile=False):
    """
    read 2D slice files and assemble an animation in a mpg movie.

    Quickly written from the example at
    http://matplotlib.sourceforge.net/faq/howto_faq.html

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

    import pylab as plt
    import matplotlib.patches as patches

    datadir = os.path.expanduser(datadir)
    if proc < 0:
        filename = datadir + '/slice_' + field + '.' + extension
    else:
        filename = datadir + '/proc' + \
            str(proc) + '/slice_' + field + '.' + extension

    # Read the global dimensions.
    dim = read_dim(datadir, proc)
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
    plane = np.zeros((vsize, hsize), dtype=precision)

    infile = npfile(filename, endian=format)

    files = []
    fig = plt.figure(figsize=(5, 10))
    ax = fig.add_subplot(111)
    ax.add_patch(patches.Rectangle(
            (220,0),
            40,
            320,
            color='gray'
        )
    )
#
#    ax.add_patch(patches.Rectangle(
#            (220,0),
#            80,
#            240,
#            hatch='/'
#        )
#    )

    ifirst = True
    islice = 0
    while True:
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
            t = raw_data[-2]
            plane = raw_data[:-2].reshape(vsize, hsize)

        if transform:
            exec('plane = plane' + transform)

        if t > tmin and t < tmax:
            ax.cla()
            ax.imshow(plane, vmin=amin, vmax=amax)
            ax.add_patch(patches.Rectangle(
                    (220,0),
                    40,
                    320,
                    color='gray'
                )
            )
            fname = '_tmp%03d.png' % islice
            print('Saving frame' + fname)
            fig.savefig(fname)
            files.append(fname)


def animate_slices_crossflow(field='uu1', datadir='data/', proc=-1, extension='yz',
                   format='native', tmin=0., tmax=1.e38, wait=0.,
                   amin=0., amax=1., transform='', oldfile=False):
    """
    read 2D slice files and assemble an animation.

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
     wait --- pause in seconds between animation slices
    """

    import pylab as plt
    import matplotlib.patches as ptc
    import matplotlib.ticker as ticker

    datadir = os.path.expanduser(datadir)
    if proc < 0:
        filename = join(datadir, 'slice_' + field + '.' + extension)
    else:
        filename = join(datadir, 'proc' + str(proc),
                        'slice_' + field + '.' + extension)

    # Read the global dimensions.
    dim = read_dim(datadir, proc)
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
    plane = np.zeros((vsize, hsize), dtype=precision)

    infile = npfile(filename, endian=format)

    ax = plt.axes()
    ax.set_xlabel('y')
    ax.set_ylabel('x')
    ax.set_ylim


    name_list_x = ('0', '2D', '4D', '6D')
    #name_list_z = ('12D', '10D', '8D', '6D', '4D','2D','0')
    name_list_z = ( '0', '2D', '4D', '6D','8D','10D','12D',)
    pos_list_x = np.array([0,80,160,240])
    pos_list_z = np.array([0,80,160,240,320,400,480])

    name_list_x = ('0', '5D', '10D')
    name_list_z = ( '0', '5D', '10D', '15D','20D')
    pos_list_x = np.array([0,159,319])
    pos_list_z = np.array([0,159,319,479,639])
    #pos_list = np.arange(len(name_list))

#ax = plt.axes()
    #ax.xaxis.set_major_locator(ticker.FixedLocator((pos_list_x)))
    #ax.xaxis.set_major_formatter(ticker.FixedFormatter((name_list_x)))
    #ax.yaxis.set_major_locator(ticker.FixedLocator((pos_list_z)))
    #ax.yaxis.set_major_formatter(ticker.FixedFormatter((name_list_z)))
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)


    # turbulent particle sim xy-vew, Ddx=40
    #r=20
    #xy = [120,240]
#    art = ptc.Circle(xy,r,color='black')
#    ax.add_artist(art)
    # uniform flow sim xy-vew, Ddx=64
   # r=16
   # xy = [160,320]
    # Pencil-tests, Re=400, Ddx=96
    #r=48
    #xy = [480,960]
    # Pencil-tests, Re=100, Ddx=64
    r=32
    xy = [320,640]
    art = ptc.Circle(xy,r,color='black')
    ax.add_artist(art)
     

    #ax.add_patch(patches.Rectangle(
    #        (220,0),
    #        40,
    #        320,
    #        color='gray'
    #    )
    #)

    image = plt.imshow(plane, vmin=amin, vmax=amax)

    # Get the figure manager for real-time image display.
    manager = plt.get_current_fig_manager()
    manager.show()

    ifirst = True
    islice = 0
    while True:
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
            t = raw_data[-2]
            plane = raw_data[:-2].reshape(vsize, hsize)

        if transform:
            exec('plane = plane' + transform)

        if t > tmin and t < tmax:
            #title = 't = %11.3e' % t
            #ax.set_title(title)
            image.set_data(plane)
            manager.canvas.draw()

            if ifirst:
                #print "----islice----------t---------min-------max-------delta" # Python 2
                print("----islice----------t---------min-------max-------delta")
            #print "%10i %10.3e %10.3e %10.3e %10.3e" \ # Python 2
                #% (islice, t, plane.min(), plane.max(), plane.max() - plane.min()) # Python 2
            print("{0:10} {1:10.3e} {2:10.3e} {3:10.3e} {4:10.3e}".format(islice, t, plane.min(), plane.max(), plane.max() - plane.min()))

            ifirst = False
            fname = '_tmp%03d.eps' % islice
            plt.savefig(fname,fonttype=42,
                  bbox_inches = 'tight', 
                 # pad_inches = 0,
                  dpi = 600, transparant=True)
            #plt.savefig(fname,fonttype=42)
            islice += 1

            sleep(wait)

    infile.close()
