# $Id: tracers.py,v 1.4 2012-03-21 09:55:49 iomsn Exp $
#
# Reads the tracer files, composes a color map.
#
# Author: Simon Candelaresi (iomsn1@googlemail.com).
#
#

import struct
import array
import numpy as np
import os
import pencil_old as pc
import multiprocessing as mp


def tracers(traceField = 'bb', hMin = 2e-3, hMax = 2e4, lMax = 500, tol = 1e-2,
                interpolation = 'weighted', trace_sub = 1, intQ = [''], varfile = 'VAR0',
                ti = -1, tf = -1,
                integration = 'simple', datadir = 'data/', destination = 'tracers.dat', nproc = 1):
    """
    Trace streamlines from the VAR files and integrate quantity 'intQ' along them.

    call signature::
    
      tracers(field = 'bb', hMin = 2e-3, hMax = 2e2, lMax = 500, tol = 2e-3,
                interpolation = 'weighted', trace_sub = 1, intQ = '', varfile = 'VAR0',
                ti = -1, tf = -1,
                datadir = 'data', destination = 'tracers.dat', nproc = 1)
    
    Trace streamlines of the vectofield 'field' from z = z0 to z = z1 and integrate
    quantities 'intQ' along the lines. Creates a 2d mapping as in 'streamlines.f90'.
    
    Keyword arguments:
    
     *traceField*:
       Vector field used for the streamline tracing.
        
     *hMin*:
       Minimum step length for and underflow to occur.
       
     *hMax*:
       Parameter for the initial step length.
       
     *lMax*:
       Maximum length of the streamline. Integration will stop if l >= lMax.
       
     *tol*:
       Tolerance for each integration step. Reduces the step length if error >= tol.
     
     *interpolation*:
       Interpolation of the vector field.
       'mean': takes the mean of the adjacent grid point.
       'weighted': weights the adjacent grid points according to their distance.
       
     *trace_sub*:
       Number of sub-grid cells for the seeds.
       
     *intQ*:
       Quantities to be integrated along the streamlines.
     
     *varfile*:
       Varfile to be read.
       
      *integration*:
        Integration method.
        'simple': low order method.
        'RK6': Runge-Kutta 6th order.
        
      *ti*:
        Initial VAR file index for tracer time sequences. Overrides 'varfile'.
        
      *tf*:
        Final VAR file index for tracer time sequences. Overrides 'varfile'.
        
      *datadir*:
        Directory where the data is stored.
        
     *destination*:
       Destination file.
       
     *nproc*:
       Number of cores for multi core computation.
    """

    # returns the tracers for the specified starting locations
    def subTracers(q, vv, p, tracers0, iproc, hMin = 2e-3, hMax = 2e4, lMax = 500, tol = 1e-2, 
                   interpolation = 'weighted', integration = 'simple', intQ = ['']):
        
        tracers = tracers0
        mapping = np.zeros((tracers.shape[0], tracers.shape[1], 3))
        
        for ix in range(tracers.shape[0]):
            for iy in range(tracers.shape[1]):
                xx = tracers[ix, iy, 2:5].copy()
                s = pc.stream(vv, p, interpolation = interpolation, integration = integration, hMin = hMin, hMax = hMax, lMax = lMax, tol = tol, xx = xx)
                tracers[ix, iy, 2:5] = s.tracers[s.sl-1]
                tracers[ix, iy, 5] = s.l
                if (any(intQ == 'curlyA')):
                    for l in range(s.sl-1):
                        aaInt = pc.vecInt((s.tracers[l+1] + s.tracers[l])/2, aa, p, interpolation)
                        tracers[ix, iy, 6] += np.dot(aaInt, (s.tracers[l+1] - s.tracers[l]))
                
                # create the color mapping
                if (tracers[ix, iy, 4] > grid.z[-2]):
                    if (tracers[ix, iy, 0] - tracers[ix, iy, 2]) > 0:
                        if (tracers[ix, iy, 1] - tracers[ix, iy, 3]) > 0:
                            mapping[ix, iy, :] = [0,1,0]
                        else:
                            mapping[ix, iy, :] = [1,1,0]
                    else:
                        if (tracers[ix, iy, 1] - tracers[ix, iy, 3]) > 0:
                            mapping[ix, iy, :] = [0,0,1]
                        else:
                            mapping[ix, iy, :] = [1,0,0]
                else:
                    mapping[ix, iy, :] = [1,1,1]
        
        q.put((tracers, mapping, iproc))
        
    
    # multi core setup
    if (np.isscalar(nproc) == False) or (nproc%1 != 0):
        print("error: invalid processor number")
        return -1
    queue = mp.Queue()
    
    # read the data
    # make sure to read the var files with the correct magic
    if (traceField == 'bb'):
        magic = 'bb'
    if (traceField == 'jj'):
        magic = 'jj'
    if (traceField == 'vort'):
        magic = 'vort'
    
    # convert intQ string into list
    if (isinstance(intQ, list) == False):
        intQ = [intQ]
    intQ = np.array(intQ)
    
    grid = pc.read_grid(datadir = datadir, trim = True, quiet = True) 
    dim  = pc.read_dim(datadir = datadir)    
    tol2 = tol**2
    
    # check if user wants a tracer time series
    if ((ti%1 == 0) and (tf%1 == 0) and (ti >= 0) and (tf >= ti)):
        series = True
        n_times = tf-ti+1
    else:
        series = False
        n_times = 1
    
    tracers = np.zeros([int(trace_sub*dim.nx), int(trace_sub*dim.ny), n_times, 6+len(intQ)])
    mapping = np.zeros([int(trace_sub*dim.nx), int(trace_sub*dim.ny), n_times, 3])
    t = np.zeros(n_times)
    
    for tIdx in range(n_times):
        if series:
            varfile = 'VAR' + str(tIdx)
        
        # read the data
        var = pc.read_var(varfile = varfile, datadir = datadir, magic = magic, quiet = True, trimall = True)   
        grid = pc.read_grid(datadir = datadir, quiet = True, trim = True)
        t[tIdx] = var.t
        
        # extract the requested vector traceField
        vv = getattr(var, traceField)
        if (any(intQ == 'curlyA')):
            aa = var.aa
        
        # initialize the parameters
        p = pc.pClass()
        p.dx = var.dx; p.dy = var.dy; p.dz = var.dz
        p.Ox = var.x[0]; p.Oy = var.y[0]; p.Oz = var.z[0]
        p.Lx = grid.Lx; p.Ly = grid.Ly; p.Lz = grid.Lz
        p.nx = dim.nx; p.ny = dim.ny; p.nz = dim.nz
        
        # initialize the tracers
        for ix in range(int(trace_sub*dim.nx)):
            for iy in range(int(trace_sub*dim.ny)):
                tracers[ix, iy, tIdx, 0] = grid.x[0] + int(grid.dx/trace_sub)*ix
                tracers[ix, iy, tIdx, 2] = tracers[ix, iy, tIdx, 0]
                tracers[ix, iy, tIdx, 1] = grid.y[0] + int(grid.dy/trace_sub)*iy
                tracers[ix, iy, tIdx, 3] = tracers[ix, iy, tIdx, 1]
                tracers[ix, iy, tIdx, 4] = grid.z[0]
            
        # declare vectors
        xMid    = np.zeros(3)
        xSingle = np.zeros(3)
        xHalf   = np.zeros(3)
        xDouble = np.zeros(3)
        
        tmp = []
        subTracersLambda = lambda queue, vv, p, tracers, iproc: \
            subTracers(queue, vv, p, tracers, iproc, hMin = hMin, hMax = hMax, lMax = lMax, tol = tol,
                       interpolation = interpolation, integration = integration, intQ = intQ)
        proc = []
        for iproc in range(nproc):
            proc.append(mp.Process(target = subTracersLambda, args = (queue, vv, p, tracers[iproc::nproc,:,tIdx,:], iproc)))
        for iproc in range(nproc):
            proc[iproc].start()
        for iproc in range(nproc):
            tmp.append(queue.get())
        for iproc in range(nproc):
            proc[iproc].join()
        for iproc in range(nproc):
            tracers[tmp[iproc][2]::nproc,:,tIdx,:], mapping[tmp[iproc][2]::nproc,:,tIdx,:] = (tmp[iproc][0], tmp[iproc][1])
        for iproc in range(nproc):
            proc[iproc].terminate()
        
    tracers = np.copy(tracers.swapaxes(0, 3), order = 'C')
    if (destination != ''):
        f = open(datadir + destination, 'wb')
        f.write(np.array(trace_sub, dtype = 'float32'))
        # write tracers into file
        for tIdx in range(n_times):
            f.write(t[tIdx].astype('float32'))
            f.write(tracers[:,:,tIdx,:].astype('float32'))
        f.close()
        
    tracers = tracers.swapaxes(0, 3)
    tracers = tracers.swapaxes(0, 1)
    mapping = mapping.swapaxes(0, 1)

    return tracers, mapping, t


def read_tracers(datadir = 'data/', fileName = 'tracers.dat', zlim = [], head_size = 3, post = False):
    """
    Reads the tracer files and composes a color map.

    call signature::

      tracers, mapping, t = read_tracers(fileName = 'tracers.dat', datadir = 'data/', zlim = [], head_size = 3, post = False)

    Reads from the tracer files and computes the color map according to
    A R Yeates and G Hornig 2011 J. Phys. A: Math. Theor. 44 265501
    doi:10.1088/1751-8113/44/26/265501.
    Returns the tracer values, the color mapping and the times of the snapshots.
    The color mapping can be plotted with:
    pc.animate_interactive(mapping[:,::-1,:,:], t, dimOrder = (2,1,0,3))

    Keyword arguments:

      *datadir*:
        Data directory.

      *fileName*:
        Name of the tracer file.

      *zlim*:
        The upper limit for the field line mapping at which a field line is considered
        to have reached the upper boundary.

      *head_size*:
        Size of the Fortran header in binary data. Most of the time this is 3.
        For the St Andrews cluster it is 5.
      
      *post*:
        If True reads the post processed tracer file 'data/tracers.dat'.
    """
    
    class data_struct:
        def __init__(self):
            self.xi = []
            self.yi = []
            self.xf = []
            self.yf = []
            self.zf = []
            self.l = []
            self.q = []

    data = []
    data = data_struct()

    # compute the offset in order to skip Fortran's header byte
    if (post):
        head_size = 0
        off = 2
    if (head_size == 3):
        off = 2
    if (head_size == 5):
        off = 3

    # read the cpu structure
    dim = pc.read_dim(datadir = datadir)
    if (dim.nprocz > 1):
        print(": number of cores in z-direction > 1")
        return -1

    # read the parameters
    params = pc.read_param(datadir = datadir, quiet = True)

    # read the grid
    grid = pc.read_grid(datadir = datadir, quiet = True)

    # determine the file structure
    if (post):
        n_proc = 1
        tracer_file = open(datadir+fileName, 'rb')
        trace_sub = struct.unpack("f", tracer_file.read(4))[0]
        tracer_file.close()
        n_times = int((os.path.getsize(datadir+fileName)-4)/(4*7*int(dim.nx*trace_sub)*int(dim.ny*trace_sub)))
    # sub sampling of the tracers
    if (not(post)):
        n_proc = dim.nprocx*dim.nprocy
        trace_sub = params.trace_sub
        n_times = int(os.path.getsize(datadir+'proc0/'+fileName)/(4*(head_size + 7*np.floor(dim.nx*trace_sub)*np.floor(dim.ny*trace_sub)/dim.nprocx/dim.nprocy)))

    # prepare the output arrays
    tracers = np.zeros((int(dim.nx*trace_sub), int(dim.ny*trace_sub), n_times, 7))
    mapping = np.zeros((int(dim.nx*trace_sub), int(dim.ny*trace_sub), n_times, 3))

    # temporary arrays for one core
    if (post):
        tracers_core = tracers
        mapping_core = mapping
    else:
        tracers_core = np.zeros((int(int(dim.nx*trace_sub)/dim.nprocx), int(int(dim.ny*trace_sub)/dim.nprocy), n_times, 7))
        mapping_core = np.zeros((int(int(dim.nx*trace_sub)/dim.nprocx), int(np.floor(dim.ny*trace_sub)/dim.nprocy), n_times, 3))

    # set the upper z-limit to the domain boundary
    if zlim == []:
        zlim = grid.z[-dim.nghostz-1]

    # read the data from all cores
    for i in range(n_proc):
        # read the cpu structure
        if (post):
            dim_core = pc.read_dim(datadir = datadir, proc = -1)
            dim_core.ipx = 0
            dim_core.ipy = 0
        else:
            dim_core = pc.read_dim(datadir = datadir, proc = i)
        stride = int(dim_core.nx*trace_sub)*int(dim_core.ny*trace_sub)
        llen = head_size + 7*stride + post

        if (post):
            tracer_file = open(datadir+fileName, 'rb')
        else:
            tracer_file = open(datadir+'proc{0}/'.format(i)+fileName, 'rb')
        tmp = array.array('f')
        tmp.read(tracer_file, int((head_size + post + 7*int(dim_core.nx*trace_sub)*int(dim_core.ny*trace_sub))*n_times)+post)
        tracer_file.close()

        t = []

        for j in range(n_times):
            t.append(tmp[off-1+j*llen])
            data.xi = tmp[off+j*llen          : off+1*stride+j*llen]
            data.yi = tmp[off+1*stride+j*llen : off+2*stride+j*llen]
            data.xf = tmp[off+2*stride+j*llen : off+3*stride+j*llen]
            data.yf = tmp[off+3*stride+j*llen : off+4*stride+j*llen]
            data.zf = tmp[off+4*stride+j*llen : off+5*stride+j*llen]
            data.l  = tmp[off+5*stride+j*llen : off+6*stride+j*llen]
            data.q  = tmp[off+6*stride+j*llen : off+7*stride+j*llen]

            # Squeeze the data into 2d array. This make the visualization much faster.
            for l in range(len(data.xi)):
                tracers_core[l%(int(dim_core.nx*trace_sub)),int(l/(int(dim_core.nx*trace_sub))),j,:] = \
                [data.xi[l], data.yi[l], data.xf[l], data.yf[l], data.zf[l], data.l[l], data.q[l]]
                if data.zf[l] >= zlim:
                    if (data.xi[l] - data.xf[l]) > 0:
                        if (data.yi[l] - data.yf[l]) > 0:
                            mapping_core[l%(int(dim_core.nx*trace_sub)),int(l/(int(dim_core.nx*trace_sub))),j,:] = [0,1,0]
                        else:
                            mapping_core[l%(int(dim_core.nx*trace_sub)),int(l/(int(dim_core.nx*trace_sub))),j,:] = [1,1,0]
                    else:
                        if (data.yi[l] - data.yf[l]) > 0:
                            mapping_core[l%(int(dim_core.nx*trace_sub)),int(l/(int(dim_core.nx*trace_sub))),j,:] = [0,0,1]
                        else:
                            mapping_core[l%(int(dim_core.nx*trace_sub)),int(l/(int(dim_core.nx*trace_sub))),j,:] = [1,0,0]
                else:
                    mapping_core[l%(int(dim_core.nx*trace_sub)),int(l/(int(dim_core.nx*trace_sub))),j,:] = [1,1,1]

            # copy single core data into total data arrays
            if (not(post)):
                tracers[np.round(dim_core.ipx*int(dim_core.nx*trace_sub)):np.round((dim_core.ipx+1)*np.floor(dim_core.nx*trace_sub)), \
                        np.round(dim_core.ipy*int(dim_core.ny*trace_sub)):np.round((dim_core.ipy+1)*np.floor(dim_core.ny*trace_sub)),j,:] = \
                        tracers_core[:,:,j,:]
                mapping[np.round(dim_core.ipx*int(dim_core.nx*trace_sub)):np.round((dim_core.ipx+1)*np.floor(dim_core.nx*trace_sub)), \
                        np.round(dim_core.ipy*int(dim_core.ny*trace_sub)):np.round((dim_core.ipy+1)*np.floor(dim_core.ny*trace_sub)),j,:] = \
                        mapping_core[:,:,j,:]

    # swap axes for post evaluation
    tracers = tracers.swapaxes(0, 1)
    mapping = mapping.swapaxes(0, 1)

    return tracers, mapping, t


def tracer_movie(datadir = 'data/', tracerFile = 'tracers.dat',
                 fixedFile = 'fixed_points.dat', zlim = [],
                 head_size = 3, hm = 1,
                 imageDir = './', movieFile = 'fixed_points.mpg',
                 fps = 5.0, bitrate = 1800):
    """
    Plots the color mapping together with the fixed points.
    Creates a movie file.

    call signature::

      tracer_movie(datadir = 'data/', tracerFile = 'tracers.dat',
                 fixedFile = 'fixed_points.dat', zlim = [],
                 head_size = 3, hm = 1,
                 imageDir = './', movieFile = 'fixed_points.mpg',
                 fps = 5.0, bitrate = 1800)

    Plots the field line mapping together with the fixed points and creates
    a movie file.

      *datadir*:
        Data directory.

      *tracerFile*:
        Name of the tracer file.

      *fixedFile*:
        Name of the fixed points file.

      *zlim*:
        The upper limit for the field line mapping at which a field line is considered
        to have reached the upper boundary.

      *head_size*:
        Size of the fortran header in binary data. Most of the time this is 3.
        For the St Andrews cluster it is 5.

      *hm*:
        Header multiplication factor in case Fortran's binary data writes extra large
        header. For most cases hm = 1 is sufficient. For the cluster in St Andrews use hm = 2.

      *imageDir*:
        Directory with the temporary png files.

      *movieFile*:
        Output file for the movie. Ending should be 'mpg', since the compression
        format is mpg.

      *fps*:
        Frames per second of the animation.

      *bitrate*:
        Bitrate of the movie file. Set to higher value for higher quality.
    """
    
    import pylab as plt

    # read the mapping and the fixed point positions
    tracers, mapping, t = pc.read_tracers(datadir = datadir, fileName = tracerFile, zlim = zlim, head_size = head_size)
    fixed = pc.read_fixed_points(datadir = datadir, fileName = fixedFile, hm = hm)

    # read the parameters for the domain boundaries
    params = pc.read_param(quiet = True)
    domain = [params.xyz0[0], params.xyz1[0], params.xyz0[1], params.xyz1[1]]

    # determine the how much faster the fixed pints have been written out than the color mapping
    advance = np.ceil(float(len(fixed.t))/len(mapping[0,0,:,0]))

    # determine the colors for the fixed points
    colors = np.zeros(np.shape(fixed.q) + (3,))
    colors[:,:,:] = 0.
    print(np.shape(colors))
    for j in range(len(colors[:,0,0])):
        for k in range(len(colors[0,:,0])):
            if fixed.q[j,k] >= 0:
                colors[j,k,1] = colors[j,k,2] = (1-fixed.q[j,k]/np.max(np.abs(fixed.q[:,k])))
                colors[j,k,0] = fixed.q[j,k]/np.max(np.abs(fixed.q[:,k]))
            else:
                colors[j,k,0] = colors[j,k,1] = (1+fixed.q[j,k]/np.max(np.abs(fixed.q[:,k])))
                colors[j,k,2] = -fixed.q[j,k]/np.max(np.abs(fixed.q[:,k]))

    # prepare the plot
    width = 6
    height = 6
    plt.rc("figure.subplot", left=(60/72.27)/width)
    plt.rc("figure.subplot", right=(width-20/72.27)/width)
    plt.rc("figure.subplot", bottom=(50/72.27)/height)
    plt.rc("figure.subplot", top=(height-20/72.27)/height)
    figure = plt.figure(figsize=(width, height))

    for k in range(len(fixed.x[0,:])):
        dots = plt.plot(fixed.x[0,k], fixed.y[0,k], 'o', c = colors[0,k,:])
    image = plt.imshow(zip(*mapping[:,::-1,0,:]), interpolation = 'nearest', extent = domain)
    j = 0
    frameName = imageDir + 'images%06d.png'%j
    imageFiles = []
    imageFiles.append(frameName)
    figure.savefig(frameName)

    for j in range(1,len(fixed.t)):
        #time.sleep(0.5)
        figure.clear()
        for k in range(len(fixed.x[j,:])):
            dots = plt.plot(fixed.x[j,k], fixed.y[j,k], 'o', c = colors[j,k,:])
        image = plt.imshow(zip(*mapping[:,::-1,np.floor(j/advance),:]), interpolation = 'nearest', extent = domain)
        frameName = imageDir + 'images%06d.png'%j
        imageFiles.append(frameName)
        figure.savefig(frameName)

    # convert the images into a mpg file
    mencodeCommand = "mencoder 'mf://"+imageDir+"images*.png' -mf type=png:fps="+np.str(fps)+" -ovc lavc -lavcopts vcodec=mpeg4:vhq:vbitrate="+np.str(bitrate)+" -ffourcc MP4S -oac copy -o "+movieFile
    os.system(mencodeCommand)

    # remove the image files
    for fname in imageFiles:
        os.remove(fname)
