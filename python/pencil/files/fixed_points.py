# $Id: fixed_points.py
#
# Creates and reads the fixed point values.
#
# Author: Simon Candelaresi (iomsn1@googlemail.com).
#
#

import struct
import numpy as np
import pencil as pc
import math as m


def read_fixed_points(datadir = 'data/', fileName = 'fixed_points.dat', hm = 1):
    """
    Reads the fixed points files.

    call signature::

      fixed = read_fixed_points(datadir = 'data/', fileName = 'fixed_points.dat', hm = 1)

    Reads from the fixed points files. Returns the fixed points positions.

    Keyword arguments:

      *datadir*:
        Data directory.

      *fileName*:
        Name of the fixed points file.

      *hm*:
        Header multiplication factor in case Fortran's binary data writes extra large
        header. For most cases hm = 1 is sufficient. For the cluster in St Andrews use hm = 2.
    """


    class data_struct:
        def __init__(self):
            self.t = []
            self.fidx = [] # number of fixed points at this time
            self.x = []
            self.y = []
            self.q = []

    # read the cpu structure
    dim = pc.read_dim(datadir = datadir)
    if (dim.nprocz > 1):
        print "error: number of cores in z-direction > 1"

    data = []

    # read the data
    fixed_file = open(datadir+fileName, 'rb')
    tmp = fixed_file.read(4*hm)

    data = data_struct()
    eof = 0
    # length of longest array of fixed points
    fixedMax = 0
    if tmp == '':
        eof = 1
    while (eof == 0):
        data.t.append(struct.unpack("<"+str(hm+1)+"f", fixed_file.read(4*(hm+1)))[0])
        n_fixed = int(struct.unpack("<"+str(2*hm+1)+"f", fixed_file.read(4*(2*hm+1)))[1+hm/2])

        x = list(np.zeros(n_fixed))
        y = list(np.zeros(n_fixed))
        q = list(np.zeros(n_fixed))
        for j in range(n_fixed):
            x[j] = struct.unpack("<"+str(hm+1)+"f", fixed_file.read(4*(hm+1)))[-1]
            y[j] = struct.unpack("<f", fixed_file.read(4))[0]
            q[j] = struct.unpack("<"+str(hm+1)+"f", fixed_file.read(4*(hm+1)))[0]
        data.x.append(x)
        data.y.append(y)
        data.q.append(q)
        data.fidx.append(n_fixed)

        tmp = fixed_file.read(4*hm)
        if tmp == '':
            eof = 1

        if fixedMax < len(x):
            fixedMax = len(x)

    fixed_file.close()

    # add NaN to fill up the times with smaller number of fixed points
    fixed = data_struct()
    for i in range(len(data.t)):
        annex = list(np.zeros(fixedMax - len(data.x[i])) + np.nan)
        fixed.t.append(data.t[i])
        fixed.x.append(data.x[i] + annex)
        fixed.y.append(data.y[i] + annex)
        fixed.q.append(data.q[i] + annex)
        fixed.fidx.append(data.fidx[i])

    fixed.t = np.array(fixed.t)
    fixed.x = np.array(fixed.x)
    fixed.y = np.array(fixed.y)
    fixed.q = np.array(fixed.q)
    fixed.fidx = np.array(fixed.fidx)

    return fixed


#
# under construction
#
#def fixed_points(datadir = 'data/', fileName = 'fixed_points_post.dat', varfile = 'VAR0', ti = -1, tf = -1,
                 #traceField = 'bb', hMin = 2e-3, hMax = 2e4, lMax = 500, tol = 1e-2,
                 #interpolation = 'mean', trace_sub = 1, integration = 'simple'):
    #"""
    #Reads the fixed points files.

    #call signature::

      #fixed = fixed_points(datadir = 'data/', fileName = 'fixed_points_post.dat')

    #Reads from the fixed points files. Returns the fixed points positions.

    #Keyword arguments:

      #*datadir*:
        #Data directory.

      #*fileName*:
        #Name of the fixed points file.

     #*varfile*:
       #Varfile to be read.
       
      #*ti*:
        #Initial VAR file index for tracer time sequences. Overrides 'varfile'.
        
      #*tf*:
        #Final VAR file index for tracer time sequences. Overrides 'varfile'.        

     #*traceField*:
       #Vector field used for the streamline tracing.
        
     #*hMin*:
       #Minimum step length for and underflow to occur.
       
     #*hMax*:
       #Parameter for the initial step length.
       
     #*lMax*:
       #Maximum length of the streamline. Integration will stop if l >= lMax.
       
     #*tol*:
       #Tolerance for each integration step. Reduces the step length if error >= tol.
     
     #*interpolation*:
       #Interpolation of the vector field.
       #'mean': takes the mean of the adjacent grid point.
       #'weighted': weights the adjacent grid points according to their distance.
       
     #*trace_sub*:
       #Number of sub-grid cells for the seeds for the initial mapping.
       
     #*intQ*:
       #Quantities to be integrated along the streamlines.
     
      #*integration*:
        #Integration method.
        #'simple': low order method.
        #'RK6': Runge-Kutta 6th order.
    #"""


    #class data_struct:
        #def __init__(self):
            #self.t = []
            #self.fidx = [] # number of fixed points at this time
            #self.x = []
            #self.y = []
            #self.q = []
    
    
    ## Computes rotation along one edge.
    #def edege(vv, p, sx, sy, diff1, diff2, phiMin, rec):
        #dtot = m.atan2(diff1[0]*diff2[1] - diff2[0]diff1[1], diff1[0]*diff2[0] - diff1[1]diff2[1])
        #if ((abs(dtot) > phiMin) and (rec < 4)):
            #xm = 0.5*(sx[0]+sx[1])
            #ym = 0.5*(sy[0]+sy[1])
            ## trace intermediate field line
            #s = pc.stream(vv, p, lMax = 500, xx = np.array([xm,ym,p.Oz]))
            #tracer = s.tracer[s.sl-1,:]
            ## discard any streamline which does not converge or hits the boundary
            #if ((tracer[5] >= 500) or (tracer(4) < p.Oz+p.Lz-p.dz)):
                #dtot = 0.
            #else:
                #diffm = np.array([tracer[2] - tracer[0], tracer[3] - tracer[1]])
                #if (sum(diffm**2) != 0):
                    #diffm = diffm / np.sqrt(sum(diffm**2))
                #dtot = edge(vv, p, [sx[0], xm], [sy[0], ym], diff1, diffm, phiMin, rec+1) +
                       #edge(vv, p, [xm, sx[1]], [ym, sy[1]], diffm, diff2, phiMin, rec+1)
        #else:
             #return dtot
        
        
    ## Finds the Poincare index of this grid cell.
    #def pIndex(vv, p, sx, sy, diff, phiMin, poincare):
        #poincare = 0
        #poincare += edge(vv, p, [sx[0], sx[1]], [sy[0], sy[0]], diff[0,:], diff[1,:], phiMin, 0)
        #poincare += edge(vv, p, [sx[1], sx[1]], [sy[0], sy[1]], diff[1,:], diff[2,:], phiMin, 0)
        #poincare += edge(vv, p, [sx[1], sx[0]], [sy[1], sy[1]], diff[2,:], diff[3,:], phiMin, 0)
        #poincare += edge(vv, p, [sx[0], sx[0]], [sy[1], sy[0]], diff[3,:], diff[0,:], phiMin, 0)        
        #return poincare     
           
           
    #phiMin = np.pi/8.
    
    ## make sure to read the var files with the correct magic
    #if (traceField == 'bb'):
        #magic = 'bb'
    #if (traceField == 'jj'):
        #magic = 'jj'
    #if (traceField == 'vort'):
        #magic = 'vort'
        
    ## read the cpu structure
    #dim = pc.read_dim(datadir = datadir)
    #if (dim.nprocz > 1):
        #print "error: number of cores in z-direction > 1"

    #var = pc.read_var(varfile = varfile, datadir = datadir, magic = magic, quiet = True, trimall = True)
    #grid = pc.read_grid(datadir = datadir, quiet = True, trim = True)
    #vv = getattr(var, traceField)
    
    ## initialize the parameters
    #p = pc.pClass()
    #p.dx = var.dx; p.dy = var.dy; p.dz = var.dz
    #p.Ox = var.x[0]; p.Oy = var.y[0]; p.Oz = var.z[0]
    #p.Lx = grid.Lx; p.Ly = grid.Ly; p.Lz = grid.Lz
    #p.nx = dim.nx; p.ny = dim.ny; p.nz = dim.nz
        
    ## create the initial mapping
    #tracers = pc.tracers(traceField = 'bb', hMin = hMin, hMax = hMax, lMax = lMax, tol = tol,
                         #interpolation = interpolation, trace_sub = trace_sub, varfile = varfile,
                         #integration = integration, datadir = datadir, destination = '')
    
    #diff = np.zeros((4,2))
    ## find fixed points
    #for ix in range(var.nx*trace_sub-1):
        #for iy in range(var.ny*trace_sub-1):
            ## compute Poincare index around this point (!= 0 for potential fixed point)
            #diff[1,:] = tracers[iy, ix, 0, 2:4] - tracers[iy, ix, 0, 0:2]
            #diff[2,:] = tracers[iy, ix+1, 0, 2:4] - tracers[iy, ix+1, 0, 0:2]
            #diff[3,:] = tracers[iy+1, ix+1, 0, 2:4] - tracers[iy+1, ix+1, 0, 0:2]
            #diff[4,:] = tracers[iy+1, ix, 0, 2:4] - tracers[iy+1, ix, 0, 0:2]
            #if (sum(np.sum(diff**2, axis = 1) != 0) == True):
                #diff = np.swapaxes(np.swapaxes(diff, 0, 1) / np.sqrt(np.sum(diff**2, axis = 1)), 0, 1)
            #poincare = pIndex(vv, p, tracers[iy, ix:ix+2, 0, 0], tracers[iy:iy+2, ix, 0, 1], diff, phiMin)
            
            #if (abs(poincare) > 5): # use 5 instead of 2pi to account for rounding errors
                ## subsampe to get starting point for iteration
                #nt = 4
                #xmin = tracers[iy, ix, 0, 0]
                #ymin = tracers[iy, ix, 0, 1]
                #xmax = tracers[iy, ix+1, 0, 0]
                #ymax = tracers[iy+1, ix, 0, 1]
                #xx = np.zeros((nt**2,3))
                #tracersSub = np.zeros((nt**2,5))
                #i1 = 0
                #for j1 in range(nt):
                    #for k1 in range(nt):
                        #xx[i1,0] = xmin + j1/(nt-1.)*(xmax - xmin)
                        #xx[i1,1] = ymin + k1/(nt-1.)*(ymax - ymin)
                        #xx[i1,2] = p.Oz
                        #i1 += 1
                #for it1 in range(nt**2):
                    #s = pc.stream(vv, p, lMax = 500, xx = xx[it1,:])
                    #tracersSub[it1,0:2] = xx[it1,0:2]
                    #tracersSub[it1,2:] = s.tracers[s.sl-1,:]
                #min2 = 1e6
                #minx = xmin
                #miny = ymin
                #i1 = 0
                #for j1 in range(nt):
                    #for k1 in range(nt):
                        #diff2 = (tracersSub[i1, 2] - tracers2[i1, 0])**2 + (tracersSub[i1, 3] - tracers2[i1, 1])**2
                        #if (diff2 < min2):
                            #min2 = diff2
                            #minx = xmin + j1/(nt-1.)*(xmax - xmin)
                            #miny = ymin + k1/(nt-1.)*(ymax - ymin)
                        #it1 += 1
                
                ## get fixed point from this starting position using Newton's method
                #dl = np.min(var.dx, dar.dy)/100.    # step-size for calculatin the Jacobian by finite differences
                #it = 0
                ## tracers used to find the fixed point
                #tracersNull = np.zeros((4,3))
                #point = np.array([minx, miny])
                #while True:
                    ## trace field lines at original point and for Jacobian:
                    ## (second order seems to be enough)
                    #xx = np.zeros((5,3))
                    #xx[0,:] = np.array([xmin, ymin, p.Oz])
                    #xx[1,:] = np.array([xmin-dl, ymin, p.Oz])
                    #xx[2,:] = np.array([xmin+dl, ymin, p.Oz])
                    #xx[3,:] = np.array([xmin, ymin-dl, p.Oz])
                    #xx[4,:] = np.array([xmin, ymin+dl, p.Oz])
                    #for it1 in range(5):
                        #s = pc.stream(vv, p, lMax = 500, xx = xx[it1,:])
                        #tracersNull[:2] = xx[it1]
                        #tracersNull[2:] = s.tracers[s.sl-1,0:2]
                    #break
                
            
    #fixed.t = np.array(fixed.t)
    #fixed.x = np.array(fixed.x)
    #fixed.y = np.array(fixed.y)
    #fixed.q = np.array(fixed.q)
    #fixed.fidx = np.array(fixed.fidx)

    #return fixed
