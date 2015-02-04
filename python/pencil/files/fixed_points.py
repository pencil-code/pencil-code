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
    #def dege():
        
    ## Finds the Poincare index of this grid cell.
    #def pIndex():
        
        
        
        
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
    #vv = getattr(var, traceField)
    
    ## create the initial mapping
    #tracers = pc.tracers(traceField = 'bb', hMin = hMin, hMax = hMax, lMax = lMax, tol = tol,
                         #interpolation = interpolation, trace_sub = trace_sub, varfile = varfile,
                         #integration = integration, datadir = datadir, destination = '')
    
    #diff = np.zeros()
    ## find fixed points
    #for ix in range(var.nx*trace_sub):
        #for iy in range(var.ny*trace_sub):
            ## compute Poincare index around this point (!= 0 for potential fixed point)
            #diff[0,:] = tracer[ix+
         #diff(1,:) = (/(tracers2(j+(l-1)*(nx*trace_sub+addx),3)-tracers2(j+(l-1)*(nx*trace_sub+addx),1)) , &
            #(tracers2(j+(l-1)*(nx*trace_sub+addx),4)-tracers2(j+(l-1)*(nx*trace_sub+addx),2))/)
   
    #fixed.t = np.array(fixed.t)
    #fixed.x = np.array(fixed.x)
    #fixed.y = np.array(fixed.y)
    #fixed.q = np.array(fixed.q)
    #fixed.fidx = np.array(fixed.fidx)

    #return fixed
