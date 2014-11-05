# $Id: streamlines.py iomsn Exp $
#
# Traces streamlines of a vector field from z0 to z1, similar to 'streamlines.f90'.
#
# Author: Simon Candelaresi (simon.candelaresi@gmail.com, iomsn1@googlemail.com).
# 
#

import numpy as np
import pencil as pc
import vtk as vtk
from vtk.util import numpy_support as VN
import struct


def tracers(traceField = 'bb', hMin = 2e-3, hMax = 2e4, lMax = 500, tol = 1e-2,
                interpolation = 'mean', sub = 1, intQ = [''], varfile = 'VAR0',
                integration = 'simple', datadir = 'data/', destination = 'tracersPost.dat'):
    """
    Trace streamlines and integrate quantity 'intQ' along them.

    call signature::
    
      tracers(field = 'bb', hMin = 2e-3, hMax = 2e2, lMax = 500, tol = 2e-3,
                interpolation = 'mean', sub = 1, intQ = '', varfile = 'VAR0',
                datadir = 'data', destination = 'tracersPost.dat')
    
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
       
     *sub*:
       Number of sub-grid cells for the seeds.
       
     *intQ*:
       Quantities to be integrated along the streamlines.
     
     *varfile*:
       Varfile to be read.
       
      *integration*:
        Integration method.
        'simple': low order method.
        'RK6': Runge-Kutta 6th order.
        
      *datadir*:
        Directory where the data is stored.
        
     *destination*:
       Destination file.
    """

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
    
    tracers = np.zeros([6+len(intQ), sub*dim.nx, sub*dim.ny])
    mapping = np.zeros([sub*dim.nx, sub*dim.ny,3])
    
    # read the data
    var = pc.read_var(varfile = varfile, datadir = datadir, magic = magic, quiet = True, trimall = True)   
    grid = pc.read_grid(datadir = datadir, quiet = True, trim = True)
    
    # extract the requested vector traceField
    vv = getattr(var, traceField)
    if (any(intQ == 'curlyA')):
        aa = var.aa
    
    # initialize the parameters
    p = pClass()
    p.dx = var.dx; p.dy = var.dy; p.dz = var.dz
    p.Ox = var.x[0]; p.Oy = var.y[0]; p.Oz = var.z[0]
    p.Lx = grid.Lx; p.Ly = grid.Ly; p.Lz = grid.Lz
    p.nx = dim.nx; p.ny = dim.ny; p.nz = dim.nz
    
    # initialize the tracers
    for ix in range(int(sub*dim.nx)):
        for iy in range(int(sub*dim.ny)):
            tracers[0, ix, iy] = grid.x[0] + grid.dx/sub*ix
            tracers[2, ix, iy] = tracers[0, ix, iy]
            tracers[1, ix, iy] = grid.y[0] + grid.dy/sub*iy
            tracers[3, ix, iy] = tracers[1, ix, iy]
            tracers[4, ix, iy] = grid.z[0]
        
    # declare vectors
    xMid    = np.zeros(3)
    xSingle = np.zeros(3)
    xHalf   = np.zeros(3)
    xDouble = np.zeros(3)
    
    for ix in range(int(sub*dim.nx)):
        for iy in range(int(sub*dim.ny)):
            outside = False     # True if the streamlines hits the physical boundary
            dh = np.sqrt(hMax*hMin)
            
            # start the streamline tracing
            xx = tracers[2:5, ix, iy].copy()
            if (integration == 'simple'):
                while ((outside == False) and (abs(dh) >= hMin) and (tracers[5, ix, iy] < lMax)):
                    # (a) single step (midpoint method)    
                    xMid = xx + 0.5*dh*vecInt(xx, vv, p, interpolation)
                    xSingle = xx + dh*vecInt(xMid, vv, p, interpolation)
                
                    # (b) two steps with half stepsize
                    xMid = xx + 0.25*dh*vecInt(xx, vv, p, interpolation)
                    xHalf = xx + 0.5*dh*vecInt(xMid, vv, p, interpolation)
                    xMid = xHalf + 0.25*dh*vecInt(xHalf, vv, p, interpolation)
                    xDouble = xHalf + 0.5*dh*vecInt(xMid, vv, p, interpolation)
                
                    # (c) check error (difference between methods)
                    dist2 = np.sum((xSingle-xDouble)**2)
                    if (dist2 > tol2):
                        dh = dh/2
                        if (abs(dh) < hMin):
                            print "Error: stepsize underflow"
                            break
                    else:
                        tracers[5, ix, iy] +=  np.sqrt(np.dot(xx - xDouble, xx - xDouble))
                        # integrate the requested quantity along the field line
                        if (any(intQ == 'curlyA')):
                            aaInt = vecInt((xDouble + xx)/2, aa, p, interpolation)
                            tracers[6, ix, iy] += np.dot(aaInt, xDouble - xx)
                        xx = xDouble.copy()
                        if (dh < hMin):
                            dh = 2*dh
                        if ((dh > hMax) or (np.isnan(dh))):
                            dh = hMax
                        # check if this point lies outside the domain
                        if ((xx[0] < p.Ox-p.dx) or (xx[0] > p.Ox+p.Lx+p.dx) or (xx[1] < p.Oy-p.dy) or (xx[1] > p.Oy+p.Ly+p.dy) or (xx[2] < p.Oz) or (xx[2] > p.Oz+p.Lz)):
                            outside = True

            tracers[2:5, ix, iy] = xx.copy()
            if (tracers[4, ix, iy] > grid.z[-2]):
                if (tracers[0, ix, iy] - tracers[2, ix, iy]) > 0:
                    if (tracers[1, ix, iy] - tracers[3, ix, iy]) > 0:
                        mapping[ix, iy, :] = [0,1,0]
                    else:
                        mapping[ix, iy, :] = [1,1,0]
                else:
                    if (tracers[1, ix, iy] - tracers[3, ix, iy]) > 0:
                        mapping[ix, iy, :] = [0,0,1]
                    else:
                        mapping[ix, iy, :] = [1,0,0]
            else:
                mapping[ix, iy, :] = [1,1,1]
                        

    # write tracers into file
    f = open(datadir + destination, 'wb')
    f.write(np.array(sub, dtype = 'float32'))
    f.write(var.t)
    f.write(tracers.astype('float32'))
    f.close()
    
    return tracers, mapping, var.t


class stream:
    """
    stream -- Holds the traced streamline.
    """
    
    def __init__(self, fileName = 'animation.vtk', interpolation = 'weighted', integration = 'simple', hMin = 2e-3, hMax = 2e4, lMax = 500, tol = 1e-2, iterMax = 1e3, xx = np.array([0,0,0])):
        """
        Creates, and returns the traced streamline.
        
        call signature:
        
          streamInit(fileName = 'animation.vtk', interpolation = 'weighted', integration = 'simple', hMin = 2e-3, hMax = 2e4, lMax = 500, tol = 1e-2, iterMax = 1e3, xx = np.array([0,0,0]))
        
        Trace magnetic streamlines.
        
        Keyword arguments:
        
         *fileName*:
            Name of the file with the field information.
            
         *interpolation*:
            Interpolation of the vector field.
            'mean': takes the mean of the adjacent grid point.
            'weighted': weights the adjacent grid points according to their distance.
       
         *integration*:
            Integration method.
            'simple': low order method.
            'RK6': Runge-Kutta 6th order.
       
         *hMin*:
            Minimum step length for and underflow to occur.
        
         *hMax*:
            Parameter for the initial step length.
        
         *lMax*:
            Maximum length of the streamline. Integration will stop if l >= lMax.
        
         *tol*:
            Tolerance for each integration step. Reduces the step length if error >= tol.
         
         *iterMax*:
            Maximum number of iterations.     
         
         *xx*:
            Initial seed.
        """
        
        # read the data
        reader = vtk.vtkStructuredPointsReader()
        reader.SetFileName(fileName)
        reader.ReadAllVectorsOn()
        reader.Update()
        
        data = reader.GetOutput()
        
        output = reader.GetOutput()
        dim = output.GetDimensions()
        origin = output.GetOrigin()
        spacing = output.GetSpacing()
        data = output.GetPointData()
        
        bb = VN.vtk_to_numpy(data.GetArray('bfield'))
        bb = bb.reshape((dim[0],dim[1],dim[2],3))
        bb = np.swapaxes(bb, 0, 2)
        
        p = pClass()
        p.dx = spacing[0]; p.dy = spacing[1]; p.dz = spacing[2]
        p.Ox = origin[0]; p.Oy = origin[1]; p.Oz = origin[2]
        p.Lx = (dim[0]+1)*spacing[0]; p.Ly = (dim[1]+1)*spacing[1]; p.Lz = (dim[2]+1)*spacing[2]
        p.nx = dim[0]; p.ny = dim[1]; p.nz = dim[2]
        
        self.tracers = np.zeros([iterMax, 3], dtype = 'float32')  # tentative streamline length
        
        tol2 = tol**2
        dh   = np.sqrt(hMax*hMin) # initial step size
        
        # declare vectors
        xMid    = np.zeros(3)
        xSingle = np.zeros(3)
        xHalf   = np.zeros(3)
        xDouble = np.zeros(3)
        
        # initialize the coefficient for the 6th order adaptive time step RK
        a = np.zeros(6); b = np.zeros((6,5)); c = np.zeros(6); cs = np.zeros(6)
        k = np.zeros((6,3))
        a[1] = 0.2; a[2] = 0.3; a[3] = 0.6; a[4] = 1; a[5] = 0.875
        b[1,0] = 0.2;
        b[2,0] = 3/40.; b[2,1] = 9/40.
        b[3,0] = 0.3; b[3,1] = -0.9; b[3,2] = 1.2
        b[4,0] = -11/54.; b[4,1] = 2.5; b[4,2] = -70/27.; b[4,3] = 35/27.
        b[5,0] = 1631/55296.; b[5,1] = 175/512.; b[5,2] = 575/13824.
        b[5,3] = 44275/110592.; b[5,4] = 253/4096.
        c[0] = 37/378.; c[2] = 250/621.; c[3] = 125/594.; c[5] = 512/1771.
        cs[0] = 2825/27648.; cs[2] = 18575/48384.; cs[3] = 13525/55296.
        cs[4] = 277/14336.; cs[5] = 0.25
    
        # do the streamline tracing
        self.tracers[0,:] = xx
        outside = False
        sl = 0
        l = 0
                
        if (integration == 'simple'):
            while ((l < lMax) and (sl < iterMax-1) and (not(np.isnan(xx[0]))) and (outside == False)):
                # (a) single step (midpoint method)                    
                xMid = xx + 0.5*dh*vecInt(xx, bb, p, interpolation)
                #print "xx.shape = ", xx.shape
                #print "xMid.shape = ", xMid.shape
                xSingle = xx + dh*vecInt(xMid, bb, p, interpolation)
            
                # (b) two steps with half stepsize
                xMid = xx + 0.25*dh*vecInt(xx, bb, p, interpolation)
                xHalf = xx + 0.5*dh*vecInt(xMid, bb, p, interpolation)
                xMid = xHalf + 0.25*dh*vecInt(xHalf, bb, p, interpolation)
                xDouble = xHalf + 0.5*dh*vecInt(xMid, bb, p, interpolation)
            
                # (c) check error (difference between methods)
                dist2 = np.sum((xSingle-xDouble)**2)
                if (dist2 > tol2):
                    dh = 0.5*dh
                    if (abs(dh) < hMin):
                        print "Error: stepsize underflow"
                        break
                else:
                    l += np.sqrt(np.sum((xx-xDouble)**2))
                    xx = xDouble.copy()
                    if (abs(dh) < hMin):
                        dh = 2*dh
                    sl += 1
                    self.tracers[sl,:] = xx.copy()
                    if ((dh > hMax) or (np.isnan(dh))):
                        dh = hMax
                    # check if this point lies outside the domain
                    if ((xx[0] < p.Ox-p.dx) or (xx[0] > p.Ox+p.Lx+p.dx) or (xx[1] < p.Oy-p.dy) or (xx[1] > p.Oy+p.Ly+p.dy) or (xx[2] < p.Oz) or (xx[2] > p.Oz+p.Lz)):
                        outside = True
                        
        if (integration == 'RK6'):
            while ((l < lMax) and (sl < iterMax-1) and (not(np.isnan(xx[0]))) and (outside == False)):
                k[0,:] = dh*vecInt(xx, bb, p, interpolation)                            
                k[1,:] = dh*vecInt(xx + b[1,0]*k[0,:], bb, p, interpolation)
                k[2,:] = dh*vecInt(xx + b[2,0]*k[0,:] + b[2,1]*k[1,:], bb, p, interpolation)
                k[3,:] = dh*vecInt(xx + b[3,0]*k[0,:] + b[3,1]*k[1,:] + b[3,2]*k[2,:], bb, p, interpolation)
                k[4,:] = dh*vecInt(xx + b[4,0]*k[0,:] + b[4,1]*k[1,:] + b[4,2]*k[2,:] + b[4,3]*k[3,:], bb, p, interpolation)
                k[5,:] = dh*vecInt(xx + b[5,0]*k[0,:] + b[5,1]*k[1,:] + b[5,2]*k[2,:] + b[5,3]*k[3,:] + b[5,4]*k[4,:], bb, p, interpolation)

                xNew  = xx + c[0]*k[0,:]  + c[1]*k[1,:]  + c[2]*k[2,:]  + c[3]*k[3,:]  + c[4]*k[4,:]  + c[5]*k[5,:]
                xNewS = xx + cs[0]*k[0,:] + cs[1]*k[1,:] + cs[2]*k[2,:] + cs[3]*k[3,:] + cs[4]*k[4,:] + cs[5]*k[5,:]

                delta2 = np.dot((xNew-xNewS), (xNew-xNewS))
                delta = np.sqrt(delta2)

                if (delta2 > tol2):
                    dh = dh*(0.9*abs(tol/delta))**0.2
                    if (abs(dh) < hMin):
                        print "Error: step size underflow"
                        break
                else:
                    l += np.sqrt(np.sum((xx-xNew)**2))
                    xx = xNew                        
                    if (abs(dh) < hMin):
                        dh = 2*dh
                    sl += 1
                    self.tracers[sl,:] = xx
                    if ((dh > hMax) or (np.isnan(dh))):
                        dh = hMax
                    # check if this point lies outside the domain
                    if ((xx[0] < p.Ox-p.dx) or (xx[0] > p.Ox+p.Lx+p.dx) or (xx[1] < p.Oy-p.dy) or (xx[1] > p.Oy+p.Ly+p.dy) or (xx[2] < p.Oz) or (xx[2] > p.Oz+p.Lz)):
                        outside = True
                if ((dh > hMax) or (delta == 0) or (np.isnan(dh))):
                    dh = hMax
        
        self.tracers = np.resize(self.tracers, (sl, 3))
        self.l = l
        self.sl = sl
        self.p = p
        
        
def vecInt(xx, vv, p, interpolation = 'weighted'):
    """
    Interpolates the field around this position.
    
    call signature:
    
        vecInt(xx, vv, p, interpolation = 'weighted')
    
    Keyword arguments:
    
    *xx*:
      Position vector around which will be interpolated.
    
    *vv*:
      Vector field to be interpolated.
    
    *p*:
      Parameter struct.
    
    *interpolation*:
      Interpolation of the vector field.
      'mean': takes the mean of the adjacent grid point.
      'weighted': weights the adjacent grid points according to their distance.
    """
    
    # find the adjacent indices
    i  = (xx[0]-p.Ox)/p.dx
    if (i < 0):
        i = 0
    if (i > p.nx-1):
        i = p.nx-1
    ii = np.array([int(np.floor(i)), \
                    int(np.ceil(i))])
    
    j  = (xx[1]-p.Oy)/p.dy    
    if (j < 0):
        j = 0
    if (j > p.ny-1):
        j = p.ny-1
    jj = np.array([int(np.floor(j)), \
                    int(np.ceil(j))])
    
    k  = (xx[2]-p.Oz)/p.dz
    if (k < 0):
        k = 0
    if (k > p.nz-1):
        k = p.nz-1
    kk = np.array([int(np.floor(k)), \
                    int(np.ceil(k))])
    
    vv = np.swapaxes(vv, 1, 3)
    # interpolate the field
    if (interpolation == 'mean'):
        return np.mean(vv[:,ii[0]:ii[1]+1,jj[0]:jj[1]+1,kk[0]:kk[1]+1], axis = (1,2,3))
    if(interpolation == 'weighted'):
        if (ii[0] == ii[1]): w1 = np.array([1,1])
        else: w1 = (i-ii[::-1])
        if (jj[0] == jj[1]): w2 = np.array([1,1])
        else: w2 = (j-jj[::-1])
        if (kk[0] == kk[1]): w3 = np.array([1,1])
        else: w3 = (k-kk[::-1])
        weight = abs(w1.reshape((2,1,1))*w2.reshape((1,2,1))*w3.reshape((1,1,2)))
        return np.sum(vv[:,ii[0]:ii[1]+1,jj[0]:jj[1]+1,kk[0]:kk[1]+1]*weight, axis = (1,2,3))/np.sum(weight)
        

# class containing simulation parameters
class pClass:
    def __init__(self):
        self.dx = 0; self.dy = 0; self.dz = 0
        self.Ox = 0; self.Oy = 0; self.Oz = 0
        self.Lx = 0; self.Ly = 0; self.Lz = 0
        self.nx = 0; self.ny = 0; self.nz = 0
        
