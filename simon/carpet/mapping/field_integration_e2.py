# fiel_integration_e2.py
#
# Integraates queentities along field lines for e2.

import numpy as np
import pylab as plt
import multiprocessing as mp
import time as tm


# Create the color mapping.
def mapping(xi, yi, xf, yf):
    nx = xi.shape[0]
    ny = yi.shape[1]
    m = np.zeros((nx, ny, 3))
    for i in range(nx):
        for j in range(ny):
            if ((xf[i,j] > xi[i,j]) and (yf[i,j] > yi[i,j])):
                m[i,j,:] = [0,1,0]
            if ((xf[i,j] > xi[i,j]) and (yf[i,j] < yi[i,j])):
                m[i,j,:] = [1,1,0]
            if ((xf[i,j] < xi[i,j]) and (yf[i,j] > yi[i,j])):
                m[i,j,:] = [0,0,1]
            if ((xf[i,j] < xi[i,j]) and (yf[i,j] < yi[i,j])):
                m[i,j,:] = [1,0,0]
    return m


# Integrate JJ along the filed line.
def jpInt(jj, bb, s, p):
    jp = 0
    
    for i in range(s.tracers.shape[0]-1):
        dx = (s.tracers[i+1,:] - s.tracers[i,:])
        ds = np.sqrt(np.dot(dx, dx))
        xx = (s.tracers[i+1,:] + s.tracers[i,:])/2
        jInt = pc.vecInt(xx, jj, p)
        #bInt = pc.vecInt(xx, bb, p)
        #jp += np.dot(jInt, bInt)*ds/np.sqrt(np.dot(bInt, bInt))
        jp += np.dot(jInt, dx)
        
    return jp


# Fiel line integration.
def intField(q, x, y, z, bb, jj, p, gx, gy, gz):
    jp = np.zeros((len(x), len(y)))
    xi = np.zeros((len(x), len(y)))
    yi = np.zeros((len(x), len(y)))
    xf = np.zeros((len(x), len(y)))
    yf = np.zeros((len(x), len(y)))

    #tm.sleep(np.random.random(1))
    
    for i in range(len(x)):
        for j in range(len(y)):
            xyz0 = np.array([x[i], y[j], z])
            s = pc.stream(bb[:,gz:p.nz+gz,gy:p.ny+gy,gx:p.nx+gx], p, xx = xyz0)
            jp[i, j] = jpInt(jj, bb, s, p)
            #xi[i, j] = x[i]
            #yi[i, j] = y[j]
            xi[i, j] = s.tracers[0,0].copy()
            yi[i, j] = s.tracers[0,1].copy()
            xf[i, j] = s.tracers[-1,0].copy()
            yf[i, j] = s.tracers[-1,1].copy()
            
    q.put((jp, xi, yi, xf, yf))


# Domain parameters.
gx = 3; gy = 3; gz = 3 # number of ghost cells
p = pc.pClass()
p.Lx = 8.; p.Ly = 8.; p.Lz = 32.
p.Ox = -4.; p.Oy = -4.; p.Oz = -16.
p.nx = 64; p.ny= p.nx; p.nz = 4*p.nx
p.dx = p.Lx/p.nx; p.dy = p.Ly/p.ny; p.dz = p.Lz/p.nz

# Domain arrays.
x = np.linspace(p.Ox-p.dx*gx, p.Ox+p.Lx+p.dx*gx, p.nx+2*gx)
y = np.linspace(p.Oy-p.dy*gy, p.Oy+p.Ly+p.dy*gy, p.ny+2*gy)
z = np.linspace(p.Oz-p.dz*gz, p.Oz+p.Lz+p.dz*gz, p.nz+2*gz)
xyz = np.meshgrid(x, y, z)
xx = xyz[0].swapaxes(0,1)
yy = xyz[1].swapaxes(0,1)
zz = xyz[2].swapaxes(0,1)
xi = np.zeros((len(x)-2*gx, len(y)-2*gy))
yi = np.zeros((len(x)-2*gx, len(y)-2*gy))
xf = np.zeros((len(x)-2*gx, len(y)-2*gy))
yf = np.zeros((len(x)-2*gx, len(y)-2*gy))
jp = np.zeros((len(x)-2*gx, len(y)-2*gy))

# Field parameters.
B0 = 1
e2_xc = np.array([1, -1, 1, -1], dtype = 'float32')
e2_yc = np.array([0, 0, 0, 0], dtype = 'float32')
e2_zc = np.array([-12, -4, 4, 12], dtype = 'float32')
e2_k = np.array([1, -1, 1, -1], dtype = 'float32')
e2_a = np.sqrt(np.array([2, 2, 2, 2], dtype = 'float32'))
e2_l = np.array([2, 2, 2, 2], dtype = 'float32')

# Proc parameters.
nproc = 4

# Create the magnetic field and compute the current.
bb = np.zeros((3, p.nx+2*gx, p.ny+2*gy, p.nz+2*gz))
for i in range(len(e2_xc)):
    bb[0,...] += 2*B0*e2_k[i]*np.exp((-(xx-e2_xc[i])**2-(yy-e2_yc[i])**2)/e2_a[i]**2-(zz-e2_zc[i])**2/e2_l[i]**2) \
                 *(-(yy-e2_yc[i]))
    bb[1,...] += 2*B0*e2_k[i]*np.exp((-(xx-e2_xc[i])**2-(yy-e2_yc[i])**2)/e2_a[i]**2-(zz-e2_zc[i])**2/e2_l[i]**2) \
                 *(xx-e2_xc[i])
             
# Add the background field.
bb[2,...] += B0
# Pencil code array structure izyx
bb = bb.swapaxes(1, 3)
# Compute the electric current density.
jj = pc.curl(bb, p.dx, p.dy, p.dz)


# Find the field line mapping and the integrated quantities.
queue = [mp.Queue() for iproc in range(nproc)]
tmp = []
#inFieldLambda = lambda queue, 
proc = []
for iproc in range(nproc):
    proc.append(mp.Process(target = intField, args = (queue[iproc], x[gx+iproc:p.nx+gx:nproc], y[gy:p.ny+gy], z[gz], bb, jj, p, gx, gy, gz)))
for iproc in range(nproc):
    proc[iproc].start()
for iproc in range(nproc):
    tmp.append(queue[iproc].get())
for iproc in range(nproc):
    proc[iproc].join()
for iproc in range(nproc):
    jp[iproc::nproc,:] = tmp[iproc][0].copy()
    xi[iproc::nproc,:] = tmp[iproc][1].copy()
    yi[iproc::nproc,:] = tmp[iproc][2].copy()
    xf[iproc::nproc,:] = tmp[iproc][3].copy()
    yf[iproc::nproc,:] = tmp[iproc][4].copy()
for iproc in range(nproc):
    proc[iproc].terminate()
    
#for i in range(gx, len(x)-gx):
    #for j in range(gy, len(y)-gy):
        #xyz0 = np.array([x[i], y[j], z[gz]])
        #s = pc.stream(bb[:,gz:p.nz+gz,gy:p.ny+gy,gx:p.nx+gx], p, xx = xyz0)
        #jp[i-gx, j-gy] = jpInt(jj, bb, s, p)
        #xi[i-gx, j-gy] = s.tracers[0,0].copy()
        #yi[i-gx, j-gy] = s.tracers[0,1].copy()
        #xf[i-gx, j-gy] = s.tracers[-1,0].copy()
        #yf[i-gx, j-gy] = s.tracers[-1,1].copy()

# Compute the color map.
m = mapping(xi, yi, xf, yf)

# Save arrays into file.
np.save('e2_jp.dat', jp)
np.save('e2_xi.dat', xi)
np.save('e2_yi.dat', yi)
np.save('e2_xi.dat', xf)
np.save('e2_yf.dat', yf)
np.save('e2_m.dat', m)
