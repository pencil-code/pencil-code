import pencil as pc
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable

ff=pc.read_var(trimall=True,ivar=3)
dim=pc.read_dim()

rad,phi=np.meshgrid(ff.x,ff.y)

fig,((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3,figsize=(10,5))

ax1.set_aspect("equal")

ax1.contourf(rad*cos(phi),rad*sin(phi),ff.rho[0,...],256)

ax2.contourf(rad,phi,ff.rho[0,...],256)

ax3.plot(ff.x,mean(ff.rho[0,...],axis=0))

ip1=39
ip2=40

ax4.plot(ff.y/pi,.5*(ff.rho[0,:,ip1]+ff.rho[0,:,ip2]))

# planet is at pi=ff.t. shift it to zero
phi = ff.t
while phi > 2*pi:
    phi = phi - 2*pi
print phi
phi_opposite = phi + pi 
if (phi_opposite > pi):
    phi_opposite = phi_opposite - 2*pi
if (phi_opposite < -pi):
    phi_opposite = phi_opposite + 2*pi

dphi = ff.y[1]-ff.y[0]
iphi = np.rint((phi_opposite + pi)/dphi)
if (iphi==dim.ny):
    iphi=0

#ax4.plot(ff.x,ff.rho[0,iphi,:])

ts=pc.read_ts()
ax6.plot(ts.t,ts.torqint_1)
ax6.plot(ts.t,ts.torqext_1)
ax6.plot(ts.t,ts.torqint_1+ts.torqext_1)


plt.show()

