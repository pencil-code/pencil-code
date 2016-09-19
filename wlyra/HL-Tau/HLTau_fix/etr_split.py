import pencil as pc
import matplotlib.pyplot as plt 
import numpy as np
from pylab import *
import sys
#
fp=pc.read_pvar()
#
pdim = pc.read_pdim()
#
# ff=pc.read_var(trimall=True)
# dim=pc.read_dim()
# rho = ff.rho[0,0:dim.ny,0:dim.nx]
# rhop=ff.rhop[0,0:dim.ny,0:dim.nx]
#
npar_species = 4 
#
npar = len(fp.xp)
#
nsplit = pdim.npar/npar_species
#
ipar = fp.ipars
#
x1=[]
x2=[]
x3=[]
x4=[]
y1=[]
y2=[]
y3=[]
y4=[]
for k in range(npar):
    kk=ipar[k]
    if (kk <= nsplit-1):
        #print 'first species',ipar[k],fp.xp[k]
        x1.append(fp.xp[k])
        y1.append(fp.yp[k])
    if (kk >   nsplit and kk <= 2*nsplit-1):
        x2.append(fp.xp[k])
        y2.append(fp.yp[k])
        #print 'second species',ipar[k],fp.xp[k]
    if (kk > 2*nsplit and kk <= 3*nsplit-1):
        x3.append(fp.xp[k])
        y3.append(fp.yp[k])
        #print 'third species',ipar[k],fp.xp[k]  
    if (kk > 3*nsplit and kk <= 4*nsplit-1):
        x4.append(fp.xp[k])
        y4.append(fp.yp[k])
        #print 'fourth species',ipar[k],fp.xp[k]

#
epsi=1e-4
#
fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,12))

ax1.plot(x1*cos(y1),x1*sin(y1),',')
ax2.plot(x2*cos(y2),x2*sin(y2),',')
ax3.plot(x3*cos(y3),x3*sin(y3),',')
ax4.plot(x4*cos(y4),x4*sin(y4),',')

ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')
ax4.set_aspect('equal')

#ax3.plot(a,e,'.')
#ax3.set_xlim([0.4,2.5])
#ax3.set_ylim([0.,1.])
plt.show()
