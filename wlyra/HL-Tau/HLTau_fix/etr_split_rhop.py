def find_index_bisect(qpar,q):

    jl=0
    ju=len(q)-1
    
    while ((ju-jl)>1):
        jm=(ju+jl)/2
        if (qpar > q[jm]):
            jl=jm
        else:
            ju=jm

    if (qpar-q[jl] <= q[ju]-qpar):
        iq0=jl
    else:
        iq0=ju

    return iq0
#
def particles_to_density(xxp,yyp,zzp,dim,x,y,z): 

    npar=len(xxp)

    nnp = np.zeros([dim.mz,dim.my,dim.mx])
    
    for k in range(npar):

        xp=xxp[k]
        yp=yyp[k]
        zp=zzp[k]
        
        ix0=find_index_bisect(xp,x)
        iy0=find_index_bisect(yp,y)
        iz0=find_index_bisect(zp,z)
        
        ixx0=ix0-1
        ixx1=ix0+1

        iyy0=iy0-1
        iyy1=iy0+1

        izz0=iz0-1
        izz1=iz0+1

        if (dim.nx > 1):
            dx=x[ix0]-x[ixx0]
        else:
            dx=1.0
           
        if (dim.ny > 1):
            dy=y[iy0]-y[iyy0]
        else:
            dy=1.0
           
        if (dim.nz > 1):
            dz=z[iz0]-z[izz0]
        else:
            dz=1.0

        dx_1=1.0/dx
        dy_1=1.0/dy
        dz_1=1.0/dz
#       
        dx_2=1.0/dx**2
        dy_2=1.0/dy**2
        dz_2=1.0/dz**2
        
        for ixx in np.arange(ixx0,ixx1+1):
            for iyy in np.arange(iyy0,iyy1+1):
                for izz in np.arange(izz0,izz1+1):

                    if ( ((ixx-ix0) == -1) or ((ixx-ix0) == +1) ):
                        weight_x = 1.125 - 1.5*abs(xp-x[ixx])  *dx_1 + 0.5*abs(xp-x[ixx])**2*dx_2
                    else:
                        if (dim.nx != 1): weight_x = 0.75 - (xp-x[ixx])**2*dx_2

                    if ( ((iyy-iy0) == -1) or ((iyy-iy0) == +1) ):
                        weight_y = 1.125 - 1.5*abs(yp-y[iyy])  *dy_1 + 0.5*abs(yp-y[iyy])**2*dy_2
                    else:
                        if (dim.ny != 1): weight_y = 0.75 - (yp-y[iyy])**2*dy_2

                    if ( ((izz-iz0) == -1) or ((izz-iz0) == +1) ):
                        weight_z = 1.125 - 1.5*abs(zp-z[izz])  *dz_1 + 0.5*abs(zp-z[izz])**2*dz_2
                    else:
                        if (dim.nz != 1): weight_z = 0.75 - (zp-z[izz])**2*dz_2

                    weight=1.0
                    if (dim.nx != 1): weight=weight*weight_x
                    if (dim.ny != 1): weight=weight*weight_y
                    if (dim.nz != 1): weight=weight*weight_z
#
                    nnp[izz,iyy,ixx]=nnp[izz,iyy,ixx] + weight
    
    return nnp            
#
#---------------------------------------------------------
#---------------------------------------------------------
#---------------------------------------------------------
#
import pencil as pc
import matplotlib.pyplot as plt 
import numpy as np
from pylab import *
import sys
#
fp=pc.read_pvar()
ff=pc.read_var()
par=pc.read_param()
#
pdim = pc.read_pdim()
dim = pc.read_dim()
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

z1=[]
z2=[]
z3=[]
z4=[]

x=ff.x
y=ff.y
z=ff.z

for k in range(npar):
    kk=ipar[k]
    xp=fp.xp[k]
    yp=fp.yp[k]
    zp=fp.zp[k]
    if (kk <= nsplit-1):
        #print 'first species',ipar[k],fp.xp[k]
        x1.append(xp)
        y1.append(yp)
        z1.append(zp)
    if (kk >   nsplit and kk <= 2*nsplit-1):
        #print 'second species',ipar[k],fp.xp[k]
        x2.append(xp)
        y2.append(yp)
        z2.append(zp)
    if (kk > 2*nsplit and kk <= 3*nsplit-1):
        #print 'third species',ipar[k],fp.xp[k]
        x3.append(xp)
        y3.append(yp)
        z3.append(zp)
    if (kk > 3*nsplit and kk <= 4*nsplit-1):
        #print 'fourth species',ipar[k],fp.xp[k]
        x4.append(xp)
        y4.append(yp)
        z4.append(zp)
#

rhop1 = par.rhop_swarm * particles_to_density(x1,y1,z1,dim,x,y,z)
rhop2 = par.rhop_swarm * particles_to_density(x2,y2,z2,dim,x,y,z)
rhop3 = par.rhop_swarm * particles_to_density(x3,y3,z3,dim,x,y,z)
rhop4 = par.rhop_swarm * particles_to_density(x4,y4,z4,dim,x,y,z)




epsi=1e-4
#
fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(12,12))

ncolors=10

rad = ff.x
theta = ff.y
rad2d,theta2d = np.meshgrid(rad,theta)

x2d = rad2d*np.cos(theta2d)
y2d = rad2d*np.sin(theta2d)

ax1.contourf(x2d,y2d,rhop1[0,:,:],ncolors)
ax2.contourf(x2d,y2d,rhop2[0,:,:],ncolors)
ax3.contourf(x2d,y2d,rhop3[0,:,:],ncolors)
ax4.contourf(x2d,y2d,rhop4[0,:,:],ncolors)

ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')
ax4.set_aspect('equal')

plt.show()

