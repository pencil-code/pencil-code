import pencil as pcn
import numpy as np

varfile = 'var.dat' # or specific snaphot as required 'VAR?'
var=pcn.read.var(varfile,magic=['tt'],trimall=True,quiet=True)
param=pcn.read.param(quiet=True)

filename='init_ism.in'
f = open(filename, 'w')
print(var.rho[:].size)
#smooth and ensure symmetric about midplane - assumes centred
#convert to cgs - so units can be applied independently
rho = (var.rho[:,0,0]+var.rho[::-1,0,0])/2*param.unit_density
tt  = (var.tt[:,0,0]+var.tt[::-1,0,0])/2*param.unit_temperature
#f.write('#--rho-------------------TT-------------\n')
for i in range(0,var.rho[:,0,0].size):
    f.write(str(rho[i])+'    '+str(tt[i])+'\n')
f.closed

def plot_ism(varfiles=[]):
    import matplotlib.pyplot as plt
    import glob
    fig1Dk=[15,3.9288]
    if len(varfiles) == 0:
        varfiles=glob.glob('data/proc0/VAR*')
        print(varfiles)
    for var in varfiles:
        varfile=var.split('/')[-1]
        globals()[varfile]=pcn.read.var(varfile,
                                       magic=['tt'],
                                       trimall=True,
                                       quiet=True
                                      )
    fig, ax = plt.subplots(1,3,figsize=fig1Dk)
    for ivar in range(0,len(varfiles)):
        for var in varfiles:
            if 'VAR'+str(ivar)==var.split('/')[-1]:
                varfile=var.split('/')[-1]
                ax[0].semilogy(globals()[varfile].z[:],globals()[varfile].rho[:,0,0],':',
                            )
                ax[1].semilogy(globals()[varfile].z[:],globals()[varfile].tt[:,0,0]*param.unit_temperature,':',
                            )
                ax[2].plot(globals()[varfile].z[:],globals()[varfile].uu[2,:,0,0],':',
                            )
    ax[0].semilogy(globals()[varfile].z,globals()[varfile].rho[:,0,0],'k-',label=str(round(globals()[varfile].t,2))+' Gyr')
    ax[1].semilogy(globals()[varfile].z,globals()[varfile].tt[:,0,0]*param.unit_temperature, 'k-',label=str(round(globals()[varfile].t,2))+' Gyr')
    ax[2].plot(globals()[varfile].z,globals()[varfile].uu[2,:,0,0], 'k-',label=str(round(globals()[varfile].t,2))+' Gyr')
    ax[0].set_ylabel(r'$\rho(z)$')
    ax[1].set_ylabel(r'$T(z)$')
    ax[2].set_ylabel(r'$u_z(z)$')
    ax[0].set_xlim([-2,2])
    ax[1].set_xlim([-2,2])
    ax[2].set_xlim([-2,2])
    #ax[2].set_ylim([globals()[varfile].uu[2,:,0,0].min(),globals()[varfile].uu[2,:,0,0].max()])
        #ax[imod,1].set_ylim(r'$p(u(k))$')
    ax[0].set_xlabel(r'$z$ [kpc]')        
    ax[1].set_xlabel(r'$z$ [kpc]')        
    ax[2].set_xlabel(r'$z$ [kpc]')        
    ax[0].legend(loc='upper left',framealpha=0.5)
    ax[1].legend(loc='lower left',framealpha=0.5)
    ax[2].legend(loc='upper right',framealpha=0.5)
    fig.tight_layout()
    plt.show() 

plot_ism()
