import numpy as np
import pencil as pc
import os
import re
import h5py
import matplotlib.pyplot as plt
from plot2d import plot_cart

def get_profiles(nt, endt, sn, hf
                ):
        
    E0=sn.EE
    rho0=sn.rho
    xi=2.026 #McKee/Ostriker 1988
    sedov_time=sn.t_sedov
    time   =np.linspace(0,endt,nt)+sedov_time
    M0=10 # (solar masses)
    m_u=1.66053886
    n0=rho0/m_u
    vej=(400*E0/M0)**0.5
    
    #Woltier transition to momentum conservation
    vtrans=2.e2*n0**(2./17.)*E0**(1./17.)
    ttrans=(2./(2.+dims)*(xi*E0/rho0)**(1./(dims+2.))*
                  vtrans**(-1))**(1./(1-2./(dims+2.)))
    rtrans=(xi*E0/rho0)**(1./(dims+2.))*ttrans**(2./(dims+2.))
    
    #pressure driven snowplough transition
    tpds=(3.61e-5*E0**(3./14.)/n0**(4./7.))/2.718281828
    rpds=(xi*E0/rho0)**(1./(dims+2.))*tpds**(2./(dims+2.))
    vpds=(2./(2+dims))*(xi*E0/rho0)**(2./(2+dims)-1)
    
    #momentum conserving snowplough transition
    tmcs=61*vej**3/n0**(3./7.)/E0**(3./14.)*tpds
    rmcs=rpds*(4./3.*tmcs/tpds-1./3.)**(3./10.)
    
    print 'tpds, tmcs, ttrans, endt=',tpds, tmcs, ttrans, endt
    #Sedov-Taylor
    rst = (xi*E0/rho0)**(1./(dims+2.))*time**(2./(dims+2.))
    vst = 2./(dims+2.)*(xi*E0/rho0)**(1./(dims+2.))*time**(2./(dims+2.)-1)
    #snowplough
    rsnpl=np.zeros(nt)
    vsnpl=np.zeros(nt)
    rsnpl[:]=rst
    vsnpl[:]=vst
    isnpl=np.where(time>ttrans)[0]
    rsnpl[isnpl]=rtrans*(8./(dims+2.)*time[isnpl]/ttrans-3./(dims+2.))**(1./4.)
    vsnpl[isnpl]=vtrans*(8./(dims+2.)*time[isnpl]/ttrans-3./(dims+2.))**(-3./4.)
    #Cioffi et al 
    rcioffi=np.zeros(nt)
    vcioffi=np.zeros(nt)
    rcioffi[:]=rst
    vcioffi[:]=vst
    icioffi=np.where(time>tpds)[0]
    rcioffi[icioffi]=rpds*(4./3.*time[icioffi]/tpds-1./3.)**(3./10.)
    vcioffi[icioffi]=0.3*rpds*4./3./tpds*(4./3.*time[icioffi]/tpds
                                         -1./3.)**(-7./10.)
    jcioffi=np.where(time>tmcs)[0]
    rcioffi[jcioffi]=rpds*(4.66*(time[jcioffi]-tmcs)/tpds
                               *(1.-0.779/(tmcs/tpds)**0.17)+(rmcs/rpds)**4
                          )**0.25
    vcioffi[jcioffi]=0.25*rpds*4.66/tpds*(1.-0.779/(tmcs/tpds)**0.17)*(
                               4.66*(time[jcioffi]-tmcs)/tpds
                              *(1.-0.779/(tmcs/tpds)**0.17)+(rmcs/rpds)**4
                               )**(-0.75)
    #save analytic profiles to file
    if not hf.__contains__('analytic'):
        grp=hf.create_group('analytic')
    else:
        grp=hf['analytic']
    grp.create_dataset('time',    (nt,),    data=time)
    agrp=grp.create_group('sedov-taylor')
    bgrp=grp.create_group('snowplough')
    cgrp=grp.create_group('cioffi')
    agrp.create_dataset('radius',  (nt,),    data=rst)
    agrp.create_dataset('speed',   (nt,),    data=vst)
    bgrp.create_dataset('radius',  (nt,),    data=rsnpl)
    bgrp.create_dataset('speed',   (nt,),    data=vsnpl)
    cgrp.create_dataset('radius',  (nt,),    data=rcioffi)
    cgrp.create_dataset('speed',   (nt,),    data=vcioffi)
    

nt=2000
home=os.path.expanduser('~')
datatopdir=os.getcwd()
models=['sedov']
imod = len(models)
os.chdir('data')
figsdir = os.getcwd() 
figsdir = re.sub('\/data\/*$','',figsdir) + '/video_slices/' # name for dir saving figures
if not os.path.exists(figsdir):
    os.makedirs(figsdir)
os.chdir(datatopdir) 
sn=pc.read_sn()
sedov_time=sn.t_sedov
f=open('data/tsnap.dat','r')
nvar= int(str.rsplit(f.readline())[1])
#nvar= 12
print nvar
var=pc.read_var(ivar=nvar-1,quiet=True,proc=0) 
endt = var.t
param=pc.read_param(quiet=True)
tokms = param.unit_velocity/1e5

dim=pc.read_dim()
dims=3
if dim.nxgrid==1:
    dims -= 1
if dim.nygrid==1:
    dims -= 1
if dim.nzgrid==1:
    dims -= 1
print 'dims ', dims

hf = h5py.File(datatopdir+'/data/'+models[imod]+'_sedov.h5', 'w')
get_profiles(nt, endt, sn, hf)
lvar=True
nx=dim.nxgrid/dim.nprocx
print 'nx =',nx
proc=dim.nprocx*dim.nprocy*dim.nprocz/2+dim.nprocy/2
#print 'proc {}, ix {}, iy {}, iz {}'.format(proc, ix, iy, iz)

TT      =np.empty([nvar,nx])
rho     =np.empty([nvar,nx])
ux      =np.empty([nvar,nx])
pp      =np.empty([nvar,nx])
shock   =np.empty([nvar,nx])
netheat =np.empty([nvar,nx])
cooling =np.empty([nvar,nx])

dsnap   =0.000025
radius  =np.empty(nvar)
rspeed  =np.empty([1,1,nvar+6])
speed   =np.empty(nvar)
ntime   =np.empty(nvar)
cspeed  =np.empty(nvar)
maxshock=np.empty(nvar)
for ivar in range(0,nvar):
    if lvar:
        varfile='VAR'+"%i"%ivar
        var=pc.read_var(varfile,
                        trimall=True,quiet=True,
                        proc=proc,
                        magic=['tt','pp']
                       )
        izz, iyy, ixx = np.array(np.where(var.rho==var.rho.max()))
        rtmp=np.sqrt(var.x[ixx]**2+var.y[iyy]**2+var.z[izz]**2).mean()
        if rtmp>var.x[nx-2]:
            oldproc = proc
            newproc = proc+1
            vartry=pc.read_var(varfile,
                            trimall=True,quiet=True,
                            proc=newproc,
                            magic=['tt','pp']
                           )
            if vartry.rho.max()>var.rho.max():
                var=vartry
                proc=newproc
                izz, iyy, ixx = np.array(np.where(var.rho==var.rho.max()))
                rtmp=np.sqrt(var.x[ixx]**2+var.y[iyy]**2+var.z[izz]**2).mean()
            else:
                proc=oldproc
        radius[ivar]=rtmp
        print 'loaded '+varfile+', radius {}, time {}'.format(radius[ivar], 
                                                            var.t+sedov_time)
        izz, iyy, ixx = np.array(np.where(var.shock==var.shock.max()))
        maxshock[ivar]=np.sqrt(var.x[ixx]**2+var.y[iyy]**2+var.z[izz]**2).mean()
        ntime  [ivar]  =var.t + sedov_time
        TT     [ivar,:]=var.tt[0,0,:]
        rho    [ivar,:]=var.rho[0,0,:]
        ux     [ivar,:]=var.ux[0,0,:]
        pp     [ivar,:]=var.pp[0,0,:]
        shock  [ivar,:]=var.shock[0,0,:]
        netheat[ivar,:]=var.netheat[0,0,:]
        cooling[ivar,:]=var.cooling[0,0,:]
        var=pc.read_var(varfile,
                        trimall=True,quiet=True,
                        magic=['tt','pp']
                       )
        cspeed[ivar]=np.sqrt(param.cp*var.tt[-1,-1,-1]*(param.gamma-1))
        if np.mod(ivar,10)==0:
            vidfile = 'density_SN'+"%.3i"%ivar+'.png'
            cbar_label=r'$ \rho $ [cm$^{-3}$]'
            vslice=1/1.6728*var.rho[dim.nzgrid/2]
            plot_cart(vslice,
                      var.x,var.y, figsdir+vidfile,
                      cbar_label=cbar_label,
                      x_unit=' [kpc]', y_unit=' [kpc]',
                      time_stamp="%.2f"%1000*ntime[ivar], t_unit=' Myr',
                      cmin=vslice.min(), cmax=vslice.max()
                     ) 
            vidfile = 'temperature_SN'+"%.3i"%ivar+'.png'
            cbar_label=r'$T$ [10$^6$K]'
            vslice=param.unit_temperature/1e6*var.tt[dim.nzgrid/2]
            plot_cart(vslice,
                      var.x,var.y, figsdir+vidfile,
                      cbar_label=cbar_label,
                      x_unit=' [kpc]', y_unit=' [kpc]',
                      time_stamp="%.2f"%1000*ntime[ivar], t_unit=' Myr',
                      cmin=vslice.min(), cmax=vslice.max()
                     ) 
if not os.path.exists('data/proc0/VAR'+"%i"%(ivar+1)): #avoid error after last var file
    lvar=False
else:
    lvar=True
maxshock[0]=radius[0]
dt=(ntime[-1]-ntime[0])/(nvar-1)
rspeed[0,0,3:-3]=radius
for i in range(0,3):
    rspeed[0,0,2-i] = 2*rspeed[0,0,3-i]-rspeed[0,0,4-i]
    rspeed[0,0,nvar+3+i] = 2*rspeed[0,0,nvar+2+i]-rspeed[0,0,nvar+1+i]
print 'rspeed[0,0,:6] ',rspeed[0,0,:6]
print 'rspeed[0,0,-6:] ',rspeed[0,0,-6:]
aspeed = pc.xder(rspeed,dt)[0,0,:]
for i in range(0,nvar):
    ismooth = np.array(np.where(np.abs(np.arange(nvar+6)-i-3)<(2+i/5)))
    speed[i] = aspeed[ismooth].mean()

if not hf.__contains__('numeric'):
    ngrp=hf.create_group('numeric')
else:
    ngrp=hf['numeric']
mgrp=ngrp.create_group(models[imod])
mgrp.create_dataset('time',          (nvar,),    data=ntime   )
mgrp.create_dataset('radius',        (nvar,),    data=radius  )
mgrp.create_dataset('speed',         (nvar,),    data=speed *tokms  )
mgrp.create_dataset('sound-speed',   (nvar,),    data=cspeed*tokms  )
mgrp.create_dataset('max-shock',     (nvar,),    data=maxshock)
mgrp.create_dataset('temperature',   (nvar,nx),  data=TT      )
mgrp.create_dataset('density',       (nvar,nx),  data=rho     )
mgrp.create_dataset('velocity',      (nvar,nx),  data=ux      )
mgrp.create_dataset('pressure',      (nvar,nx),  data=pp      )
mgrp.create_dataset('shock',         (nvar,nx),  data=shock   )
mgrp.create_dataset('net-cooling',   (nvar,nx),  data=netheat )
mgrp.create_dataset('cooling',       (nvar,nx),  data=cooling )
mgrp.create_dataset('x',             (nx,),      data=var.x   )
    
    
figname='sedov-taylor_radius.png'
plt.figure()
plt.plot(
         hf['analytic']['time'][:],
         hf['analytic']['sedov-taylor']['radius'][:], 'k:',
                  label='sedov-taylor'
        )
plt.plot(
         hf['analytic']['time'][:],
         hf['analytic']['snowplough']['radius'][:], 'b:',
                  label='snowplough'
        )
plt.plot(
         hf['analytic']['time'][:],
         hf['analytic']['cioffi']['radius'][:], 'g:',
                  label='cioffi'
        )
for imod in range(0,len(models)):
    plt.plot(
             hf['numeric'][models[imod]]['time'],
             hf['numeric'][models[imod]]['radius'], '--',
                     label=models[imod]+'radius'
        )
plt.plot(
         hf['numeric'][models[imod]]['time'],
         hf['numeric'][models[imod]]['max-shock'], 'c-.',
                 label=models[imod]+'max-shock'
        )
plt.legend(loc='lower right')
plt.savefig(figname)

figname='sedov-taylor_speed.png'
plt.figure()
plt.plot(
         hf['analytic']['time'],
         hf['analytic']['sedov-taylor']['speed'], 'k:',
               label='sedov-taylor'
        )
plt.plot(
         hf['analytic']['time'],
         hf['analytic']['snowplough']['speed'], 'b:',
               label='snowplough'
        )
plt.plot(
         hf['analytic']['time'],
         hf['analytic']['cioffi']['speed'], 'g:',
               label='cioffi'
        )
for imod in range(0,len(models)):
    plt.plot(
             hf['numeric'][models[imod]]['time'],
             hf['numeric'][models[imod]]['speed'], '--',
                  label=models[imod]+'speed'
        )
plt.plot(
         hf['numeric'][models[imod]]['time'],
         hf['numeric'][models[imod]]['sound-speed'], 'c:',
              label=models[imod]+'sound-speed'
        )                
plt.legend(loc='upper right')
plt.gca().set_yscale('log',subsy=[5,10])
plt.savefig(figname)

hf.close
