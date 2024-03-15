@parameters
@data/index
@data/pc_constants
!p.multi=0
;
default,ivar,-1
default,iread,0
if iread eq 0 then begin
  pc_read_dim,obj=dim
  pc_read_param,/param2,obj=param
  if ivar lt 0 then pc_read_var,/bb,/trimall,obj=var else pc_read_var,/bb,/trimall,obj=var,ivar=ivar
  ;if ivar lt 0 then pc_read_var,/trimall,obj=var else pc_read_var,/trimall,obj=var,ivar=ivar
  xxx=var.x
  yyy=var.y
  zzz=var.z
;  B0=param.B_ext(2)
  iread=1
endif
nz=n_elements(zzz)
;
stop
!p.multi=[0,1,3]
!p.charsize=2
;
lnTT=haver(var.lnTT)
lnrho=haver(var.lnrho)
sss=haver(var.sss)
ppp=haver(var.ppp)
rho=exp(haver(var.lnrho))
TT=exp(haver(var.lnTT))
kapparho=haver(var.kapparho)
yh=haver(var.yh)
uz=haver(var.uu(*,2))
;
;plot,zzz,uz,yr=[-1,1]*10.,ytit='uz',ps=-1
;plot_io,zzz,rho,yr=[1e-10,2e-3],ytit='rho'
;plot_io,zzz,sss,ytit='!8s!6'
;plot_io,zzz,kapparho,ytit="kapparho"
;plot,zzz,yH,ytit="!8y!6!dH!n",yr=[0,1]
;
;  compute tau
;
circ_sym,1.3,1
tau_top=param.tau_top
tau=var.kapparho*0.
tau=integr(var.kapparho,/rev,x=var.z)+tau_top
iz_tau1=findex(1.,tau,/rev)
iz_tau1=iz_tau1<nz-1
print,'iz_tau1=',iz_tau1
zzz_tau1=zzz(iz_tau1)
TT_tau1=TT(iz_tau1)
plot,zzz,TT,ytit="!8T!6 [K]",yst=0
oplot,zzz_tau1*[1,1],TT_tau1*[1,1],ps=8,col=122,thick=5
;
plot_io,zzz,rho,ytit="!7q!6 [g/cm!u3!n]",yst=0
oplot,zzz(iz_tau1)*[1,1],rho(iz_tau1)*[1,1],ps=8,col=122,thick=5
;
;  radiative nabla
;
cp=var.cp
gam=var.gamma
grav=-param.gravz
frad=var.kr_frad(*,2)/var.kapparho
Teff=(frad/sigmaSB)^.25
K=16.*sigmaSB*TT^3/(3.*var.kapparho)
nab_ad=1.-1./gam
nab_rad=nab_ad*cp*frad/(grav*K)
plot_io,zzz,nab_rad,yr=[.1,max(nab_rad)*1.5],ytit='!9G!6' & oplot,zzz,zzz*0+nab_ad
if n_elements(inabad) ne 0 then oplot,zzz,var.nabad,col=122,thick=5
if n_elements(inabad) ne 0 then nab_ad=var.nabad
;
;  superadiabatic gradient
;
;save,file='var.sav',zzz,rho,TT,kapparho,yh,tau,sss
;
!p.multi=0
print,'write stratification file?'
stop
file='stratification.dat'
wtable3,file,zzz,lnrho,lnTT
;device, /close
;
;  compute effective KK and effective Prandtl number
;
print,'compute chi?'
stop
cp=haver(var.cp)
gam=haver(var.gamma)
KK=16.*sigmaSB*TT^3/(3.*var.kapparho)
chi=KK/(cp*rho)
nu=param.nu
Pr=nu/chi
kappa_cst=param.kappa_cst
save,file='var.sav',zzz,rho,TT,kapparho,yh,tau,sss,cp,gam,chi,nab_ad,nab_rad,kappa_cst,zzz_tau1,TT_tau1,frad
plot_io, zzz, KK
;
xHe=param.xHe
YY=1./(1.+1./(4.*param.xHe))
muY=1./(1.-YY)
mu=muY/(1.+var.yH+xHe)
Hp=var.cp*exp(var.lnTT)/(-mu*param.gravz)
;
print,'write stratification file?'
stop
file='stratification.dat'
wtable3,file,zzz,lnrho,lnTT
;
itau1=findex(1.,tau)
Teff=(reform(TT(itau1)))(0)
flux=var.kr_frad(*,2)/var.kapparho
flux_top=sigmaSB*Teff^4
grav=-param.gravz
nab_ad=1.-1./gam
nab_rad=nab_ad*cp*flux/grav*KK
help,nab_ad,cp,flux,grav,KK
print,Teff,flux
!p.multi=0
END
