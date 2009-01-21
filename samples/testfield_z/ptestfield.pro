;$Id$
;
;  plot alpa and eta results from time series
;
pc_read_ts,o=o
;
;  set minimum time after which averaging begins
;
t2=1e30
spawn,'touch parameters.pro'
@parameters
@data/index.pro
@data/testfield_info.dat
default,alpmax,.05
default,gammax,.05
default,delmax,.05
default,etamax,.01
default,run,''
;
;  introduce abbreviations
;
tt=o.t
urms=o.urms
;
;  alpha tensor
;
alp11=o.alp11
alp21=o.alp21
;
;  eta tensor as in EMF_i = ... -eta_ij J_j
;  so these new etas are also referred to as eta^*
;
eta12=-o.eta11
eta22=-o.eta21
;
;  read extra fields (if itestfield eq 'B11-B22'
;  as opposed to just itestfield eq 'B11-B21')
;
if itestfield eq 'B11-B22' then begin
  alp12=o.alp12
  alp22=o.alp22
  eta11=+o.eta12
  eta21=+o.eta22
endif
;
;  range of time where to do the analysis
;
tmin=min(tt)
tmax=max(tt)
default,t1,(tmin+tmax)/2.
good=where(tt gt t1 and tt lt t2)
kf=5.
;
!p.charsize=1.6
!p.multi=[0,2,2]
!x.title='!6'
;
;  give modified alpmax values in parameters.pro file
;
yralp=[-.2,1.]*alpmax
plot,tt,-alp11,yr=yralp,ytit='!7-a!6'
if itestfield eq 'B11-B22' then begin
  oplot,tt,-alp22,li=2
  alp=.5*(alp11+alp22)
endif else begin
  alp=alp11
endelse
oplot,tt,urms/3.,li=1
pc_error_range,tt(good),-alp(good),/oplot,mean=alpm,error=alp_error
oplot,tt(good),-accum(alp(good)),col=55
oplot,tt,tt-tt,li=3
;
;  off-diagonals of alpha tensor
;
yrgam=[-.5,1.]*gammax
plot,tt,-alp12,yr=yrgam,ytit='!7c!6'
if itestfield eq 'B11-B22' then begin
  gam=-.5*(alp12-alp21)
endif else begin
  gam=-.5*alp12
endelse
oplot,tt,urms/3.,li=1
pc_error_range,tt(good),gam(good),mean=gamm,error=gam_error
pc_error_range,tt(good),-alp12(good),/oplot,mean=alp12m,error=alp12_error
pc_error_range,tt(good),+alp21(good),/oplot,mean=alp21m,error=alp21_error
oplot,tt(good),-accum(alp12(good)),col=55
oplot,tt(good),+accum(alp21(good)),col=55
oplot,tt,tt-tt,li=3
;
;  off-diagonals of eta tensor
;  Warning: plot here +eta*_12 and -eta*_21.
;  The both components correspond to delta, in EMF=...+delta x J = delta*(-Jy,+Jx,0)
;
!x.title='!8t!6'
yrdel=[-.5,1.]*delmax
plot,tt,eta12,yr=yrdel,ytit='!7-d!6'
if itestfield eq 'B11-B22' then begin
  oplot,tt,-eta21,li=2
  del=.5*(eta12-eta21)
endif else begin
  del=.5*eta12
endelse
oplot,tt,urms/3.,li=1
pc_error_range,tt(good),del(good),mean=delm,error=del_error
pc_error_range,tt(good),+eta12(good),/oplot,mean=eta12m,error=eta12_error
pc_error_range,tt(good),-eta21(good),/oplot,mean=eta21m,error=eta21_error
oplot,tt(good),+accum(eta12(good)),col=55
oplot,tt(good),-accum(eta21(good)),col=55
oplot,tt,tt-tt,li=3
;
;  give modified etamax values in parameters.pro file
;
yreta=[-.2,1.]*etamax
plot,tt,eta22,yr=yreta
if itestfield eq 'B11-B22' then begin
  oplot,tt,eta11,li=2
  eta=.5*(eta11+eta22)
endif else begin
  eta=eta22
endelse
oplot,tt,urms/(3.*kf),li=1
pc_error_range,tt(good),eta(good),/oplot,mean=etam,error=eta_error
oplot,tt(good),accum(eta(good)),col=55
oplot,tt,tt-tt,li=3
;oplot,tt,etat,col=122,li=0
;
pc_error_range,tt(good),urms(good),/oplot,mean=urmsm,error=urmsm_error
;
;  print on screen
;
print
fo='(2(a,f6.4))'
fo2='(2(a,f8.5))'
fo2='(2(a,e9.2))'
print,' eta_t=',etam,' +/-',eta_error,fo=fo2
print,' del_t=',delm,' +/-',del_error,fo=fo2
print,'+eta12=',eta12m,' +/-',eta12_error,fo=fo2
print,'-eta21=',eta21m,' +/-',eta21_error,fo=fo2
print,'  urms=',urmsm,' +/-',urmsm_error,fo=fo
print
!p.multi=0
;
fo='(2e8.1,i6,8e10.2,"  ",a)'
pc_read_param,/param2,obj=param
if iuu eq 0 then nu=0. else nu=param.nu
print,nu,param.etatest,fix(max(tt)),etam,eta_error,eta12m,eta12_error,eta21m,eta21_error,urmsm,urmsm_error,run,fo=fo
print
END
