; $Id: $
;
; 07-01-2011 PJK
;
; Generates piecewise polytropic stratification and heat conductivity 
; profiles for spherical convection runs. The upper layer between 
; r2/R and r/R=1. is isothermal.
;
; nr         : Number of grid points in radius. Must match nxgrid in 
;              src/cparam.local
; Tbot       : Temperature at r=rbot
; rho_bot    : Density at r=rbot
; rbot       : Radial position of the bottom of the domain in units of 
;              the stellar radius R, must match xyz0[1] in start.in
; r1         : Base of the convectively unstable region in units of the
;              stellar radius R
; r2         : Top of the convectively unstable region in units of the 
;              stellar radius R
; GM         : Gravity constant times mass of the star, must match 
;              gravx in start.in/run.in
; d          : Depth of transition layers at stable/unstable interfaces
; Luminosity : Luminosity of the star
;
; The outer boundary is assumed to be at r/R=1.
;
nr=32 & Tbot=2.418 & rho_bot=10. & rbot=0.6 & r1=0.70 & r2=0.97 & GM=3. & d=0.015 & Luminosity=0.0404275
;
r=rbot+findgen(nr)/(nr-1)*(1.-rbot) & dr=r(1)-r(0) ; Radius
;
gamma=5./3. & cp=1. & cv=cp/gamma ; Heat capacities
;
gr=-GM/r^2 ; Gravity proportional to r^(-2)
;
rho=fltarr(nr) & T=fltarr(nr) & lnr=fltarr(nr) & m=fltarr(nr)
;
; Define a depth dependent polytropic index that describes the
; stratification. For the corresponding equations for dT/dr and
; drho/dr see Käpylä et al. (2004), Astron. Astrophys., 422, 793-816
;
m=-1.*tanh((r-r1)/d)+2.
;
; Temperature stratification
;
dTdr =gr/(cv*(gamma-1.)*(m+1.))     ; T-gradient for thermal stratification
dTdrc=gr/(cv*(gamma-1.)*(m+1.))     ; T-gradient for heat conductivity
good=where(r ge r2) & dTdr(good)=0. ; Isothermal region
;
T(0)=Tbot
for i=1,nr-1 do begin
  T(i)=T(i-1)+dTdr(i-1)*dr
endfor
;
cs2=cV*T*(gamma-1.)*gamma ; Sound speed squared
;
dlnTdr=deriv(r,T)/T
dlnrdr=gr/(T*cv*(gamma-1.))-dlnTdr
;
; Density stratification
;
lnr(0)=alog(rho_bot)
for i=1,nr-1 do begin
  lnr(i)=lnr(i-1)+dlnrdr(i-1)*dr
endfor
;
; Thermal conductivity
;
kappa=-Luminosity/(4.*!pi*r^2*dTdrc) ; Heat conductivity
dlogkappa=deriv(r,kappa)             ; Gradient of heat conductivity
Fbot=-kappa(0)*dTdr(0)               ; Energy flux at the base
Ftop=-kappa(nr-1)*dTdrc(nr-1)        ; Energy flux at the surface
;
rho=exp(lnr)                         ; Density
p=rho*cv*T*(gamma-1.)                ; Pressure
s=alog(p)/gamma-lnr                  ; Entropy
;
chi=kappa/(cp*rho)                   ; Thermal diffusivity
;
; Brunt-Väisalä frequency
;
BV=r*(deriv(r,p)/(gamma*p)-deriv(r,rho)/rho)
;
plot,r,rho/rho_bot,yr=[-0.2,1.05],ys=3,xtit='!8r!6/!8R!6',tit='!7q!8, T, p, !6B-V freq.'
oplot,r,T/Tbot,li=2
oplot,r,p/p(0),li=3
oplot,r,r*0,li=1
oplot,r,BV,li=4
;
print,''
print,'Writing stratification.dat...'
openw,1,'stratification.dat'
for i=0,nr-1 do printf,1,alog(rho(i)),alog(rho(i)),alog(T(i))
close,1
;
print,'Writing hcond_glhc.dat...'
openw,1,'hcond_glhc.dat'
for i=0,nr-1 do printf,1,kappa(i),dlogkappa(i)
close,1
print,''
print,'Place the files stratification.dat and hcond_glhc.dat in the run'
print,'directory of your simulation and insert the following lines in '
print,'run.in:'
print,''
print,'&density_run_pars:'
print,'  cs2top=',cs2(nr-1)
print,''
print,'&entropy_run_pars:'
print,'  Fbot=',Fbot
print,'  cs2cool=',cs2(nr-1)
print,''
;
end
