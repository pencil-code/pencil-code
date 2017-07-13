;$Id: pspiegel.pro,v 1.7 2007/05/12 07:43:28 brandenb Exp $
;
n=1000
;
sigmaSB=5.67d-5 ;(erg/cm2/s/K)
G=6.67d-8       ;(cm3/g/s2)
c=3d10          ;(cm/s)
kappa_es=3.4d-1
Rgas=8.31434d7  ;erg/mol/K
;
Lsun=3.9e33 ;erg/s
Msun=1.99e33 ;(g)
Rsun=7e10 ;(cm)
;
; http://en.wikipedia.org/wiki/Main_sequence
;
;  F0 star
;
L=Lsun*9.0
R=Rsun*1.5
M=Msun*1.6
;
;  scaled solar like star
;
L=Lsun*2d4
R=Rsun
M=Msun
;
;  B0 star
;
L=Lsun*16d3
R=Rsun*5.7
M=Msun*16.
;
grav=G*M/R^2
Flux=L/(4.*!pi*R^2)
Teff=(Flux/sigmaSB)^.25
grav_rad=kappa_es*flux/c
grav_eff=grav-grav_rad
arad=4.*sigmaSB/c
;
;  define Theta range
;
That=0.9*Teff
;That=1e5
Theta=exp(grange(1e-4,1.3,n))
T=Theta*That
;
;  Theta - .5*atan(Theta) - .5*acoth(Theta)
;
rhs=Theta-.5*atan(Theta)-.5*atanh(1./Theta)
z=-rhs*4.*2.*Rgas*That/grav_eff
;
Mm=1d8
plot,z/Mm,T,ps=-1
;
E=arad*T^4
E0=arad*That^4
p=c*grav_eff/(3.*kappa_es*flux)*(E-E0)
plot,z/Mm,p,ps=-1
;
;  factor 2 from mu=1/2 for full ionization
;
rho=p/(2.*Rgas*T)
pgas=2.*Rgas*T/rho
plot_io,z/Mm,rho,ps=-1
;
nz=256 & mz=nz+6
z1=-15d0 & z2=15d0
z1=-90d0 & z2=30d0
z1=-90d0 & z2=-30d0
z1=-30d0 & z2=30d0
default,z1,-30d0
default,z2,+30d0
@parameters
dz=(z2-z1)/(nz-1.)
;
zi=grange(z1-3.*dz,z2+3.*dz,mz)
Ti=interpol(T,z/Mm,zi)
rhoi=interpol(rho,z/Mm,zi)
oplot,zi,rhoi,ps=-5,col=122
;
;------------ blob-like perturbation ------------------
w=2
zi0=-15.
zi0=+15.
ampl=.02
ampl=.0
blob=exp(-(zi-zi0)^2/w^2)
ptot=2.*Rgas*rhoi*Ti+arad*Ti^4/3.
;
;  perturb
;
Ti=Ti*(1.+ampl*blob)
rhoi=(ptot-arad*Ti^4/3.)/(2.*Rgas*Ti)
;
;--------- end of blob-like perturbation --------------
;
;
lnTTi=alog(Ti)
lnrhoi=alog(rhoi/1e-6)
;
;  calculate mean temperature and mean density
;
lnTTm=alog(mean(exp(lnTTi(3:nz+2))))
lnrhom=alog(mean(exp(lnrhoi(3:nz+2))))
;
fo='(3f13.3)'
fo='(3f13.6)'
close,1
openw,1,'data/proc0/stratification.ascii'
for i=0,mz-1 do begin
  print,zi(i),lnrhoi(i),lnTTi(i),fo=fo
  printf,1,zi(i),lnrhoi(i),lnTTi(i),fo=fo
  ;printf,1,zi(i),lnrhom,lnTTm,fo=fo
endfor
close,1
plot,zi,lnTTi,ps=-1
;
;  print average values
;
length_unit=1e8
density_unit=1e-6
velocity_unit=1e5
flux_unit=density_unit*velocity_unit^3
grav_unit=velocity_unit^2/length_unit
print,'flux,flux/flux_unit=',flux,flux/flux_unit
print,'grav,grav/grav_unit=',grav,grav/grav_unit
print,'grav_rad,grav_rad/grav_unit=',grav_rad,grav_rad/grav_unit
;
;  pressure
;
lnppi=alog(2.*Rgas)+lnrhoi+lnTTi
;
END
