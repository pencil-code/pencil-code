;--------------
; Main program
;--------------

COMMON constants,m_e,m_p,m_H,chiH,lnrho_ion,ss_ion,TT_ion,twothirds

hbar_cgs=1.0545726663d-27 ; [erg*s]
k_B_cgs=1.38065812d-16    ; [erg/K]
m_p_cgs=1.672623110d-24   ; [g]
m_e_cgs=9.109389754d-28   ; [g]
eV_cgs=1.602177250d-12    ; [erg]


unit_mass=unit_density*unit_length^3
unit_energy=unit_mass*unit_velocity^2
unit_time=unit_length/unit_velocity
unit_temperature=1.

IF (unit_system='cgs') THEN BEGIN
  hbar=hbar_cgs/(unit_energy*unit_time)
  k_B=k_B_cgs/(unit_energy/unit_temperature)
  m_p=m_p_cgs/unit_mass
  m_e=m_e_cgs/unit_mass
  eV=eV_cgs/unit_energy
  ENDIF ELSE IF (unit_system='SI') THEN BEGIN
  print,'units of length,velocity,density are given in SI'
  k_B=1e-7*k_B_cgs/(unit_energy/unit_temperature)
  m_p=m_p_cgs*1e-3/unit_mass
  m_e=m_e_cgs*1e-3/unit_mass
  eV=eV_cgs*1e-7/unit_energy
ENDIF

m_H=m_p+m_e
chiH=13.6*eV
lnrho_ion=1.5*alog((m_e/hbar)*(chiH/hbar)/(2.*!pi))+alog(m_H)
ss_ion=k_B/m_H
TT_ion=chiH/k_B
twothirds=2./3.

xxx=x(l1:l2)
yyy=y(m1:m2)
zzz=z(n1:n2)
uuu=uu(l1:l2,m1:m2,n1:n2,*)
sss=ss(l1:l2,m1:m2,n1:n2)
llnrho=lnrho(l1:l2,m1:m2,n1:n2)
rho=exp(llnrho)
yyH=yH(l1:l2,m1:m2,n1:n2)

lnTT_=twothirds*((sss/ss_ion-1.5*(1.-yyH)*alog(m_H/m_e) $
                  -1.5*yyH*alog(m_p/m_e)+(1.-yyH)*alog(1.-yyH) $
                  +2.*yyH*alog(yyH))/(1.+yyH)+llnrho-lnrho_ion-2.5)
f=lnrho_ion-lnrho+1.5*lnTT_-exp(-lnTT_)+alog(1.-yyH)-2.*alog(yyH)
dlnTT_dy=(alog(m_H/m_p)-gamma1*(f+exp(-lnTT_))-1.)/(1.+yyH)
dfdy=dlnTT_dy*(1.5+exp(-lnTT_))-1./(1.-yyH)-2./yyH
dlnTT_dlnrho=gamma1
dfdlnrho=gamma1*exp(-lnTT_)
dydlnrho=-dfdlnrho/dfdy
dlnPdlnrho=1.+dydlnrho/(1.+yyH)+dlnTT_dy*dydlnrho+dlnTT_dlnrho
dlnTT_dss=gamma1/((1.+yyH)*ss_ion)
dfdss=(1.+dfdlnrho)/((1.+yyH)*ss_ion)
dydss=-dfdss/dfdy
dlnPdss=dydss/(1.+yyH)+dlnTT_dy*dydss+dlnTT_dss
;kappa=.25*exp(lnrho-lnrho_ion)*(TT_ion_/TT)^1.5 $
;*exp(TT_ion_/TT)*yyH*(1.-yyH)*kappa0

TTT=exp(lnTT_)*TT_ion
PPP=(1.+yyH)*rho*ss_ion*TTT
cs2=(1.+yyH)*ss_ion*TT*dlnPdlnrho
cp1tilde=dlnPdss/dlnPdlnrho
END
