;--------------
; Main program
;--------------

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
m_He=3.97153*m_H
fHe=.1
mu=1.+3.97153*fHe
chiH=13.6*eV
lnrho_ion=1.5*alog((m_e/hbar)*(chiH/hbar)/(2.*!pi))+alog(m_H)+alog(mu)
ss_ion=k_B/m_H/mu
TT_ion=chiH/k_B
twothirds=2./3.

xxx=x(l1:l2)
yyy=y(m1:m2)
zzz=z(n1:n2)
uuu=uu(l1:l2,m1:m2,n1:n2,*)
sss=ss(l1:l2,m1:m2,n1:n2)
llnrho=lnrho(l1:l2,m1:m2,n1:n2)
rho=exp(llnrho)
yH=yyH(l1:l2,m1:m2,n1:n2)
TTT=TT(l1:l2,m1:m2,n1:n2)

lnTT_=alog(TTT/TT_ion)
f=lnrho_ion-llnrho+1.5*lnTT_-exp(-lnTT_)+alog(1.-yH)-2.*alog(yH)
dlnTT_dy=(alog(m_H/m_p)-gamma1*(f+exp(-lnTT_))-1.)/(1.+yH+fHe)
dfdy=dlnTT_dy*(1.5+exp(-lnTT_))-1./(1.-yH)-2./yH

dlnTT_dlnrho=gamma1
dfdlnrho=gamma1*exp(-lnTT_)
dydlnrho=-dfdlnrho/dfdy
dlnPdlnrho=1.+dydlnrho/(1.+yH+fHe)+dlnTT_dy*dydlnrho+dlnTT_dlnrho

cs2=ss_ion*TTT*dlnPdlnrho

PPP=(1.+yH+fHe)*rho*ss_ion*TTT

END
