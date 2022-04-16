kms=1e5
u=30.*kms
sigmaSB=5.67d-5	;(erg/cm2/s/K4)
Rgas=8.31434d7  ;erg/mol/K
mu=.6
gam=5./3.
nabad=1.-1./gam
cp=Rgas/mu/nabad
print,'cp=',cp
print,'sigmaSB=',sigmaSB
;
rho=sigmaSB*u^5/cp^4
print,'[rho]=',rho
;
rho=4e-4
sigmaSB_=cp^4*rho/u^5
;
print,'sigmaSB_/1e5=',sigmaSB_/1e5
print,'sigmaSB_/sigmaSB=',sigmaSB_/sigmaSB
;
yr=3e7
g=274e2
tauKH=rho*cp^4/(sigmaSB*g*u^4)
print,'tauKH [yr]=',tauKH/yr
;
END
