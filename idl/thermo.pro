;;;;;;;;;;;;;;;;;;;;;;
;;;   thermo.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   25-Feb-2002
;;;
;;;  Description:
;;;   Calculate all relevant thermodynamical variables from lam and
;;;   ent.

pp  = exp(gamma*(ent+lam))
cs2 = gamma * exp(gamma*ent+gamma1*lam)
TT  = cs2/gamma1

;; initial profiles of temperature, density and entropy for solar
;; convection runs
entinit = (laminit = (Tinit = 0.*z))

top = where(z ge z2)
if (top[0] ge 0) then begin
  Tref = cs0^2/gamma1
  lamref = alog(rho0)
  entref = (2*alog(cs0) - gamma1*alog(rho0)-alog(gamma))/gamma
  zo = [z2,ztop]
  beta = gamma/gamma1*gravz/(mpoly2+1)
  Tinit[top] = Tref + beta*(z[top]-ztop)
  laminit[top] = lamref + mpoly2*alog(1.+beta*(z[top]-ztop)/Tref)
  entinit[top] = entref $
                 + (1-mpoly2*(gamma-1))/gamma $
                   * alog(1.+beta*(z[top]-ztop)/Tref)
endif
;
stab = where((z le z2) and (z ge z1))
if (stab[0] ge 0) then begin
  lamref = lamref + mpoly2*alog(1.+beta*(z2-ztop)/Tref)
  entref = entref $
           + (1-mpoly2*(gamma-1))/gamma * alog(1.+beta*(z2-ztop)/Tref)
  Tref = Tref + beta*(z2-ztop)
  beta = gamma/gamma1*gravz/(mpoly0+1)
  Tinit[stab] = Tref + beta*(z[stab]-z2)
  laminit[stab] = lamref + mpoly0*alog(1.+beta*(z[stab]-z2)/Tref)
  entinit[stab] = entref $
                  + (1-mpoly0*(gamma-1))/gamma $
                    * alog(1.+beta*(z[stab]-z2)/Tref)
endif
;
unstab = where(z le z1)
if (unstab[0] ge 0) then begin
  lamref = lamref + mpoly0*alog(1.+beta*(z1-z2)/Tref)
  entref = entref $
                 + (1-mpoly0*(gamma-1))/gamma $
                   * alog(1.+beta*(z1-z2)/Tref)
  Tref = Tref + beta*(z1-z2)
  beta = gamma/gamma1*gravz/(mpoly1+1)
  Tinit[unstab] = Tref + beta*(z[unstab]-z1)
  laminit[unstab] = lamref + mpoly1*alog(1.+beta*(z[unstab]-z1)/Tref)
  entinit[unstab] = entref $
                    + (1-mpoly1*(gamma-1))/gamma $
                      * alog(1.+beta*(z[unstab]-z1)/Tref)
endif

; anything else?


end
; End of file thermo.pro
