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

;; initial temperature profile for solar convection runs
Tinit = 0.*z
top = where(z ge z2)
if (top[0] ge 0) then begin
  Tref = cs0^2/gamma1
  zo = [z2,ztop]
  beta = gamma/gamma1*gravz/(mpoly2+1)
  Tinit[top] = Tref + beta*(z[top]-ztop)
endif
;
stab = where((z le z2) and (z ge z1))
if (stab[0] ge 0) then begin
  Tref = Tref + beta*(z2-ztop)
  beta = gamma/gamma1*gravz/(mpoly0+1)
  Tinit[stab] = Tref + beta*(z[stab]-z2)
endif
;
unstab = where(z le z1)
if (unstab[0] ge 0) then begin
  Tref = Tref + beta*(z1-z2)
  beta = gamma/gamma1*gravz/(mpoly1+1)
  Tinit[unstab] = Tref + beta*(z[unstab]-z1)
endif

; anything else?


end
; End of file thermo.pro
