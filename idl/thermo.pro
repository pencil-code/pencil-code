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

; anything else?


end
; End of file thermo.pro
