;;;;;;;;;;;;;;;;;;;;;;
;;;   thermo.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   25-Feb-2002
;;;
;;;  Description:
;;;   Calculate all relevant thermodynamical variables from lnrho and
;;;   ss.

pp  = cs20*rho0/gamma*exp(gamma*(ss+lnrho-lnrho0))
cs2 = cs20 * exp(gamma*ss+gamma1*(lnrho-lnrho0))
TT  = cs2/gamma1

;; initial profiles of temperature, density and entropy for solar
;; convection runs
ssinit = (lnrhoinit = (Tinit = 0.*z))

top = where(z ge z2)
Tref = cs20/gamma1
lnrhoref = alog(rho0)
ssref = 0.
zo = [z2,ztop]
zint = zref                     ; position of top/unstable layer interface

if (top[0] ge 0) then begin
  zint = z2 < ztop
  if (isothtop eq 0) then begin ; polytropic top layer
    beta = gamma/gamma1*gravz/(mpoly2+1)
    Tinit[top] = Tref + beta*(z[top]-zref)
    lnrhoinit[top] = lnrhoref + mpoly2*alog(1.+beta*(z[top]-zref)/Tref)
    ssinit[top] = ssref $
                  + (1-mpoly2*(gamma-1))/gamma $
                    * alog(1.+beta*(z[top]-zref)/Tref)
;;;    zint = zrefz2 < ztop
    lnrhoref = lnrhoref + mpoly2*alog(1.+beta*(zint-zref)/Tref)
    ssref = ssref $
            + (1-mpoly2*(gamma-1))/gamma * alog(1.+beta*(zint-zref)/Tref)
    Tref = Tref + beta*(zint-zref)
  endif else begin              ; isothermal top layer
    beta = 0.
    Tinit[top] = Tref
    lnrhoinit[top] = lnrhoref + gamma/gamma1*gravz*(z[top]-zref)/Tref
    ssinit[top] = ssref - gravz*(z[top]-zref)/Tref
    lnrhoref = lnrhoref + gamma/gamma1*gravz*(z2-zref)/Tref
    ssref = ssref  - gravz*(zint-zref)/Tref
    Tref = Tref
  endelse
endif

;
stab = where((z le z2) and (z ge z1))
if (stab[0] ge 0) then begin
  beta = gamma/gamma1*gravz/(mpoly0+1)
  Tinit[stab] = Tref + beta*(z[stab]-zint)
  lnrhoinit[stab] = lnrhoref + mpoly0*alog(1.+beta*(z[stab]-zint)/Tref)
  ssinit[stab] = ssref $
                 + (1-mpoly0*(gamma-1))/gamma $
                   * alog(1.+beta*(z[stab]-zint)/Tref)
endif
;
unstab = where(z le z1)
if (unstab[0] ge 0) then begin
  lnrhoref = lnrhoref + mpoly0*alog(1.+beta*(z1-z2)/Tref)
  ssref = ssref $
          + (1-mpoly0*(gamma-1))/gamma $
            * alog(1.+beta*(z1-z2)/Tref)
  Tref = Tref + beta*(z1-z2)
  beta = gamma/gamma1*gravz/(mpoly1+1)
  Tinit[unstab] = Tref + beta*(z[unstab]-z1)
  lnrhoinit[unstab] = lnrhoref + mpoly1*alog(1.+beta*(z[unstab]-z1)/Tref)
  ssinit[unstab] = ssref $
                   + (1-mpoly1*(gamma-1))/gamma $
                     * alog(1.+beta*(z[unstab]-z1)/Tref)
endif

; anything else?


end
; End of file thermo.pro
