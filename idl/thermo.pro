;;;;;;;;;;;;;;;;;;;;;;
;;;   thermo.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   25-Feb-2002
;;;  $Id: thermo.pro,v 1.10 2006-03-27 17:05:47 ngrs Exp $
;;;
;;;  Description:
;;;   Calculate all relevant thermodynamical variables from lnrho and
;;;   ss.
;;;  
;;;  Modified 24-3-2006, for use with cp /= 0 runs, ngrs & neamvonk

cp = par.cp
mpoly = par.mpoly
beta = gamma/gamma1/cp*gravz/(mpoly+1)       ; gamma1*cp/gamma=R_{*}
TT0 = beta*z2
lnTT0=alog(TT0)

rho=exp(lnrho)
pp  = rho/gamma * cs20*exp(gamma*ss/cp + gamma1*(lnrho-lnrho0))       ; use cs2=rho pp/gamma
cs2 = cs20*exp(gamma*ss/cp + gamma1*(lnrho-lnrho0))
TT  = cs2/gamma1/cp
;;TT = exp(lnTT0 + gamma*ss/cp + gamma1*(lnrho-lnrho0))  ;; alternative formula, also correct for general cp

;; initial profiles of temperature, density and entropy for solar
;; convection runs
ssinit = (lnrhoinit = (Tinit = 0.*z))

top = where(z ge z2)
Tref = cs20/gamma1/cp
lnrhoref = alog(rho0)
ssref = 0.
zo = [z2,ztop]
zint = zref                     ; position of top/unstable layer interface

if (top[0] ge 0) then begin
  zint = z2 < ztop
  if (isothtop eq 0) then begin ; polytropic top layer
    Tinit[top] = Tref + beta*(z[top]-zref)
    lnrhoinit[top] = lnrhoref + mpoly2*alog( (1.+beta*(z[top]-zref)/Tref) > 1.e-5 )
    ssinit[top] = ssref $
                  + (1-mpoly2*(gamma-1))/gamma $
                    * alog( (1.+beta*(z[top]-zref)/Tref) > 1.e-5 )
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
  ggamma = 1. + (1./mpoly0)
  zinf = zref + (mpoly0+1.)*cs20/gamma/(-gravz)
  phi = (z[stab] - zinf)*(-gravz)
  Tinit[stab] = Tref + beta*(z[stab]-zint)
  lnrhoinit[stab] = lnrhoref + (mpoly0 * alog(-gamma*phi/(mpoly0+1.)/cs20))
  ssinit[stab] = ssref + ((ggamma/gamma -1) * mpoly0 * alog(-gamma*phi/(mpoly0+1.)/cs20))*cp

endif
;
unstab = where(z le z1)
if (unstab[0] ge 0) then begin
  lnrhoref = lnrhoref + mpoly0*alog(1.+beta*(z1-z2)/Tref)
  ssref = ssref $
          + (1-mpoly0*(gamma-1))/gamma $
            * alog(1.+beta*(z1-z2)/Tref)
  Tref = Tref + beta*(z1-z2)
  Tinit[unstab] = Tref + beta*(z[unstab]-z1)
  lnrhoinit[unstab] = lnrhoref + mpoly1*alog(1.+beta*(z[unstab]-z1)/Tref)
  phi = (z[unstab] - zinf)*(-gravz)
  ssinit[unstab] = ((ggamma/gamma -1) * mpoly0 * alog(-gamma*phi/(mpoly0+1)/cs20))*cp
  
endif

; anything else?

end
; End of file thermo.pro
