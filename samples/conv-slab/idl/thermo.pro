;;;;;;;;;;;;;;;;;;;;;;
;;;   thermo.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   25-Feb-2002
;;;  $Id$
;;;
;;;  Description:
;;;   Calculate all relevant thermodynamical variables from lnrho and
;;;   ss for box convection runs.
;;;  Usage:
;;;   cd ${PENCIL_HOME}/samples/conv-slab; make; start.csh; run.csh; idl
;;;   .r start
;;;   .r rall
;;;   .r thermo
;;;   .r pvert

;; Get cp (make sure par.cp exists and default to 1)
default, STRUCT=par, 'cp', 1.
cp = par.cp

pp  = cs20*rho0/gamma*exp(gamma*(ss/cp + lnrho-lnrho0))
cs2 = cs20*exp(gamma*ss/cp + gamma_m1*(lnrho-lnrho0))
TT  = cs2/gamma_m1/cp

z1 = par.z1
z2 = par.z2
ztop = z0 + Lz
zref = par.zref
gravz = par.gravz

;; initial profiles of temperature, density and entropy for solar
;; convection runs
ssinit = (lnrhoinit = (Tinit = 0.*z))

top = where(z ge z2)
Tref = cs20/gamma_m1/cp
lnrhoref = alog(rho0)
ssref = 0.
zo = [z2,ztop]
zint = zref                     ; position of top/unstable layer interface

if (top[0] ge 0) then begin
  zint = z2 < ztop
  if (isothtop eq 0) then begin ; polytropic top layer
    beta = gamma/gamma_m1/cp*gravz/(mpoly2+1)
    Tinit[top] = (Tref + beta*(z[top]-zref)) > 1.e-5*Tref ; (or rather 1.e-20?)
    lnrhoinit[top] = lnrhoref + mpoly2*alog(Tinit[top]/Tref)
    ssinit[top] = ssref $
                  + (1-mpoly2*(gamma-1))/gamma $
                    * alog(Tinit[top]/Tref)
;;;    zint = zrefz2 < ztop
    lnrhoref = lnrhoref + mpoly2*alog(1.+beta*(zint-zref)/Tref)
    ssref = ssref $
            + (1-mpoly2*(gamma-1))/gamma * alog(1.+beta*(zint-zref)/Tref)
    Tref = Tref + beta*(zint-zref)
  endif else begin              ; isothermal top layer
    beta = 0.
    Tinit[top] = Tref
    lnrhoinit[top] = lnrhoref + gamma/gamma_m1*gravz*(z[top]-zref)/Tref
    ssinit[top] = ssref - gravz*(z[top]-zref)/Tref
    lnrhoref = lnrhoref + gamma/gamma_m1*gravz*(z2-zref)/Tref
    ssref = ssref  - gravz*(zint-zref)/Tref
    Tref = Tref
  endelse
endif

;
stab = where((z le z2) and (z ge z1))
if (stab[0] ge 0) then begin
  beta = gamma/gamma_m1/cp*gravz/(mpoly0+1)
  Tinit[stab] = Tref + beta*(z[stab]-zint)
  lnrhoinit[stab] = lnrhoref + mpoly0*alog(Tinit[stab]/Tref)
  ssinit[stab] = ssref $
                 + (1-mpoly0*(gamma-1))/gamma $
                   * alog(Tinit[stab]/Tref)
endif
;
unstab = where(z le z1)
if (unstab[0] ge 0) then begin
  lnrhoref = lnrhoref + mpoly0*alog(1.+beta*(z1-z2)/Tref)
  ssref = ssref $
          + (1-mpoly0*(gamma-1))/gamma $
            * alog(1.+beta*(z1-z2)/Tref)
  Tref = Tref + beta*(z1-z2)
  beta = gamma/gamma_m1/cp*gravz/(mpoly1+1)
  Tinit[unstab] = Tref + beta*(z[unstab]-z1)
  lnrhoinit[unstab] = lnrhoref + mpoly1*alog(Tinit[unstab]/Tref)
  ssinit[unstab] = ssref $
                   + (1-mpoly1*(gamma-1))/gamma $
                     * alog(Tinit[unstab]/Tref)
endif


end
; End of file thermo.pro
