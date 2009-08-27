;;;;;;;;;;;;;;;;;;;;;;
;;;   strati.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   23-Nov-2001
;;;  $Id$
;;;
;;;  Description:
;;;   Plot vertical profile of p, rho, s, etc for a stratified
;;;   atmosphere with linear entropy profile.

Nz = 100
gamma = 5.D/3.
cs0 = 1.
rho0 = 1.
s0 = -alog(gamma)/gamma
dsdz0 = -0.2D0
gz = -0.1D0

gam1 = gamma-1
cs20 = cs0^2
z = linspace(-3.6D0,3.6,Nz)
ent = s0 + dsdz0*z
lnrho = alog(rho0) $
        - dsdz0*z $
        + 1./gam1*alog(1 + gam1*gz/dsdz0/cs20*(1.-exp(-dsdz0*z)))
pp = exp(gamma*(ent + lnrho))
cs2 = gamma*pp*exp(-lnrho)

; ----

save_state

!p.multi = [0,3,2]
!p.charsize = 2
!y.title = '!8z!X'

plot, ent, z, XTITLE='!8s!X', TITLE='!6Entropy!X'
;
plot, lnrho, z, XTITLE='!6ln !7r!X', TITLE='!6Density!X'
;
plot, cs2, z, XTITLE='!8c!6!Ds!N!U2!N!X', TITLE='!6Sound speed!X'
;
;
plot, pp, z, XTITLE='!8p!X', TITLE='!6Pressure!X'
;
var = -deriv(z,pp)*exp(-lnrho) ; + gz
plot, var, z, XTITLE='!8-(dp/dz) / !7r!X', TITLE='!6Pressure term!X', $
    XRANGE=minmax([var,0.,-gz]), XSTYLE=3
opvline, -gz, LINE=2, COLOR=120
opvline
;
var = var+gz
var2 = -cs2*(deriv(z,ent) + deriv(z,lnrho)) + gz ; different form of press. term
plot, var, z, $
    TITLE='!6Error in pressure term!X', $
    XRANGE=[-1.1,1.2]*max(abs([var,var2])), XSTYLE=3
oplot, var2, z, LINE=2
opvline

restore_state

; ----

end
; End of file strati.pro
