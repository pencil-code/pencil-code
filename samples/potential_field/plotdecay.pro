;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   plotdecay.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ucalgary.ca)
;;;  Date:   11-Oct-2006
;;;
;;;  Description:
;;;    Plot decay of magnetic field etc.
;;;  Usage:
;;;    @plotdecay


pc_read_param, /PARAM2, OBJECT=par2, /QUIET
eta = par2.eta

pc_read_ts, OBJECT=ts

; kx =  1. & mu_tilde =  3.6731944
; kx =  1. & mu_tilde =  9.6316846
; kx =  1. & mu_tilde = 15.834105
; kx =  1. & mu_tilde = 22.081660
; kx =  1. & mu_tilde = 28.344864
;
; kx =  2. & mu_tilde =  2.02875785
; kx =  4. & mu_tilde =  1.14446488
; kx =  8. & mu_tilde =  0.64260789
  kx = 12. & mu_tilde =  0.452743291

gamma = eta*kx^2*(1 + mu_tilde^2)

N = n_elements(ts.t)
timefact = exp(-gamma*(ts.t-ts.t[N-1]))

save_state

!p.multi = [0, 2, 2]
!x.title = '!8t!X'
!x.range = [0, 1]
!y.range = 1 + [-1,1]*0.3
xr = minmax(ts.t)
idx=[0,1,2,3,4]

var = ts.brms/mean(ts.brms[idx]/timefact)
plot, ts.t, var/timefact, PSYM=-4, TITLE='brms'
oplot, xr, [1,1], color=120

var = ts.jrms/mean(ts.jrms[idx]/timefact)
plot, ts.t, var/timefact, PSYM=-4, TITLE='jrms'
oplot, xr, [1,1], color=120

var = ts.bmax/mean(ts.bmax[idx]/timefact)
plot, ts.t, var/timefact, PSYM=-4, TITLE='bmax'
oplot, xr, [1,1], color=120

var = ts.jmax/mean(ts.jmax[idx]/timefact)
plot, ts.t, var/timefact, PSYM=-4, TITLE='jmax'
oplot, xr, [1,1], color=120

restore_state
; End of file plotdecay.pro
