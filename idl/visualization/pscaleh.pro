;;;;;;;;;;;;;;;;;;;;;;;
;;;   pscaleh.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   25-Feb-2002
;;;  $Id$
;;;
;;;  Description:
;;;   Plot pressure and pressure scale height as function of z.
;;;   H_p is calculated in two ways:
;;;     H_p = 1/(d ln p/dz)
;;;   and
;;;     H_p = cs^2/(gamma |g_z|)
;;;

save_state

!p.multi=[0,2,1]
!y.title = '!8z!X'

pp = exp(gamma*(ent+lam))
cs2 = gamma * exp(gamma*ent+gamma_m1*lam)

plot_binned, alog(pp), zz, $
    PSYM=1,  TITLE='!6Pressure!X', XTITLE='!6ln !8p!X'
ophline,[-0.68,0.,1,z[nz-4]]

plot_binned, -1./zder(alog(pp)), zz, $
    PSYM=1, TITLE='!6Pressure scale height!X', $
    XTITLE='!8dz/d!6!E !Nln!E !N!8p!X'
plot_binned, /OVER, cs2/gamma/abs(gravz), zz, $
    PSYM=4
ophline,[z0,z1,z2,ztop]

opvline,0.45 & ophline,0.5

restore_state

end
; End of file pscaleh.pro
