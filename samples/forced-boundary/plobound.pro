;;;;;;;;;;;;;;;;;;;;;;;;
;;;   plobound.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   29-Apr-2004
;;;
;;;  Description:
;;;   Plot magnetic field and velocity at lower boundary

save_state

!p.multi=[0,2,1]
!x.title='!8x!X'
!y.title='!8y!X'
bb = curl(aa)

iz1 = n1                        ; vertical level to plot

plot_3d_vect, uu[l1:l2,m1:m2,iz1,0], uu[l1:l2,m1:m2,iz1,1], $
              bb[l1:l2,m1:m2,iz1,2], $
    x[l1:l2], y[l1:l2], /KEEP, $
    POS=aspect_pos(Ly/Lx), XSTYLE=1, YSTYLE=1, $
    TITLE='!8B!Dz!N!6 and !17u!X'

plot_3d_vect, bb[l1:l2,m1:m2,iz1,*], $
    x[l1:l2], y[l1:l2], $
    POS=aspect_pos(Ly/Lx), XSTYLE=1, YSTYLE=1, $
    TITLE='!17B!X'

restore_state



end
; End of file plobound.pro
