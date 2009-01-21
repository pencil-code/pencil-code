; $Id$
;
; Plot entropy and velocity vectors for 2d convection runs
;

save_state

;!p.charsize=2
!p.charsize=1.4
;!p.charthick=2
;!p.thick=2
;!x.thick=2
;!y.thick=2

bgcol = !d.table_size-1
fgcol = 0
pos = aspect_pos(Lz/Lx, MARGIN=0.08)
contourfill, ss[l1:l2,3,n1:n2], x[l1:l2], z[n1:n2], $
    BACK=bgcol, COLOR=fgcol, $
    POSITION=pos, XSTYLE=1, YSTYLE=1, $
    XTITLE='!8x!X', YTITLE='!8z!X'

;velovect,reform(uu(l1:l2,3,n1:n2,0)),reform(uu(l1:l2,3,n1:n2,2)),x(l1:l2),z(n1:n2),len=2,/over
;!p.charthick=1 & !p.thick=1 & !x.thick=1 & !y.thick=1
vel_a, reform(uu[l1:l2,3,n1:n2,0]), reform(uu[l1:l2,3,n1:n2,2]), $
    x[l1:l2],z[n1:n2], $
    /OVER, NVEC=1800, COLOR=bgcol

restore_state

; print,'import Hurlburt84-Ra1e5.png'

end
