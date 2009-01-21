; $Id$
;
; scp Hurlburt84-Ra1e5.png $ns1/f90/pencil-code/www/samples/2d-tests
;
window,xs=900,ys=360
!p.charsize=2
!p.charthick=2 & !p.thick=2 & !x.thick=2 & !y.thick=2
!x.title='!8x!6'
!y.title='!8z!6'
contour,-reform(ss(l1:l2,3,n1:n2)),x(l1:l2),z(n1:n2),/fil,nlev=60,back=255,col=0
;velovect,reform(uu(l1:l2,3,n1:n2,0)),reform(uu(l1:l2,3,n1:n2,2)),x(l1:l2),z(n1:n2),len=2,/over
!p.charthick=1 & !p.thick=1 & !x.thick=1 & !y.thick=1
vel_a,reform(uu(l1:l2,3,n1:n2,0)),reform(uu(l1:l2,3,n1:n2,2)),x(l1:l2),z(n1:n2),/over,nvec=1800,col=255
print,'import Hurlburt84-Ra1e5.png'
END
