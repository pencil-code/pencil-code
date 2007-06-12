; $Id: r_sph.pro,v 1.1 2007-06-12 16:16:57 dhruba Exp $

;;;;;;;;;;;;;;;;;;;
;;;  r_sph.pro  ;;;
;;;;;;;;;;;;;;;;;;;

;;; command for spherical polar coordinate system
;;; to be use after r.pro or rall.pro



rr = spread(x, [1,2], [my,mz])
rr_1 = spread(1/x,[1,2],[my,mz])
theta = yy
phi = zz
cotth = spread(cos(y)/sin(y),[0,2],[mx,mz]) 

END

; End of file
