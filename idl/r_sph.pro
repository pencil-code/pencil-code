; $Id$

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
sinth = spread(sin(y),[0,2],[mx,mz])
sinth1 = 1./sinth
END

; End of file
