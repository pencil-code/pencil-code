;;;;;;;;;;;;;;;
;;;  r.pro  ;;;
;;;;;;;;;;;;;;;

;;; Read the data produced by one processor with the fortran minimal
;;; hydro code.
;;; Assumes you have run `start.pro' once before.

;
;  read data
;
default, datadir, 'tmp'
default, file, 'var.dat'
;
uu  = fltarr(nx,ny,nz,3)*one
lam = fltarr(nx,ny,nz)*one
ent = fltarr(nx,ny,nz)*one
aa  = fltarr(nx,ny,nz,3)*one
;
close,1
openr,1, datadir+'/'+file, /F77
readu,1, uu, lam, ent, aa
readu,1, t, x, y, z
close,1
;
xx = spread(x, [1,2], [ny,nz])
yy = spread(y, [0,2], [nx,nz])
zz = spread(z, [0,1], [nx,ny])
rr = sqrt(xx^2+yy^2+zz^2)
;
print, ' var        minval       maxval'
for j=0,2 do print, 'uu_', strtrim(j,2), ' =', minmax(uu(*,*,*,j))
print, 'lam  =', minmax(lam)
print, 'ent  =', minmax(ent)
for j=0,2 do print, 'aa_', strtrim(j,2), ' =', minmax(aa(*,*,*,j))
;
print,'t = ',t
;
END

; End of file
