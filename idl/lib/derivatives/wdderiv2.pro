;;;;;;;;;;;;;;;;;;;;;;;;
;;;   wdderiv2.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   $Date: 2004-05-05 16:42:57 $
;;;
;;;  Description:
;;;   Second derivative function like IDL's deriv().

function wdderiv2, arg1, arg2
COMPILE_OPT IDL2,HIDDEN

  n1 = n_elements(arg1)
  n2 = n_elements(arg2)

  if (n2 gt 0) then begin
    if (n1 ne n2) then message, "x and y must be of same size"
    x = arg1
    y = arg2
  endif else if (n1 gt 0) then begin
    x = indgen(n1)
    y = arg1
  endif else begin
    message, "Usage: res = wdderiv2([x,] y)"
  endelse

  x_l = shift(x, 1)
  x_r = shift(x,-1)
  d1_l = (y-shift(y, 1))/(x-x_l)
  d1_r = (y-shift(y,-1))/(x-x_r)

  ;; Inner points
  der2 = 2/(x_r-x_l)*(d1_r-d1_l)

  ;; Boundary points by extrapolation
  der2[0]  = der2[1] + (x[0]-x[1])*(der2[2]-der2[1])/(x[2]-x[1])
  der2[n1-1] = der2[n1-2] $
               + (x[n1-1]-x[n1-2])*(der2[n1-3]-der2[n1-2])/(x[n1-3]-x[n1-2])

  return, der2;

end
; End of file wdderiv2.pro
