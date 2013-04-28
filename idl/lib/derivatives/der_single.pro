;;
;;  $Id$
;;
;;  Calculate single component of firstst derivative; the second argument
;;  gives the direction.
;;
;;  20-mar-04/axel: coded
;;
function der_single,f,j,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
;
;  Default values.
;
  default, ghost, 0
;
  if (j eq 0) then w = xder(f)
  if (j eq 1) then w = yder(f)
  if (j eq 2) then w = zder(f)
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
