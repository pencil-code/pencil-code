;;
;;  $Id$
;;
;;  Calculate single component of second derivative matrix.
;;
;;  20-mar-04/axel: coded
;;
function derij_single,f,i,j,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
;
;  Default values.
;
  default, ghost, 0
;
  if ((i eq 0) and (j eq 0)) then w = xder2(f)
  if ((i eq 0) and (j eq 1)) then w = xderyder(f)
  if ((i eq 0) and (j eq 2)) then w = xderzder(f)
  if ((i eq 1) and (j eq 0)) then w = yderxder(f)
  if ((i eq 1) and (j eq 1)) then w = yder2(f)
  if ((i eq 1) and (j eq 2)) then w = yderzder(f)
  if ((i eq 2) and (j eq 0)) then w = zderxder(f)
  if ((i eq 2) and (j eq 1)) then w = zderyder(f)
  if ((i eq 2) and (j eq 2)) then w = zder2(f)
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
