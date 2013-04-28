;;
;;  $Id$
;;
;;  Calculate del6 of f.
;;
;;  03-jun-08/anders: coded
;;
function del6,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
;  Default values.
;
  default, ghost, 0
;
  w = xder6(f) + yder6(f) + zder6(f)
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
