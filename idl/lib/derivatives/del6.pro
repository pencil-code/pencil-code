;;
;;  $Id$
;;
;;  Calculate del6 of f.
;;
;;  03-jun-08/anders: coded
;;
function del6,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t, $
              ignoredx=ignoredx
  COMPILE_OPT IDL2,HIDDEN
;
;  Default values.
;
  default, ghost, 0
;
  w = xder6(f,ignoredx=ignoredx) + yder6(f,ignoredx=ignoredx) + zder6(f,ignoredx=ignoredx)
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
