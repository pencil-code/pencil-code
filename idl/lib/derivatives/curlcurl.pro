;;
;;  $Id$
;;
;;  Calculate two consecutive curls of a 3-D vector field.
;;
;;  18-jan-08/anders: coded
;;
function curlcurl,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2, HIDDEN
;
;  Default values.
;
  default, ghost, 0
;
  w = graddiv(f) - del2(f)
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
