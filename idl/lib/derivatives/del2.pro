;;
;;  $Id$
;;
;;  Calculate the Laplacian of f, i.e.
;;    div(grad(f)) if f is a scalar field, or
;;    grad(div(f)) - curl(curl(f)) if f is a vector field.
;;
;;  18-jan-08/anders: coded
;;
function del2,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
;  Default values.
;
  default, ghost, 0
;
  w = xder2(f) + yder2(f) + zder2(f)
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
