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
  common cdat_coords, coord_system
;
;  Default values.
;
  default, ghost, 0
;
  if (coord_system ne 'cartesian') then message, $
      "del2 not yet implemented for coord_system='" + coord_system + "'"
;
  s=size(f)
;
  if ((s[0] ge 2) and (s[0] le 4)) then begin
;
    w=make_array(size=s)
    w=xder2(f)+yder2(f)+zder2(f)
;
  endif else begin
    print, 'error: del2 not implemented for arrays of size ', s
    message, 'no point in continuing'
  endelse
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
