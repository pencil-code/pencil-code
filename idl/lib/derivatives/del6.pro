;;
;;  $Id: del6.pro,v 1.3 2008-06-10 17:24:36 ajohan Exp $
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
  common cdat_coords, coord_system
;
  if (coord_system ne 'cartesian') then message, $
      "del6 not yet implemented for coord_system='" + coord_system + "'"
;
  s=size(f)
;
  if (s[0] eq 3) then begin
;
    w=make_array(n_elements(f[*,0,0]),n_elements(f[0,*,0]),n_elements(f[0,0,*]),3)
    w=xder6(f)+yder6(f)+zder6(f)
;
  endif else if (s[0] eq 4) then begin
;
    w=make_array(n_elements(f[*,0,0,0]),n_elements(f[0,*,0,0]),n_elements(f[0,0,*,0]),3)
    w=xder6(f)+yder6(f)+zder6(f)
;
  endif else begin
    print, 'error: del6 not implemented for arrays of size ', s
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
