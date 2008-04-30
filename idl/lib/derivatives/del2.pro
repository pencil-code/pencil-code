;
;  $Id: del2.pro,v 1.4 2008-04-30 15:00:10 ajohan Exp $
;
;  Calculate the Laplacian of f, i.e.
;    div(grad(f)) if f is a scalar field, or
;    grad(div(f)) - curl(curl(f)) if f is a vector field.
;
;  18-jan-08/anders: coded
;
function del2,f
  COMPILE_OPT IDL2,HIDDEN
  common cdat_coords, coord_system
;
  if (coord_system ne 'cartesian') then message, $
      "del2 not yet implemented for coord_system='" + coord_system + "'"
;
  s=size(f)
;
  if (s[0] eq 3) then begin
;
    w=make_array(n_elements(f[*,0,0]),n_elements(f[0,*,0]),n_elements(f[0,0,*]),3,/nozero)
    w=xder2(f)+yder2(f)+zder2(f)
;
  endif else if (s[0] eq 4) then begin
;
    w=make_array(n_elements(f[*,0,0,0]),n_elements(f[0,*,0,0]),n_elements(f[0,0,*,0]),3,/nozero)
    w=xder2(f)+yder2(f)+zder2(f)
;
  endif else begin
    print, 'error: del2 not implemented for arrays of size ', s
    message, 'no point in continuing'
  endelse
;
  return, w
;
end
