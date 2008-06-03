;
;  $Id: del6.pro,v 1.1 2008-06-03 13:49:16 ajohan Exp $
;
;  Calculate del6 of f.
;
;  03-jun-08/anders: coded
;
function del6,f
  COMPILE_OPT IDL2,HIDDEN
  common cdat_coords, coord_system
;
  if (coord_system ne 'cartesian') then message, $
      "del6 not yet implemented for coord_system='" + coord_system + "'"
;
  s=size(f)
;
  if (s[0] eq 3) then begin
;
    w=make_array(n_elements(f[*,0,0]),n_elements(f[0,*,0]),n_elements(f[0,0,*]),3,/nozero)
    w=xder6(f)+yder6(f)+zder6(f)
;
  endif else if (s[0] eq 4) then begin
;
    w=make_array(n_elements(f[*,0,0,0]),n_elements(f[0,*,0,0]),n_elements(f[0,0,*,0]),3,/nozero)
    w=xder6(f)+yder6(f)+zder6(f)
;
  endif else begin
    print, 'error: del6 not implemented for arrays of size ', s
    message, 'no point in continuing'
  endelse
;
  return, w
;
end
