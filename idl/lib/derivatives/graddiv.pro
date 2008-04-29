;
;  $Id: graddiv.pro,v 1.2 2008-04-29 22:13:08 dobler Exp $
;
;  Calculate gradient of the divergence of a vector.
;
function graddiv,f
  COMPILE_OPT IDL2,HIDDEN
  common cdat_coords, coord_system

  if (coord_system ne 'cartesian') then $
      message, $
        "graddiv not yet implemented for coord_system='" + coord_system + "'"

  s=size(f)
;
  if (s[0] eq 4) then begin

;
    w=make_array(n_elements(f[*,0,0,0]),n_elements(f[0,*,0,0]),n_elements(f[0,0,*,0]),3,/nozero)
    w[*,*,*,0]=   xder2(f[*,*,*,0])+xderyder(f[*,*,*,1])+xderzder(f[*,*,*,2])
    w[*,*,*,1]=xderyder(f[*,*,*,0])+   yder2(f[*,*,*,1])+yderzder(f[*,*,*,2])
    w[*,*,*,2]=xderzder(f[*,*,*,0])+yderzder(f[*,*,*,1])+   zder2(f[*,*,*,2])
;
  endif else begin
    print, 'error: graddiv only implemented for 4-D arrays'
  endelse
;
  return, w
;
end
