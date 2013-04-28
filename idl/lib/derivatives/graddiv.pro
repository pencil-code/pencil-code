;;
;;  $Id$
;;
;;  Calculate gradient of the divergence of a vector.
;;
function graddiv,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat_coords, coord_system
;
;  Default values.
;
  default, ghost, 0
;
  if (coord_system ne 'cartesian') then message, $
      "graddiv not yet implemented for coord_system='" + coord_system + "'"
;
  s = size(f)
  if (s[0] ne 4) then message, "grad is only implemented for 4D arrays."
  w = make_array(size=s)
;
  w[*,*,*,0]=   xder2(f[*,*,*,0])+xderyder(f[*,*,*,1])+xderzder(f[*,*,*,2])
  w[*,*,*,1]=xderyder(f[*,*,*,0])+   yder2(f[*,*,*,1])+yderzder(f[*,*,*,2])
  w[*,*,*,2]=xderzder(f[*,*,*,0])+yderzder(f[*,*,*,1])+   zder2(f[*,*,*,2])
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
