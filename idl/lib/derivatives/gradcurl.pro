;
;  $Id$
;
;  Gradient from curl in cartesian coordinates.
;  Used to calculate grad(B) as grad(curl(A)).
;
;  20-Jan-2017/PABourdin: coded
;
function gradcurl,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat_coords, coord_system
;
;  Default values.
;
  default, ghost, 0
;
  if (coord_system ne 'cartesian') then message, $
      "gradcurl not yet implemented for coord_system='" + coord_system + "'"
;
  s = size(f)
  if (s[0] ne 4) then message, "gradcurl is only implemented for 4D arrays."
  w = make_array(size=s)
;
;  Calculate grad(curl(A)) as:
;  ( dx*dy*Az-dx*dz*Ay, dy*dx*Az-dy*dz*Ax, dz*dx*Ay-dz*dy*Ax )
;
  w[*,*,*,0] = derij_single(f[*,*,*,2],0,1) - derij_single(f[*,*,*,1],0,2)
  w[*,*,*,1] = derij_single(f[*,*,*,2],1,0) - derij_single(f[*,*,*,0],1,2)
  w[*,*,*,2] = derij_single(f[*,*,*,1],2,0) - derij_single(f[*,*,*,0],2,1)
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
