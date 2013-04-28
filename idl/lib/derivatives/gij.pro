;;
;;  $Id$
;;
;;  Calculate derivative matrix.
;;
function gij,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat_coords, coord_system
;
;  Default values.
;
  default, ghost, 0
;
  if (coord_system ne 'cartesian') then $
      message, "gij not yet implemented for coord_system='"+coord_system+"'"
;
  s = size(f)
  if (s[0] ne 4) then message, "gij is only implemented for 4D arrays."
  fmx = s[1] & fmy = s[2] & fmz = s[3]
  w = make_array(size=[5,fmx,fmy,fmz,s[4],3,s[5],fmx*fmy*fmz*s[4]*3])
;
  w[*,*,*,*,0] = xder(f)
  w[*,*,*,*,1] = yder(f)
  w[*,*,*,*,2] = zder(f)
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
