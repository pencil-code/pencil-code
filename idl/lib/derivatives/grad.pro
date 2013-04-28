;;
;;  $Id$
;;
;;  Calculate the gradient of a 3-D scalar field.
;;
function grad,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
;  Default values.
;
  default, ghost, 0
;
  s = size(f)
  if (s[0] ne 3) then message, "grad is only implemented for 3D arrays."
  fmx = s[1] & fmy = s[2] & fmz = s[3]
  w = make_array(size=[4,fmx,fmy,fmz,3,s[4],fmx*fmy*fmz*3])
;
  w[*,*,*,0] = xder(f)
  w[*,*,*,1] = yder(f)
  w[*,*,*,2] = zder(f)
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
