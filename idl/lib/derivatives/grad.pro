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
  s=size(f) & nx=s[1] & ny=s[2] & nz=s[3]
  if s[0] eq 2 then nz=1
  w=fltarr(nx,ny,nz,3)
;
  w[*,*,*,0]=xder(f)
  w[*,*,*,1]=yder(f)
  w[*,*,*,2]=zder(f)
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
