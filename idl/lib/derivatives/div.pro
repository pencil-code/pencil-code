;;
;;  $Id$
;;
;;  Calculate divergence of vector field.
;;
function div,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat,x,y,z
  common cdat_coords,coord_system
;
;  Default values.
;
  default, ghost, 0
;
  tmp=size(f)
  mx=tmp[1] & my=tmp[2] & mz=tmp[3]
;
  if (coord_system eq 'cartesian') then corr = 0.
;
  if (coord_system eq 'cylindric') then begin
    xx=spread(x,[1,2],[my,mz])
    corr = f[*,*,*,0]/xx
  endif
;
  if (coord_system eq 'spherical') then begin
    xx=spread(x,[1,2],[my,mz])
    yy=spread(y,[0,2],[mx,mz])

    cotth=cos(yy)/sin(yy)      
    threshold = 100.*num_model(yy, /EPSILON) ; sinth_min=1e-5/2e-14
    i_sin=where(abs(sin(yy)) lt threshold)
    if (i_sin[0] ne -1) then cotth[i_sin]=0.
    corr = (2.*f[*,*,*,0]+cotth*f[*,*,*,1])/xx
  endif
;
  w=xder(f[*,*,*,0]) + yder(f[*,*,*,1]) + zder(f[*,*,*,2]) + corr
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
