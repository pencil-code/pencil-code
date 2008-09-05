;;
;;  $Id$
;;
;;  Calculate the curl of a 3-D vector field
;;
function curlx,f,coord_system,xx,yy
;
  COMPILE_OPT IDL2,HIDDEN
;
  if (coord_system eq 'cartesian') then corr=0.
  if (coord_system eq 'cylindric') then corr=0.
  if (coord_system eq 'spherical') then begin
    cotth=cos(yy)/sin(yy)      
    i_sin=where(abs(sin(yy)) lt 1e-5) ;sinth_min=1e-5
    if (i_sin[0] ne -1) then cotth[i_sin]=0.
    corr=f[*,*,*,2]*cotth/xx
  endif
  return,yder(f[*,*,*,2])-zder(f[*,*,*,1])+corr
end
;;
function curly,f,coord_system,xx
;
  COMPILE_OPT IDL2,HIDDEN
;
  if (coord_system eq 'cartesian') then corr=0.
  if (coord_system eq 'cylindric') then corr=0.
  if (coord_system eq 'spherical') then corr=-f[*,*,*,2]/xx
;
  return,zder(f[*,*,*,0])-xder(f[*,*,*,2])+corr
;
end
;;
function curlz,f,coord_system,xx
;
  COMPILE_OPT IDL2,HIDDEN
;
  if (coord_system eq 'cartesian') then corr=0.
  if (coord_system eq 'cylindric') then corr=f[*,*,*,1]/xx
  if (coord_system eq 'spherical') then corr=f[*,*,*,1]/xx
;
  return,xder(f[*,*,*,1])-yder(f[*,*,*,0])+corr
end
;;
function curl,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat,x,y
  common cdat_coords,coord_system
;
;  Default values.
;
  default, ghost, 0
;
  w=make_array(size=size(f))
;
  tmp=size(f) & mx=tmp[1] & my=tmp[2] &  mz=tmp[3]
  xx=spread(x,[1,2],[my,mz])
  yy=spread(y,[0,2],[mx,mz])
;
  w[*,*,*,0]=curlx(f,coord_system,xx,yy)
  w[*,*,*,1]=curly(f,coord_system,xx)
  w[*,*,*,2]=curlz(f,coord_system,xx)
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
