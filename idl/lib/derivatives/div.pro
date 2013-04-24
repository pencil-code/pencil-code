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
  corr=fltarr(mx,my,mz)
;
  if (coord_system eq 'cartesian') then corr = 0.
;
  if (coord_system eq 'cylindric') then begin
    for m=m1,m2 do begin
      for n=n1,n2 do begin   
        corr[l1:l2,m,n] = f[l1:l2,m,n,0]/x[l1:l2]
      endfor
    endfor  
  endif
;
  if (coord_system eq 'spherical') then begin
    cotth=cos(y)/sin(y)      
    threshold = 100.*num_model(y, /EPSILON) ; sinth_min=1e-5/2e-14
    i_sin=where(abs(sin(y)) lt threshold)
    if (i_sin[0] ne -1) then cotth[i_sin]=0.
    for m=m1,m2 do begin
      for n=n1,n2 do begin   
        corr[l1:l2,m,n] = (2.*f[l1:l2,m,n,0]+cotth*f[l1:l2,m,n,1])/x[l1:l2]
      endfor
    endfor 
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
