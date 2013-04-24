;;
;;  $Id$
;;
;;  Calculate the curl of a 3-D vector field
;;
function curlx,f,coord_system,x,y,d
;
  COMPILE_OPT IDL2,HIDDEN
;
  if (coord_system eq 'cartesian') then corr=0.
  if (coord_system eq 'cylindric') then corr=0.
  if (coord_system eq 'spherical') then begin
    cotth=cos(y)/sin(y)      
    i_sin=where(abs(sin(y)) lt 1e-5) ;sinth_min=1e-5
    if (i_sin[0] ne -1) then cotth[i_sin]=0.
    l1=d.l1 & l2=d.l2 & m1=d.m1
    m2=d.m2 & n1=d.n1 & n2=d.n2
    tmp=size(f) & mx=tmp[1] & my=tmp[2] &  mz=tmp[3]
    corr=fltarr(mx,my,mz)
    for m=m1,m2 do begin
      for n=n1,n2 do begin   
        corr[l1:l2,m,n]=f[l1:l2,m,n,2]*cotth[m]/x[l1:l2]
      endfor
    endfor  
  endif
  return,yder(f[*,*,*,2])-zder(f[*,*,*,1])+corr
end
;;
function curly,f,coord_system,x,d
;
  COMPILE_OPT IDL2,HIDDEN
;
  if (coord_system eq 'cartesian') then corr=0.
  if (coord_system eq 'cylindric') then corr=0.
  if (coord_system eq 'spherical') then begin
    l1=d.l1 & l2=d.l2 & m1=d.m1
    m2=d.m2 & n1=d.n1 & n2=d.n2
    tmp=size(f) & mx=tmp[1] & my=tmp[2] &  mz=tmp[3] 
    corr=fltarr(mx,my,mz)
    for m=m1,m2 do begin
      for n=n1,n2 do begin   
        corr[l1:l2,m,n]=-f[l1:l2,m,n,2]/x[l1:l2]
      endfor
    endfor  
  endif   
  return,zder(f[*,*,*,0])-xder(f[*,*,*,2])+corr
;
end
;;
function curlz,f,coord_system,x,d
;
  COMPILE_OPT IDL2,HIDDEN
;
  if (coord_system eq 'cartesian') then corr=0.
  if (coord_system eq 'cylindric' or $
      coord_system eq 'spherical') then begin
    l1=d.l1 & l2=d.l2 & m1=d.m1 
    m2=d.m2 & n1=d.n1 & n2=d.n2 
    tmp=size(f) & mx=tmp[1] & my=tmp[2] &  mz=tmp[3]
    corr=fltarr(mx,my,mz)
    for m=m1,m2 do begin
      for n=n1,n2 do begin   
        corr[l1:l2,m,n]=f[l1:l2,m,n,1]/x[l1:l2]
      endfor
    endfor  
  endif   
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
  default, nghost, 3
;
  w=make_array(size=size(f))
;
  tmp=size(f) & mx=tmp[1] & my=tmp[2] &  mz=tmp[3]
  l1=nghost & l2=mx-nghost-1
  m1=nghost & m2=my-nghost-1
  n1=nghost & n2=mz-nghost-1
;
  d={l1:l1,l2:l2,m1:m1,m2:m2,n1:n1,n2:n2,mx:mx,my:my,mz:mz}
;
  w[*,*,*,0]=curlx(f,coord_system,x,y,d)
  w[*,*,*,1]=curly(f,coord_system,x,d)
  w[*,*,*,2]=curlz(f,coord_system,x,d)
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
