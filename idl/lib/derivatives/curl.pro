;;; Calculate the curl of a 3-d vector field
;cartesian
function curlx,f,lsystem,xx,yy
COMPILE_OPT IDL2,HIDDEN
  if (lsystem eq 0) then corr=0.
  if (lsystem eq 1) then corr=0.
  if (lsystem eq 2) then begin
    cotth=cos(yy)/sin(yy)      
    i_sin=where(abs(sin(yy)) lt 1e-5) ;sinth_min=1e-5
    if (i_sin ne -1) then cotth[i_sin]=0.
    corr=f[*,*,*,2]*cotth/xx
  endif
  return,yder(f[*,*,*,2],lsystem)-zder(f[*,*,*,1],lsystem)+corr
end
function curly,f,lsystem,xx
COMPILE_OPT IDL2,HIDDEN
  if (lsystem eq 0) then corr=0.
  if (lsystem eq 1) then corr=0.
  if (lsystem eq 2) then corr=-f[*,*,*,2]/xx
  return,zder(f[*,*,*,0],lsystem)-xder(f[*,*,*,2])+corr
end
function curlz,f,lsystem,xx
COMPILE_OPT IDL2,HIDDEN
  if (lsystem eq 0) then corr=0.
  if (lsystem eq 1) then corr=f[*,*,*,1]/xx
  if (lsystem eq 2) then corr=f[*,*,*,1]/xx
  return,xder(f[*,*,*,1])-yder(f[*,*,*,0],lsystem)+corr
end
;
function curl,f,coord_system

COMPILE_OPT IDL2,HIDDEN
common cdat, x, y
;
default,coord_system,'cartesian'
;
  w=make_array(size=size(f),/nozero)
  lsystem=-1
  if (coord_system eq 'cartesian') then lsystem=0
  if (coord_system eq 'cylindric') then lsystem=1
  if (coord_system eq 'spherical') then lsystem=2
  if (lsystem eq -1) then $
    print,'coord_system= ',coord_system,' is not valid'
;
  tmp=size(f) & mx=tmp[1] & my=tmp[2] &  mz=tmp[3]
  xx=spread(x,[1,2],[my,mz])
  yy=spread(y,[0,2],[mx,mz])
;
  w[*,*,*,0]=curlx(f,lsystem,xx,yy)
  w[*,*,*,1]=curly(f,lsystem,xx)
  w[*,*,*,2]=curlz(f,lsystem,xx)
;
  return,w
;
end
