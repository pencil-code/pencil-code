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
  return,yder(f[*,*,*,2])-zder(f[*,*,*,1])+corr
end
function curly,f,lsystem,xx
COMPILE_OPT IDL2,HIDDEN
  if (lsystem eq 0) then corr=0.
  if (lsystem eq 1) then corr=0.
  if (lsystem eq 2) then corr=-f[*,*,*,2]/xx
  return,zder(f[*,*,*,0])-xder(f[*,*,*,2])+corr
end
function curlz,f,lsystem,xx
COMPILE_OPT IDL2,HIDDEN
  if (lsystem eq 0) then corr=0.
  if (lsystem eq 1) then corr=f[*,*,*,1]/xx
  if (lsystem eq 2) then corr=f[*,*,*,1]/xx
  return,xder(f[*,*,*,1])-yder(f[*,*,*,0])+corr
end
;
function curl,f

COMPILE_OPT IDL2,HIDDEN
common cdat, x, y
;
  w=make_array(size=size(f),/nozero)
  pc_read_param,obj=par,datadir=datadir,dim=dim,/quiet
  if (par.coord_system eq 'cartesian') then lsystem=0
  if (par.coord_system eq 'cylindric') then lsystem=1
  if (par.coord_system eq 'spherical') then lsystem=2
;
  xx=spread(x,[1,2],[dim.my,dim.mz])
  yy=spread(y,[0,2],[dim.mx,dim.mz])
;
  w[*,*,*,0]=curlx(f,lsystem,xx,yy)
  w[*,*,*,1]=curly(f,lsystem,xx)
  w[*,*,*,2]=curlz(f,lsystem,xx)
;
  return,w
;
end
