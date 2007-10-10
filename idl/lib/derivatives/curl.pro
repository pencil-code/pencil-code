;;; Calculate the curl of a 3-d vector field
function curlx,f
COMPILE_OPT IDL2,HIDDEN
  return,yder(f[*,*,*,2])-zder(f[*,*,*,1])
end
function curly,f
COMPILE_OPT IDL2,HIDDEN
  return,zder(f[*,*,*,0])-xder(f[*,*,*,2])
end
function curlz,f
COMPILE_OPT IDL2,HIDDEN
  return,xder(f[*,*,*,1])-yder(f[*,*,*,0])
end
function curlrad,f,xx
COMPILE_OPT IDL2,HIDDEN
  return,yder(f[*,*,*,2])/xx -zder(f[*,*,*,1])
end
function curlphi,f,xx
COMPILE_OPT IDL2,HIDDEN
  return,zder(f[*,*,*,0])-xder(f[*,*,*,2])
end
function curlzed,f,xx
COMPILE_OPT IDL2,HIDDEN
  return,xder(f[*,*,*,1])+f[*,*,*,1]/xx-yder(f[*,*,*,0])/xx
end
function curl,f

COMPILE_OPT IDL2,HIDDEN
common cdat, x
;
  w=make_array(size=size(f),/nozero)
;
  pc_read_param,obj=par,/quiet
;
  if (par.coord_system eq 'cylindric') then begin
    tmp=size(f) &  my=tmp[2] &  mz=tmp[3]
    xx=spread(x,[1,2],[my,mz])
    w[*,*,*,0]=curlrad(f,xx)
    w[*,*,*,1]=curlphi(f,xx)
    w[*,*,*,2]=curlzed(f,xx)
  endif else if (par.coord_system eq 'cartesian') then begin
    w[*,*,*,0]=curlx(f)
    w[*,*,*,1]=curly(f)
    w[*,*,*,2]=curlz(f)
  endif else begin
    print, 'error: curl not implemented for spherical polars'
  endelse
;
  return,w
;
end
