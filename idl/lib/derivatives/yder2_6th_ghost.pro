;***********************************************************************
function yder2,f
COMPILE_OPT IDL2,HIDDEN
;
common cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0 
common cdat_nonequidist,xprim,yprim,zprim,xprim2,yprim2,zprim2,lequidist
;
;  Check if we have a degenerate case (no x-extension)
;
if n_elements(lequidist) ne 3 then lequidist=[1,1,1]
if (ny eq 1) then return,fltarr(nx,ny,nz)
s=size(f) & d=make_array(size=s)
;
;  assume uniform mesh
;
m1=3 & m2=ny-4
;
if lequidist[1] then begin
  dy2=1./(180.*(y[4]-y[3])^2)
endif else begin
  d1=yder(f)
;
;  nonuniform mesh correction
;  d2f/dy2  = f"*psi'^2 + psi"f'   see also the manual
;
  dy2=spread(spread(yprim^2/180.,0,nx),2,nz)
  dd =d1*spread(spread(yprim2,0,nx),2,nz)
endelse
;
if s[0] eq 2 then begin
  print,'not implemented yet'
end else if s[0] eq 3 then begin
  d[*,m1:m2,*]=dy2*(-490.*f[*,m1:m2,*]$
                   +270.*(f[*,m1-1:m2-1,*]+f[*,m1+1:m2+1,*])$
                    -27.*(f[*,m1-2:m2-2,*]+f[*,m1+2:m2+2,*])$
                     +2.*(f[*,m1-3:m2-3,*]+f[*,m1+3:m2+3,*])$
                   )
  if not lequidist[1] then d=d+dd
  
end
;
return,d
end
