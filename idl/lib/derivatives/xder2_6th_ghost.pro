;***********************************************************************
function xder2,f
COMPILE_OPT IDL2,HIDDEN
;
common cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0 
common cdat_nonequidist,xprim,yprim,zprim,xprim2,yprim2,zprim2,lequidist
;
;  Check if we have a degenerate case (no x-extension)
;

if n_elements(lequidist) ne 3 then lequidist=[1,1,1]
if (nx eq 1) then return,fltarr(nx,ny,nz)
s=size(f) & d=make_array(size=s)

l1=3 & l2=nx-4
;
if lequidist[0] then begin
  dx2=1./(180.*(x[4]-x[3])^2)
endif else begin
  d1=xder(f)
  dx2=spread(spread(1./(180.*xprim^2),1,ny),2,nz)-d1*spread(spread(xprim2/xprim^2,1,ny),2,nz)
endelse

;
if s[0] eq 1 then begin
  print,'not implemented yet'
end else if s[0] eq 2 then begin
  print,'not implemented yet'
end else if s[0] eq 3 then begin
  d[l1:l2,*,*]=dx2*(-490.*f[l1:l2,*,*]$
                   +270.*(f[l1-1:l2-1,*,*]+f[l1+1:l2+1,*,*])$
                    -27.*(f[l1-2:l2-2,*,*]+f[l1+2:l2+2,*,*])$
                     +2.*(f[l1-3:l2-3,*,*]+f[l1+3:l2+3,*,*])$
        )
end

;
return,d
end
