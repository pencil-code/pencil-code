;***********************************************************************
function xder2,f
;
common cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0 
;
;  Check if we have a degenerate case (no x-extension)
;
if (nx eq 1) then return,fltarr(nx,ny,nz)
s=size(f) & d=make_array(size=s)
;
;  assume uniform mesh
;
dx2=1./(180.*(x(4)-x(3))^2)
l1=3 & l2=nx-4
;
if s(0) eq 1 then begin
  print,'not implemented yet'
end else if s(0) eq 2 then begin
  print,'not implemented yet'
end else if s(0) eq 3 then begin
  d(l1:l2,*,*)=dx2*(-490.*f(l1:l2,*,*)$
                   +270.*(f(l1-1:l2-1,*,*)+f(l1+1:l2+1,*,*))$
                    -27.*(f(l1-2:l2-2,*,*)+f(l1+2:l2+2,*,*))$
                     +2.*(f(l1-3:l2-3,*,*)+f(l1+3:l2+3,*,*))$
        )
end
;
return,d
end
