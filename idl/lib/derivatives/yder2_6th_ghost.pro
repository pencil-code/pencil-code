;***********************************************************************
function yder2,f
;
common cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0 
;
;  Check if we have a degenerate case (no x-extension)
;
if (ny eq 1) then return,fltarr(nx,ny,nz)
s=size(f) & d=make_array(size=s)
;
;  assume uniform mesh
;
dy2=1./(180.*(y(4)-y(3))^2)
m1=3 & m2=ny-4
;
if s(0) eq 2 then begin
  print,'not implemented yet'
end else if s(0) eq 3 then begin
  d(*,m1:m2,*)=dy2*(-490.*f(*,m1:m2,*,*)$
                   +270.*(f(*,m1-1:m2-1,*)+f(*,m1+1:m2+1,*))$
                    -27.*(f(*,m1-2:m2-2,*)+f(*,m1+2:m2+2,*))$
                     +2.*(f(*,m1-3:m2-3,*)+f(*,m1+3:m2+3,*))$
        )
end
;
return,d
end
