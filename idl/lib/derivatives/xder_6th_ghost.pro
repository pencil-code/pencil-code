;***********************************************************************
function xder,f
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
dx2=1./(60.*(x(4)-x(3)))
l1=3 & l2=nx-4
s=size(f)
;
if s(0) eq 1 then begin
end else if s(0) eq 2 then begin
end else if s(0) eq 3 then begin
  d(l1:l2,*,*)=dx2*(+45.*(f(l1+1:l2+1,*,*)-f(l1-1:l2-1,*,*))$
                     -9.*(f(l1+2:l2+2,*,*)-f(l1-2:l2-2,*,*))$
                        +(f(l1+3:l2+3,*,*)-f(l1-3:l2-3,*,*))$
        )
end
;
return,d
end
