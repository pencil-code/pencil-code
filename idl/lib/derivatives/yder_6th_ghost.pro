;***********************************************************************
function yder,f
COMPILE_OPT IDL2,HIDDEN
;
common cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0 
;
;  Check if we have a degenerate case (no y-extension)
;
if (ny eq 1) then return,fltarr(nx,ny,nz)
s=size(f) & d=make_array(size=s)
;
;  assume uniform mesh
;
dy2=1./(60.*(y[4]-y[3]))
m1=3 & m2=ny-4
;
  d[*,m1:m2,*]=dy2*(+45.*(f[*,m1+1:m2+1,*]-f[*,m1-1:m2-1,*])$
                     -9.*(f[*,m1+2:m2+2,*]-f[*,m1-2:m2-2,*])$
                        +(f[*,m1+3:m2+3,*]-f[*,m1-3:m2-3,*])$
      )
;
return,d
end
