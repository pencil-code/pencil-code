;***********************************************************************
function zder2,f
COMPILE_OPT IDL2,HIDDEN
;
common cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0 
;
;  Check if we have a degenerate case (no x-extension)
;
if (nz eq 1) then return,fltarr(nx,ny,nz)
s=size(f) & d=make_array(size=s)
;
;  assume uniform mesh
;
dz2=1./(180.*(z[4]-z[3])^2)
n1=3 & n2=nz-4
;
  d[*,*,n1:n2]=dz2*(-490.*f[*,*,n1:n2]$
                   +270.*(f[*,*,n1-1:n2-1]+f[*,*,n1+1:n2+1])$
                    -27.*(f[*,*,n1-2:n2-2]+f[*,*,n1+2:n2+2])$
                     +2.*(f[*,*,n1-3:n2-3]+f[*,*,n1+3:n2+3])$
      )
;
return,d
end
