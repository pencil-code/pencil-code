;***********************************************************************
function zder,f
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
dz2=1./(60.*(z[4]-z[3]))
n1=3 & n2=nz-4
;
if n2 gt n1 then begin
  d[*,*,n1:n2]=dz2*(+45.*(f[*,*,n1+1:n2+1]-f[*,*,n1-1:n2-1])$
                     -9.*(f[*,*,n1+2:n2+2]-f[*,*,n1-2:n2-2])$
                        +(f[*,*,n1+3:n2+3]-f[*,*,n1-3:n2-3])$
      )
endif else begin
  d[*,*,n1:n2]=0.
endelse
;
return,d
end
