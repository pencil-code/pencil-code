;***********************************************************************
function yder,f
COMPILE_OPT IDL2,HIDDEN
;
common cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0 
common cdat_nonequidist,xprim,yprim,zprim,xprim2,yprim2,zprim2,lequidist
;
;  Check if we have a degenerate case (no y-extension)
;
if n_elements(lequidist) ne 3 then lequidist=[1,1,1]
if (ny eq 1) then return,fltarr(nx,ny,nz)
s=size(f) & d=make_array(size=s)
;
;  assume uniform mesh
;
m1=3 & m2=ny-4
;
if lequidist[2] then begin
  dy2=1./(60.*(y[4]-y[3]))
endif else begin
  dy2=spread(spread(1./(60.*yprim),0,nx),2,nz)
endelse
;
  d[*,m1:m2,*]=dy2*(+45.*(f[*,m1+1:m2+1,*]-f[*,m1-1:m2-1,*])$
                     -9.*(f[*,m1+2:m2+2,*]-f[*,m1-2:m2-2,*])$
                        +(f[*,m1+3:m2+3,*]-f[*,m1-3:m2-3,*])$
      )
;
return,d
end
