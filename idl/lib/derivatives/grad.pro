;;; Calculate the gradient of a 3-d scalar field
function grad,f
COMPILE_OPT IDL2,HIDDEN
common cdat,x
s=size(f) & nx=s[1] & ny=s[2] & nz=s[3]
if s[0] eq 2 then nz=1
 w=fltarr(nx,ny,nz,3)

 w[*,*,*,0]=xder(f)
; w[*,*,*,1]=yder(f)
 w[*,*,*,2]=zder(f)

 pc_read_param,obj=par
 if (par.coord_system eq 'cylindric') then begin
     tmp=size(f) &  my=tmp[2] &  mz=tmp[3]
     xx=spread(x,[1,2],[my,mz])
     w[*,*,*,1]=yder(f)/xx
 endif else if (par.coord_system eq 'cartesian') then begin
     w[*,*,*,1]=yder(f)
 endif else begin
     print, 'error: grad not implemented for spherical polars'
 endelse
;
 return,w
end
