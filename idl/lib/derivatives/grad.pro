;;; Calculate the gradient of a 3-d scalar field
function grad,f
COMPILE_OPT IDL2,HIDDEN
;
  s=size(f) & nx=s[1] & ny=s[2] & nz=s[3]
  if s[0] eq 2 then nz=1
  w=fltarr(nx,ny,nz,3)
;
  w[*,*,*,0]=xder(f)
  w[*,*,*,1]=yder(f)
  w[*,*,*,2]=zder(f)
;
  return,w
;
end
