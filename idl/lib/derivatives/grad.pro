;;; Calculate the gradient of a 3-d scalar field
function grad,f,coord_system
COMPILE_OPT IDL2,HIDDEN
;
default,coord_system,'cartesian'
;
  s=size(f) & nx=s[1] & ny=s[2] & nz=s[3]
  if s[0] eq 2 then nz=1
  w=fltarr(nx,ny,nz,3)
;
  lsystem=-1
  if (coord_system eq 'cartesian') then lsystem=0
  if (coord_system eq 'cylindric') then lsystem=1
  if (coord_system eq 'spherical') then lsystem=2
  if (lsystem eq -1) then $
    print,'coord_system= ',coord_system,' is not valid'
;
  w[*,*,*,0]=xder(f)
  w[*,*,*,1]=yder(f,lsystem)
  w[*,*,*,2]=zder(f,lsystem)
;
  return,w
;
end
