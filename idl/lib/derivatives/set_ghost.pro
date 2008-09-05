;***********************************************************************
pro set_ghost,f,debug=debug,add=add,nghost=nghost
;
;  set ghost zones
;  currently, only periodic boundary conditions are prepared
;  This procedure is probably obsolete now that I've checked in pc_setghost.
;  I keep this in case it proves advantageous with big data sets.
;

default,nghost,3
if keyword_set(add) then begin
  forig=f
  sorig=size(f) & mx=sorig[1]+2*nghost & my=sorig[2]+2*nghost & mz=sorig[3]+2*nghost 

  
  mvar=1
  if sorig[0] eq 4 then mvar=sorig[4]

  f=fltarr(mx,my,mz,mvar)*forig[0]

  f(nghost:mx-nghost-1,nghost:my-nghost-1,nghost:mz-nghost-1,*)=forig(*,*,*,*)
endif
forig=0.

s=size(f) & mx=s[1] & my=s[2] & mz=s[3]
;
;  debug output
;
if keyword_set(debug) then begin
  print,'$Id$'
  help,f
  print,'mx,my,mz=',mx,my,mz
endif
;
;  set ghost zones on the left end
;
f(0:nghost-1,*,*,*)=f(mx-2*nghost:mx-nghost-1,*,*,*)
f(*,0:nghost-1,*,*)=f(*,my-2*nghost:my-nghost-1,*,*)
f(*,*,0:nghost-1,*)=f(*,*,mz-2*nghost:mz-nghost-1,*)
;
;  set ghost zones on the right end
;
f(mx-nghost:mx-1,*,*,*)=f(nghost:2*nghost-1,*,*,*)
f(*,my-nghost:my-1,*,*)=f(*,nghost:2*nghost-1,*,*)
f(*,*,mz-nghost:mz-1,*)=f(*,*,nghost:2*nghost-1,*)
;
end
