;***********************************************************************
function pc_setghost,f,debug=debug
;
;  set ghost zones
;  currently, only periodic boundary conditions are prepared
;
s=size(f) & mx=s[1] & my=s[2] & mz=s[3]
nghost=3
;
;  debug output
;
if keyword_set(debug) then begin
  print,'$Id: pc_setghost.pro,v 1.1 2004-07-04 11:10:05 brandenb Exp $'
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
return,f
end
