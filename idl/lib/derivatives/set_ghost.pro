;***********************************************************************
pro set_ghost,f,debug=debug
;
;  set ghost zones
;  currently, only periodic boundary conditions are prepared
;
if keyword_set(debug) then print,'$Id: set_ghost.pro,v 1.1 2004-05-23 15:46:10 brandenb Exp $'
;
help,f
s=size(f) & mx=s[1] & my=s[2] & mz=s[3]
nghost=3
print,'mx,my,mz=',mx,my,mz
;
f(0:nghost-1,*,*,*)=f(mx-2*nghost:mx-nghost-1,*,*,*)
f(*,0:nghost-1,*,*)=f(*,my-2*nghost:my-nghost-1,*,*)
f(*,*,0:nghost-1,*)=f(*,*,mz-2*nghost:mz-nghost-1,*)
;
f(mx-nghost:mx-1,*,*,*)=f(nghost:2*nghost-1,*,*,*)
f(*,my-nghost:my-1,*,*)=f(*,nghost:2*nghost-1,*,*)
f(*,*,mz-nghost:mz-1,*)=f(*,*,nghost:2*nghost-1,*)
;
end
