function cross,g,f,debug=debug
if keyword_set(debug) then print,'$Id: cross.pro,v 1.1 2002-07-28 17:12:56 brandenb Exp $'
s=size(f) & nx=s(1) & ny=s(2) & nz=s(3); & print,nx,ny,nz
w=fltarr(nx,ny,nz,3)
w(*,*,*,0)=g(*,*,*,1)*f(*,*,*,2)-g(*,*,*,2)*f(*,*,*,1)
w(*,*,*,1)=g(*,*,*,2)*f(*,*,*,0)-g(*,*,*,0)*f(*,*,*,2)
w(*,*,*,2)=g(*,*,*,0)*f(*,*,*,1)-g(*,*,*,1)*f(*,*,*,0)
return,w
end
