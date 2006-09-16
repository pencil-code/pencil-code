function cross,g,f,debug=debug
Id='$Id: cross.pro,v 1.3 2006-09-16 11:31:50 brandenb Exp $'
if keyword_set(debug) then print,Id
w=make_array(size=size(f))
w(*,*,*,0)=g(*,*,*,1)*f(*,*,*,2)-g(*,*,*,2)*f(*,*,*,1)
w(*,*,*,1)=g(*,*,*,2)*f(*,*,*,0)-g(*,*,*,0)*f(*,*,*,2)
w(*,*,*,2)=g(*,*,*,0)*f(*,*,*,1)-g(*,*,*,1)*f(*,*,*,0)
return,w
end
