function unitvector,f,debug=debug
if keyword_set(debug) then print,'$Id: unitvector.pro,v 1.1 2004-03-20 10:24:54 brandenb Exp $'
;
;  calculate unit vector
;  20-mar-04/axel: coded
;
w=make_array(size=size(f),/nozero)
;
tiny=1e-30
f1=1./(sqrt(dot2(f)) > tiny)
for j=0,2 do w(*,*,*,j)=f1*f(*,*,*,j)
;
return,w
end
