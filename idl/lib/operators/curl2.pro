function curl2,f,debug=debug
if keyword_set(debug) then print,'$Id: curl2.pro,v 1.1 2002-07-26 11:57:31 brandenb Exp $'
;
;  double curl in cartesian coordinates
;
s=size(f)
w=make_array(size=s,/nozero)
;
  w(*,*,*,0)=xder(yder(f(*,*,*,1))+zder(f(*,*,*,2)))-yder2(f(*,*,*,0))-zder2(f(*,*,*,0))
  w(*,*,*,1)=yder(zder(f(*,*,*,2))+xder(f(*,*,*,0)))-zder2(f(*,*,*,1))-xder2(f(*,*,*,1))
  w(*,*,*,2)=zder(xder(f(*,*,*,0))+yder(f(*,*,*,1)))-xder2(f(*,*,*,2))-yder2(f(*,*,*,2))
;
return,w
end
