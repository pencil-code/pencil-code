function haver,f,xvertical=xvertical,yvertical=yvertical
;
;  Horizontal average of 3-d scalar
;  $Id: haver.pro,v 1.1 2002-07-22 22:59:40 brandenb Exp $
;
	s=size(f)
	if n_elements(xvertical) ne 0 then begin
	  h=fltarr(s(1))
	  for n=0,s(1)-1 do h(n)=aver(f(n,*,*))
	end else if n_elements(yvertical) ne 0 then begin
	  h=fltarr(s(2))
	  for n=0,s(2)-1 do h(n)=aver(f(*,n,*))
	end else begin
	  h=fltarr(s(3))
	  for n=0,s(3)-1 do h(n)=aver(f(*,*,n)) 
	end
	return,h
end
