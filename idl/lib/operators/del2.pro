function del2,f,debug=debug
if keyword_set(debug) then print,'$Id: del2.pro,v 1.1 2004-05-23 15:42:07 brandenb Exp $'
;
;  laplacian in cartesian coordinates
;
s=size(f)
w=make_array(size=s,/nozero)
;
if s[0] eq 4 then begin
  w(*,*,*,0)=xder2(f(*,*,*,0))+yder2(f(*,*,*,0))+zder2(f(*,*,*,0))
  w(*,*,*,1)=yder2(f(*,*,*,1))+zder2(f(*,*,*,1))+xder2(f(*,*,*,1))
  w(*,*,*,2)=zder2(f(*,*,*,2))+xder2(f(*,*,*,2))+yder2(f(*,*,*,2))
endif else begin
  w=xder2(f)+yder2(f)+zder2(f)
endelse
;
return,w
end
