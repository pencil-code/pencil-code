function unitvector,f,debug=debug
if keyword_set(debug) then print,'$Id$'
;
;  calculate unit vector
;  20-mar-04/axel: coded
;
w=make_array(size=size(f),/nozero)
;
tiny=1e-30
f1=1./(sqrt(dot2(f)) > tiny)
;
s=size(f1)
;
  case (s[0]) of
    0: for j=0,2 do w(j)=f1*f(j)
    1: for j=0,2 do w(*,j)=f1*f(*,j)
    2: for j=0,2 do w(*,*,j)=f1*f(*,*,j)
    3: for j=0,2 do w(*,*,*,j)=f1*f(*,*,*,j)
    else: message, "dot: can't handle fields with " + strtrim (s, 2) + " dimensions."
  endcase
;
return,w
end
