function der,f,j,debug=debug
if keyword_set(debug) then print,'$Id: der.pro,v 1.1 2004-03-20 10:24:54 brandenb Exp $'
;
;  calculate 1st derivative; the second argument gives the direction
;  20-mar-04/axel: coded
;
if j eq 0 then return,xder(f)
if j eq 1 then return,yder(f)
if j eq 2 then return,zder(f)
;
end
