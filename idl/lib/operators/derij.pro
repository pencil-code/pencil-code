function derij,f,i,j,debug=debug
if keyword_set(debug) then print,'$Id: derij.pro,v 1.1 2004-03-20 10:24:54 brandenb Exp $'
;
;  calculate derivative matrix
;  20-mar-04/axel: coded
;
if i eq 0 and j eq 0 then return,xder2(f)
if i eq 0 and j eq 1 then return,xder(yder(f))
if i eq 0 and j eq 2 then return,xder(zder(f))
if i eq 1 and j eq 0 then return,yder(xder(f))
if i eq 1 and j eq 1 then return,yder2(f)
if i eq 1 and j eq 2 then return,yder(zder(f))
if i eq 2 and j eq 0 then return,zder(xder(f))
if i eq 2 and j eq 1 then return,zder(yder(f))
if i eq 2 and j eq 2 then return,zder2(f)
;
end
