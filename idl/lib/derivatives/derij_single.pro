;
;  $Id: derij_single.pro,v 1.1 2008-01-18 10:39:04 ajohan Exp $
;
;  Calculate single component of second derivative matrix.
;
;  20-mar-04/axel: coded
;
function derij_single, f, i, j
;
if (i eq 0 and j eq 0) then return, xder2(f)
if (i eq 0 and j eq 1) then return, xderyder(f)
if (i eq 0 and j eq 2) then return, xderzder(f)
if (i eq 1 and j eq 0) then return, yderxder(f)
if (i eq 1 and j eq 1) then return, yder2(f)
if (i eq 1 and j eq 2) then return, yderzder(f)
if (i eq 2 and j eq 0) then return, zderxder(f)
if (i eq 2 and j eq 1) then return, zderyder(f)
if (i eq 2 and j eq 2) then return, zder2(f)
;
end
