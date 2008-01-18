;
;  $Id: der_single.pro,v 1.1 2008-01-18 10:39:04 ajohan Exp $
;
;  Calculate single component of firstst derivative; the second argument
;  gives the direction.
;
;  20-mar-04/axel: coded
;
function der_single, f, j
;
if (j eq 0) then return, xder(f)
if (j eq 1) then return, yder(f)
if (j eq 2) then return, zder(f)
;
end
