;
;  $Id: gijcurl.pro,v 1.3 2008-04-30 15:00:10 ajohan Exp $
;
;  Derivative matrix from curl in cartesian coordinates.
;  Used to calculate B_{i,j} as bij=gijcurl(aa).
;
;  20-mar-04/axel: adapted from f90 routine bij_etc
;
function gijcurl, f
COMPILE_OPT IDL2,HIDDEN
common cdat_coords, coord_system
;
if (coord_system ne 'cartesian') then message, $
    "gijcurl not yet implemented for coord_system='" + coord_system + "'"
;
s=size(f) & mx=s[1] & my=s[2] & mz=s[3]
bij=fltarr(mx,my,mz,3,3)
;
;  Calculate B_i,j = eps_ikl A_l,jk
;
for i=0,2 do begin & for j=0,2 do begin & for k=0,2 do begin &for l=0,2 do begin
  eps=levi_civita(i,k,l)
  if (eps ne 0.0) then $
      bij[*,*,*,i,j]=bij[*,*,*,i,j]+eps*derij_single(f[*,*,*,l],j,k)
endfor & endfor & endfor & endfor
;
return, bij
;
end
