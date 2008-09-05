function curlij,f,debug=debug
if keyword_set(debug) then print,'$Id$'
;
;  derivative matrix from curl in cartesian coordinates
;  used to calculate B_{i,j} as bij=curlij(aa)
;  20-mar-04/axel: adapted from f90 routine bij_etc
;
s=size(f) & mx=s(1) & my=s(2) & mz=s(3)
bij=fltarr(mx,my,mz,3,3)
;
;  calculate B_i,j = eps_ikl A_l,jk
;
for i=0,2 do begin
for j=0,2 do begin
for k=0,2 do begin
for l=0,2 do begin
  eps=levi_civita(i,k,l)
  if eps ne 0. then bij(*,*,*,i,j)=bij(*,*,*,i,j)+eps*derij(f(*,*,*,l),j,k)
endfor
endfor
endfor
endfor
;
return,bij
end
