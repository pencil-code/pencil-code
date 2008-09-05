function der_unitvector,f,fij,i,j,debug=debug
if keyword_set(debug) then print,'$Id$'
;
;  calculate 1st derivative of a unitvector
;  f is vector, and fij is its gradient matrix (f is *not* a unitvector)
;  fhat_{i,j} = (delta_ik - fhati*fhatk) f_{k,j}
;  20-mar-04/axel: coded
;
tiny=1e-30
f1=1./(sqrt(dot2(f)) > tiny)
;
;  fhat_{i,j} = (1/f)*[f_{i,j} - fhati*fhatk*f_{k,j}]
;  
fhat=unitvector(f)
fhatij=f1*fij(*,*,*,i,j)
for k=0,2 do fhatij=fhatij-f1*fhat(*,*,*,i)*fhat(*,*,*,k)*fij(*,*,*,k,j)
;
return,fhatij
end
