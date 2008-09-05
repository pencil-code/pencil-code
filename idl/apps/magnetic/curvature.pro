;$Id$
;
;  calculate B.gradB, where B is in terms of curl(A).
;
;  B_i = eps_ijk A_j,k      ;ie B_i=(curlA)_i
;
;  Next, (B.gradB)_i = B_l B_i,l
;                    = B_l eps_ijk A_j,kl
;
BB=curl(AA)
BdotgradB=fltarr(mx,my,mz,3)
;
for l=0,2 do begin
for k=0,2 do begin
for j=0,2 do begin
for i=0,2 do begin
  BdotgradB(*,*,*,i)=BdotgradB(*,*,*,i)+ $
    BB(*,*,*,l)*levi_civita(i,j,k)*derij(AA(*,*,*,j),k,l)
endfor
endfor
endfor
endfor
;
b2=dot2(bb(l1:l2,m1:m2,n1:n2,*))            
q=1./(b2 > 1e-4)
BdotgradBq=multv(q,BdotgradB(l1:l2,m1:m2,n1:n2,*))
;
tvscl,BdotgradBq(*,*,14,1)
;  .r ~/pencil-code/idl/apps/magnetic/curvature
END
