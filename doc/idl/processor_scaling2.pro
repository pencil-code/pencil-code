N=2^findgen(5)*128
P=2^findgen(10)
nN=n_elements(N)
nP=n_elements(P)
W=fltarr(nP,nN)
for i=0,nP-1 do begin
for j=0,nN-1 do begin
  W(i,j)=(N(j)+6.)*(N(j)/sqrt(P(i))+6.)^2*P(i)/N(j)^3
endfor
endfor
;
fo="(i5,8f9.2)" & fo2="(5x,8i9)"
fo="(i5,8('  &',f7.2))" & fo2="(5x,8i9)"
print,reform(fix(N(*))),fo=fo2
print
for i=0,nP-1 do print,fix(P(i)),reform(W(i,*)),fo=fo
END
