gamma=5./3.
gamma1=gamma-1.
;
print,'calculate xx,yy,zz (comment out if there isnt enough memory!)'
xx = spread(x, [1,2], [my,mz])
yy = spread(y, [0,2], [mx,mz])
zz = spread(z, [0,1], [mx,my])
;
;  calculate extra stuff that may be of some convenience
;
oo=curl(uu)
if (iaa ne 0) then bb=curl(aa)
if (iaa ne 0) then jj=curl2(aa)
xxx=x(l1:l2)
yyy=y(m1:m2)
zzz=z(n1:n2)
uuu=uu(l1:l2,m1:m2,n1:n2,*)
ooo=oo(l1:l2,m1:m2,n1:n2,*)
if (iaa ne 0) then aaa=aa(l1:l2,m1:m2,n1:n2,*)
if (iaa ne 0) then bbb=bb(l1:l2,m1:m2,n1:n2,*)
if (iaa ne 0) then jjj=jj(l1:l2,m1:m2,n1:n2,*)
if (ilnrho ne 0) then llnrho=lnrho(l1:l2,m1:m2,n1:n2)
if (ient ne 0) then sss=ss(l1:l2,m1:m2,n1:n2)
if (ient ne 0) then cs2=exp(gamma1*llnrho+gamma*sss)
;
;  calculate magnetic energy of mean field in the 3 directions
;
if (iaa ne 0) then begin
  EMx=mean(means(means(bbb,3),2)^2)
  EMy=mean(means(means(bbb,3),1)^2)
  EMz=mean(means(means(bbb,2),1)^2)
  print,'Emx,EMy,EMz=',Emx,EMy,EMz
end
;
END
