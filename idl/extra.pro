print,'calculate xx,yy,zz (comment out if there isnt enough memory!)'
xx = spread(x, [1,2], [my,mz])
yy = spread(y, [0,2], [mx,mz])
zz = spread(z, [0,1], [mx,my])
;
;  calculate extra stuff that may be of some convenience
;
oo=curl(uu)
bb=curl(aa)
jj=curl2(aa)
xxx=x(l1:l2)
xyy=y(m1:m2)
xzz=z(n1:n2)
uuu=uu(l1:l2,m1:m2,n1:n2,*)
ooo=oo(l1:l2,m1:m2,n1:n2,*)
aaa=aa(l1:l2,m1:m2,n1:n2,*)
bbb=bb(l1:l2,m1:m2,n1:n2,*)
jjj=jj(l1:l2,m1:m2,n1:n2,*)
;sss=ss(l1:l2,m1:m2,n1:n2)
llnrho=lnrho(l1:l2,m1:m2,n1:n2)
;
;  calculate magnetic energy of mean field in the 3 directions
;
EMx=mean(means(means(bbb,3),2)^2)
EMy=mean(means(means(bbb,3),1)^2)
EMz=mean(means(means(bbb,2),1)^2)
print,'Emx,EMy,EMz=',Emx,EMy,EMz
;
END
