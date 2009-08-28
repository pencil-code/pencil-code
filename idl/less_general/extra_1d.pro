print,'calculate xx,yy,zz (comment out if there isnt enough memory!)'
xx = spread(x, [1,2], [my,mz])
yy = spread(y, [0,2], [mx,mz])
zz = spread(z, [0,1], [mx,my])
;
;  calculate extra stuff that may be of some convenience
;
xxx=x(l1:l2)
xyy=y(m1:m2)
xzz=z(n1:n2)
uuu=uu(l1:l2,m1:m2,n1:n2,*)
sss=ss(l1:l2,m1:m2,n1:n2)
llnrho=lnrho(l1:l2,m1:m2,n1:n2)
;
cs2=exp(gamma_m1*llnrho+gamma*ss)
eee=cs2/(gamma_m1*gamma)
rho=exp(llnrho)
;
END
