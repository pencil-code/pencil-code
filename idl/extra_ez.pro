;  $Id$
;
;  this is to be used in connection with r_ez.pro
;  it calculates a number of extra variables
;
gamma=5./3.
gamma_m1=gamma-1.
;
print,'calculate xx,yy,zz (comment out if there isnt enough memory!)'
xx = spread(x, [1,2], [my,mz])
yy = spread(y, [0,2], [mx,mz])
zz = spread(z, [0,1], [mx,my])
;
;  calculate extra stuff that may be of some convenience
;
xxx=x(l1:l2)
yyy=y(m1:m2)
zzz=z(n1:n2)
;
ff=f(l1:l2,m1:m2,n1:n2,*)
if (iuu ne 0) then uuu=f(l1:l2,m1:m2,n1:n2,iux-1:iuz-1)
if (ilnrho ne 0) then llnrho=f(l1:l2,m1:m2,n1:n2,ilnrho-1)
if (iss ne 0) then sss=f(l1:l2,m1:m2,n1:n2,iss-1)
if (iaa ne 0) then bb=curl(f(*,*,*,iax-1:iaz-1))
if (iaa ne 0) then jj=curl2(f(*,*,*,iax-1:iaz-1))

;ooo=oo(l1:l2,m1:m2,n1:n2,*)
if (iaa ne 0) then aaa=f(l1:l2,m1:m2,n1:n2,iax-1:iaz-1)
if (iaa ne 0) then bbb=bb(l1:l2,m1:m2,n1:n2,*)
if (iaa ne 0) then jjj=jj(l1:l2,m1:m2,n1:n2,*)
if (ilnrho ne 0) then rho=exp(llnrho)
if (iss ne 0) then cs2=exp(gamma_m1*llnrho+gamma*sss)
if (iss ne 0) then ppp=rho*cs2/gamma
;
;  calculate magnetic energy of mean field in the 3 directions
;
if (iaa ne 0) then begin
  bmx=sqrt(means(dot2_1d(means(means(bbb,3),2))))
  bmy=sqrt(means(dot2_1d(means(means(bbb,3),1))))
  bmz=sqrt(means(dot2_1d(means(means(bbb,2),1))))
  print,'bmx,bmy,bmz=',bmx,bmy,bmz
  print,'bmax=',sqrt(max(dot2(bbb)))
end
;
END
