;  $Id: extra.pro,v 1.17 2003-03-10 15:19:29 brandenb Exp $
;
;  This routine calculates a number of extra variables
;
gamma=5./3.
gamma1=gamma-1.
;
;print,'calculate xx,yy,zz (comment out if there isnt enough memory!)'
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
;
if (ilnrho ne 0) then llnrho=lnrho(l1:l2,m1:m2,n1:n2)
if (ilnrho ne 0) then rho=exp(llnrho)
if (ilnrho ne 0 and ient eq 0) then cs2=cs20*exp(gamma1*llnrho)
if (ient ne 0) then sss=ss(l1:l2,m1:m2,n1:n2)
if (ient ne 0) then cs2=cs20*exp(gamma1*llnrho+gamma*sss)
if (ient ne 0) then ppp=rho*cs2/gamma
;
if (iuu ne 0) then for j=0,2 do print,sqrt(mean(uu(*,*,*,j)^2))
if (iuu ne 0) then print
;
;  calculate magnetic energy of mean field in the 3 directions
;
if (iaa ne 0) then begin
  for j=0,2 do bbb(*,*,*,j)=bbb(*,*,*,j)+par2.b_ext(j)
  bmx=sqrt(mean(dot2_1d(means(means(bbb,3),2))))
  bmy=sqrt(mean(dot2_1d(means(means(bbb,3),1))))
  bmz=sqrt(mean(dot2_1d(means(means(bbb,2),1))))
  print,'bmx,bmy,bmz=',bmx,bmy,bmz
end
;
;  calculate vertical averages
;
if (ient ne 0) then begin
  if (nz gt 1) then begin
    cs2m=haver(cs2) & csm=sqrt(cs2m)
    rhom=haver(rho)
  end
end
;
;  in case we need it
;
hhh=cs2/gamma1  ;(enthalpy)
;
END
