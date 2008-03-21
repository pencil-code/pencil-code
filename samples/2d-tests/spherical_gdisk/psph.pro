
pc_read_ts,obj=ts
pc_read_param,obj=par
!p.charsize=1.5

t=ts.t

rp=ts.xpar1  &  phip=ts.zpar1 & thtp=ts.ypar1
rs=ts.xpar2  &  phis=ts.zpar2 & thts=ts.ypar2

xp=rp*sin(thtp)*cos(phip)
yp=rp*sin(thtp)*sin(phip)

xs=rs*sin(thts)*cos(phis)
ys=rs*sin(thts)*sin(phis)

planetm=par.pmass(0)
starm=1-par.pmass(0)
totmass=1.

;center of mass

xcm = (starm*xs + planetm*xp)/totmass
ycm = (starm*ys + planetm*yp)/totmass

!p.multi=[0,2,2]


loadct,3
nt=n_elements(ts.t)-1
aa=fltarr(1)
bb=fltarr(1)
theta=grange(-!pi,!pi,100)


plot,xs,ys  ,title='star position',xs=3,ys=3,xr=[-1.1*planetm,1.1*planetm],$
  yr=[-1.1*planetm,1.1*planetm]
  aa(0)=xs(nt)&bb(0)=ys(nt)&oplot,aa,bb,color=200,ps=2


plot,xp,yp  ,title='planet position',xs=3,ys=3,xr=[-1.1*starm,1.1*starm],$
  yr=[-1.1*starm,1.1*starm]
  aa(0)=xp(nt)&bb(0)=yp(nt)&oplot,aa,bb,color=200,ps=2

plot,xcm,ycm,ps=3,title='cm position',xs=3,ys=3
aa(0)=xcm(nt)&bb(0)=ycm(nt)&oplot,aa,bb,color=200,ps=2



end
