
;
;  Plot profiles of some variables for shock tube problem
;
;window,xs=800,ys=600
save_state
!p.charsize=1.0
!p.multi=[0,2,2]
!x.title='!8z!X'

ss_left=par.ss_left & ss_right=par.ss_right
rho_left=par.rho_left & rho_right=par.rho_right
p_left=exp(gamma*(ss_left+alog(rho_left)))/gamma
p_right=exp(gamma*(ss_right+alog(rho_right)))/gamma
;
;  exclude ghost zones
;
zzz=z(n1:n2)
uuu=uu(l1:l2,m1:m2,n1:n2,2)
sss=ss(l1:l2,m1:m2,n1:n2)
llnrho=lnrho(l1:l2,m1:m2,n1:n2)
cs2=cs20*exp(gamma1*llnrho+gamma*sss)
rho=exp(llnrho)
ppp=rho*cs2/gamma
;
;  get exact solution
;
shocktube,zzz,t,p,rho,u,[0.,p_left,rho_left],[0.,p_right,rho_right],gamma
circ_sym,.6,0
ps=8
thi=1
li=0
;
;  velocity
;
plot, zzz, u, /NODATA, $
    TITLE='!6Velocity!X',YTITLE='!8u!X', PSYM=ps, $
    YRANGE=stretchrange(minmax([minmax(u),minmax(uuu)]),.1),back=255,col=1

oplot,zzz,uuu,ps=PS,col=1
oplot,zzz,u,col=122,li=li,thi=thi
;
;  pressure
;
plot_io,zzz,ppp,TITLE='!6Pressure!X',PSYM=ps,col=1, $
    YRANGE=stretchrange(minmax(p),.1,/log)
oplot,zzz,p,col=122,li=li,thi=thi
;
;  density
;
plot_io,zzz,exp(llnrho),yst=3,TITLE='!6Density!X',PSYM=ps,col=1, $
    YRANGE=stretchrange(minmax(rho),.1,/log)
oplot,zzz,rho,col=122,li=li,thi=thi
;
;  entropy
;
s=alog(gamma*p)/gamma-alog(rho)
sss=alog(gamma*ppp)/gamma-llnrho
plot,zzz,sss,yst=3,TITLE='!6Entropy!X',PSYM=ps,col=1, $
    YRANGE=stretchrange(minmax(s),.1)
oplot,zzz,s,col=122,li=li,thi=thi

restore_state

;
;  print momentum and theoretical value
;
pz = total(exp(llnrho)*uuu[*,*,*,0])*dz
print, 'momentum p_z = ', pz, '; theor.:', (p_left-p_right)*t
print,'nu/(dz*cs)=',nu/(dz*sqrt(gamma))

;print,'import sod_10.gif'
;print,'import sod_100.gif'
;print,'import sod_1000.gif'
END
