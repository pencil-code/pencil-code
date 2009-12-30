;
;  Plot profiles of some variables for shock tube problem
;
;window,xs=800,ys=600
save_state
!p.charsize=1.0
!p.multi=[0,2,2]
!x.title='!8y!X'

pc_read_param,obj=par 
pc_read_param,obj=par2,/param2
pc_read_var,obj=data,/trimall,variables=['pp','cs2'],/magic,/add
gamma=par.gamma
gamma1=par.gamma-1
ss_left=par.ss_left & ss_right=par.ss_right
rho_left=par.rho_left & rho_right=par.rho_right
p_left=exp(gamma*(ss_left+alog(rho_left)))/gamma
p_right=exp(gamma*(ss_right+alog(rho_right)))/gamma
;
;  exclude ghost zones
;
yyy=data.y
uuu=data.uu[*,1]
sss=data.ss
if (par.ldensity_nolog) then begin
  rho=data.lnrho
  llnrho=alog(data.lnrho)
endif else begin
  rho=exp(data.lnrho)
  llnrho=data.lnrho
endelse
cs2=data.cs2
ppp=data.pp

;
;  get exact solution
;
shocktube,yyy,data.t,p,rho,u,[0.,p_left,rho_left],[0.,p_right,rho_right],gamma
circ_sym,.6,0
ps=8
thi=1
li=0
;
;  velocity
;
plot, yyy, u, /NODATA, $
    TITLE='!6Velocity!X',YTITLE='!8u!X', PSYM=ps, $
    YRANGE=stretchrange(minmax([minmax(u),minmax(uuu)]),.1),back=255,col=1

oplot,yyy,uuu,ps=PS,col=1
oplot,yyy,u,col=122,li=li,thi=thi
;
;  pressure
;
plot_io,yyy,ppp,TITLE='!6Pressure!X',PSYM=ps,col=1, $
    YRANGE=stretchrange(minmax(ppp),.1,/log)
oplot,yyy,p,col=122,li=li,thi=thi
;
;  density
;
plot_io,yyy,exp(llnrho),yst=3,TITLE='!6Density!X',PSYM=ps,col=1, $
    YRANGE=stretchrange(minmax(rho),.1,/log)
oplot,yyy,rho,col=122,li=li,thi=thi
;
;  entropy
;
s=alog(gamma*p)/gamma-alog(rho)
sss=alog(gamma*ppp)/gamma-llnrho
plot,yyy,sss,yst=3,TITLE='!6Entropy!X',PSYM=ps,col=1, $
    YRANGE=stretchrange(minmax(s),.1)
oplot,yyy,s,col=122,li=li,thi=thi

restore_state

;
;  print momentum and theoretical value
;
py = total(exp(llnrho)*uuu[*,*,*,0])*data.dy
print, 'momentum p_y = ', py, '; theor.:', (p_left-p_right)*data.t
print,'nu/(dy*cs)=',par2.nu/(data.dy*sqrt(par.gamma))

;print,'import sod_10.gif'
;print,'import sod_100.gif'
;print,'import sod_1000.gif'
END
