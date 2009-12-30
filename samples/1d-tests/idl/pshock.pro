;
;  Plot profiles of some variables for shock tube problem
;
;window,xs=800,ys=600
save_state
!p.charsize=1.0
!p.multi=[0,2,2]
!x.title='!8x!X'

pc_read_param,obj=par 
pc_read_param,obj=par2,/param2
pc_read_var,obj=data,/trimall,variables=['pp','cs2'],/magic,/add
gamma=par.gamma
gamma1=par.gamma-1
ss_left=par.ss_left & ss_right=par.ss_right
rho_left=par.rho_left(0) & rho_right=par.rho_right(0)
p_left=exp(gamma*(ss_left+alog(rho_left)))/gamma
p_right=exp(gamma*(ss_right+alog(rho_right)))/gamma
;
;  exclude ghost zones
;
xxx=data.x
uuu=data.uu
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
shocktube,xxx,data.t,p,rho,u,[0.,p_left,rho_left],[0.,p_right,rho_right],gamma
circ_sym,.6,0
ps=8
thi=1
li=0
;
;  velocity
;
plot, xxx, u, /NODATA, $
    TITLE='!6Velocity!X',YTITLE='!8u!X', PSYM=ps, $
    YRANGE=stretchrange(minmax([minmax(u),minmax(uuu)]),.1),back=255,col=1

oplot,xxx,uuu,ps=PS,col=1
oplot,xxx,u,col=122,li=li,thi=thi
;
;  pressure
;
plot_io,xxx,ppp,TITLE='!6Pressure!X',PSYM=ps,col=1, $
    YRANGE=stretchrange(minmax(p),.1,/log)
oplot,xxx,p,col=122,li=li,thi=thi
;
;  density
;
plot_io,xxx,exp(llnrho),yst=3,TITLE='!6Density!X',PSYM=ps,col=1, $
    YRANGE=stretchrange(minmax(rho),.1,/log)
oplot,xxx,rho,col=122,li=li,thi=thi
;
;  entropy
;
s=alog(gamma*p)/gamma-alog(rho)
sss=alog(gamma*ppp)/gamma-llnrho
plot,xxx,sss,yst=3,TITLE='!6Entropy!X',PSYM=ps,col=1, $
    YRANGE=stretchrange(minmax(s),.1)
oplot,xxx,s,col=122,li=li,thi=thi

restore_state

;
;  print momentum and theoretical value
;
px = total(exp(llnrho)*uuu[*,*,*,0])*data.dx
print, 'momentum p_x = ', px, '; theor.:', (p_left-p_right)*data.t
print,'nu/(dx*cs)=',par2.nu/(data.dx*sqrt(gamma))

print, 'RMS Error - Density  = ', RMS(exp(llnrho)-rho)
print, 'RMS Error - Entropy  = ', RMS(sss-s)
print, 'RMS Error - Velocity = ', RMS(uuu-u)
print, 'RMS Error - Pressure = ', RMS(ppp-p)

print,'import sod_10.gif'
print,'import sod_100.gif'
print,'import sod_1000.gif'
END
