!p.multi=[0,1,3]
!p.charsize=3
!x.title='Time [!4l!6s]'
@data/pc_constants.pro
!x.range=[0,100]


time=ts.t*1e6
pressure_pascal=ts.PPm/10
;
; Plot temperature and pressure
;
plot,time,ts.TTm
oplot,time,pressure_pascal/100,li=2
xx=25
dx=8
yy=2000
dy=-300
legend,xx,dx,yy+dy*0,0,'Temperature [K]'
legend,xx,dx,yy+dy*1,2,'Pressure [Pa/100]'
;
; Show species mass fractions
;
plot,time,ts.Y1m,yrange=[0.0,0.25],ytit='!8Y!6 [kg/kg]'
oplot,time,ts.Y2m,li=2
oplot,time,ts.Y3m,li=3
xx=25
dx=8
yy=0.17
dy=-0.018
legend,xx,dx,yy+dy*0,0,specname(0)
legend,xx,dx,yy+dy*1,2,specname(1)
legend,xx,dx,yy+dy*2,3,specname(2)
;
; Show species mass fractions (different yrange)
;
plot,time,ts.Y1m,yrange=[0.0,0.03],ytit='!8Y!6 [kg/kg]'
;oplot,time,ts.Y2m,li=2
;oplot,time,ts.Y3m,li=3
oplot,time,ts.Y4m,li=4
oplot,time,ts.Y5m,li=5
oplot,time,ts.Y6m,li=6
oplot,time,ts.Y7m,li=7
xx=25
dx=8
yy=0.02
dy=-0.0018
legend,xx,dx,yy+dy*0,0,specname(0)
;legend,xx,dx,yy+dy*1,2,specname(1)
;legend,xx,dx,yy+dy*2,3,specname(2)
legend,xx,dx,yy+dy*1,4,specname(3)
legend,xx,dx,yy+dy*2,5,specname(4)
legend,xx,dx,yy+dy*3,6,specname(5)
legend,xx,dx,yy+dy*4,7,specname(6)


END
