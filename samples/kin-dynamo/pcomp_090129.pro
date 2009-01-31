;$Id$
; 
;  plot to demonstrate that the growth rates are the same
;  before and after the pencilization of the initial condition.
;
pc_read_ts,file='./reference1000_before090129.out',o=ts0
pc_read_ts,file='./reference1000_after090129.out',o=ts1
;
!p.charsize=2
!p.multi=[0,1,2]
plot_io,ts0.t,ts0.brms,yr=[4e-6,1e-3],xtit='t',ytit='Brms'
oplot,ts1.t,ts1.brms,li=2
;
plot,ts0.t,deriv(ts0.t,alog(ts0.brms)),yr=[.09,.11],xtit='t',ytit='lambda'
oplot,ts1.t,deriv(ts1.t,alog(ts1.brms)),li=2,col=122
;
END
