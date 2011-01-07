;$Id$
;
;  Routine for computing the growth rate.
;  Use data in the range from  t1 to t2
;
t1=40 & t2=80
;
;  Read data.
;
pc_read_param,o=param,/param2
pc_read_ts,o=ts
tt=ts.t
brms=ts.brms
;
;  Isolate the range within t1 and t2.
;
good=where(tt gt t1 and tt lt t2)
lam=deriv(tt,alog(deriv(brms)))
;
;  Plot the actual growth of brms.
;
!p.multi=[0,1,2]
plot_io,tt,brms
;
;  Plot the growth rate.
;
plot,tt,lam,yr=[0,.5]
lamm=mean(lam(good))
oplot,tt(good),tt(good)*0.+lamm,col=122
;
;  Print the result.
;
fo='(2f8.3)'
print,param.tau1_emf,lamm,fo=fo
END
