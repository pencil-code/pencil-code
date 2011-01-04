;  
;  $Id$
;  for plotting brms and local growth rate
;  a summary of the results are in roberts.dat in this directory
;
t1=.4 & t2=2
t1=2. & t2=4
!p.multi=[0,1,2]
pc_read_ts,o=ts
tt=ts.t
brms=ts.brms
tmax=max(tt)
lam=deriv(tt,alog(brms))
default,t1,(1<(tmax/2.))
default,t2,(4<(tmax/2.))
default,yr,[-1,1]
@param
;
;  brms
;
plot_io,tt,brms
;
;  growth rate
;
plot,tt,lam,yr=yr
;
good=where(tt gt t1 and tt lt t2)
if good(0) ne -1 then begin
  oplot,tt(good),lam(good),col=122
  lamm=mean(lam(good))
  oplot,tt(good),tt(good)-tt(good)+lamm,li=3,col=188
  print,'lam=',lamm
endif else begin
  print,'time series still too short for getting growth rate'
endelse
;
END
