;  
;  $Id$
;  for plotting brms and local growth rate
;  a summary of the results are in roberts.dat in this directory
;
!p.multi=[0,1,2]
lam=deriv(tt,alog(brms))
default,t1,18
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
good=where(tt gt t1)
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
