;$Id: pc_fluct_stat.pro,v 1.1 2007-09-29 10:38:15 brandenb Exp $
pro pc_fluct_stat,t,f,fm,frms,frms_error,plo=plo,good=good
;
;  calculate average and rms fluctuations
;
fm=mean(f(good))
frms=sqrt(mean((f(good)-fm)^2))
;
;  plot the accumulated mean
;
if keyword_set(plo) then begin
  fac=2.
  yr=minmax([fac*fm,-fac*fm,fac*frms])
  plot,t,accum(f),li=1,yr=yr
  pc_error_range,t(good),f(good),mean=fm,error=fm_error,/accum,/oplot
  ;
  g=sqrt((f-fm)^2)
  oplot,t,accum(g)
  pc_error_range,t(good),g(good),mean=frms,error=frms_error,/accum,/oplot
  oplot,t,t*0.,li=3
endif
;
END
