;$Id: pc_fluct_stat.pro,v 1.2 2007-10-20 17:57:01 brandenb Exp $
pro pc_fluct_stat,t,f,fm,frms,frms_error,plo=plo,good=good,debug=debug
;
;  calculate average and rms fluctuations
;
fm=mean(f(good))
frms=sqrt(mean((f(good)-fm)^2))
if keyword_set(debug) then print,fm,frms
;
;  plot the accumulated mean
;
if keyword_set(plo) then begin
  fac=2.
  yr=minmax([fac*fm,-fac*fm,fac*frms])
  ;plot,t,accum(f),li=1,yr=yr
  plot,t,f,li=1,yr=yr
  pc_error_range,t(good),f(good),mean=fm,error=fm_error,/accum,/oplot
  ;
  ;  g = (f-fm)^2
  ;  g = gm +/- gm_error = gm*(1 +/- gm_error/gm)
  ;  frms = sqrt(gm)*(1 +/- .5*gm_error/gm)
  ;  frms = sqrt(gm) +/- .5*gm_error/sqrt(gm)
  ;  frms = sqrt(gm) +/- .5*gm_error/frms
  ;
  g=(f(good)-fm)^2
  oplot,t,accum(g)^.5
  pc_error_range,t(good),g,mean=gm,error=gm_error,/accum,/oplot
  frms=gm^.5
  frms_error=.5*gm_error/frms
  oplot,t,t*0.,li=3
endif
;
END
