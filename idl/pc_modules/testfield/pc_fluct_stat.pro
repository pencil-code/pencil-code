;$Id$
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
  iplo=1
endif else begin
  iplo=0
endelse
;
;  plot the accumulated mean
;
if iplo ne 0 then begin
  fac=2.
  yr=minmax([fac*fm,-fac*fm,fac*frms])
  plot,t,f,li=1,yr=yr
endif
pc_error_range2,t(good),f(good),mean=fm,error=fm_error,/accum,oplot=iplo
;
;  g = (f-fm)^2
;  g = gm +/- gm_error = gm*(1 +/- gm_error/gm)
;  frms = sqrt(gm)*(1 +/- .5*gm_error/gm)
;  frms = sqrt(gm) +/- .5*gm_error/sqrt(gm)
;  frms = sqrt(gm) +/- .5*gm_error/frms
;
g=(f(good)-fm)^2
if iplo ne 0 then begin
  plot,t,f,li=1,yr=yr
  oplot,t,accum(g)^.5
  oplot,t,t*0.,li=3
endif
;
pc_error_range2,t(good),g,mean=gm,error=gm_error,/accum,oplot=iplo
frms=gm^.5
frms_error=.5*gm_error/frms
;
END
