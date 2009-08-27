;$Id$
pro pc_fluct_stat,t,f,fm,frms,frms_error,plo=plo,good=good,debug=debug, $
  method=method
;
;  calculate average and rms fluctuations
;
fm=mean(f(good))
frms=sqrt(mean((f(good)-fm)^2))
if keyword_set(debug) then print,fm,frms
;
;  Choice of different methods:
;  Default is method 1 (splitting time series into 3 chunks
;    and determine the maximum departure)
;  Sometimes these errors are rather small and then it might be good
;    to split the time series into many shorter chunks of length 1
;    (assumed here that time is in correlation times) and determine
;    the width of the corresponding distribution.
;
if keyword_set(method) then begin
  imethod=method
endif else begin
  imethod=0
endelse
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
;
if imethod eq 0 then begin
  pc_error_range,t(good),f(good),mean=fm,error=fm_error,/accum,oplot=iplo
endif else if imethod eq 1 then begin
  pc_error_range2,t(good),f(good),mean=fm,error=fm_error,/accum,oplot=iplo
endif
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
if imethod eq 0 then begin
  pc_error_range,t(good),g,mean=gm,error=gm_error,/accum,oplot=iplo
endif else if imethod eq 1 then begin
  pc_error_range2,t(good),g,mean=gm,error=gm_error,/accum,oplot=iplo
endif
frms=gm^.5
frms_error=.5*gm_error/frms
;
END
