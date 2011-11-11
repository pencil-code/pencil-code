;$Id$
pro pc_error_range,tt,a,mean=am,error=err,oplot=oplot,accum=accum, $
  col=col,li=li,sqr=sqr
;
;  Calculate averages for each third of time series and use
;  maximum departure from full average as error estimate.
;  With /sqr keyword take square root (used for rms value).
;  In that case, the squared values are to be supplied.
;
;  Location:
;    pencil-code/idl/statistics/pc_error_range.pro
;
default,li,-1
default,col,122
;
;  determine 3 ranges for which separate errors are calculated
;
nt=n_elements(tt)
it1=0
it2=it1+(nt-it1)/3
it3=it2+(nt-it1)/3
it4=nt-1
;
;  calculate the 3 averages
;
am=total(a(it1:it4))/(it4-it1+1)
am1=total(a(it1:it2))/(it2-it1+1)
am2=total(a(it2:it3))/(it3-it2+1)
am3=total(a(it3:it4))/(it4-it3+1)
aml=am<am1<am2<am3
amu=am>am1>am2>am3
err=(amu-am) > (am-aml)
print,'am,am1,am2,am3,err=',am,am1,am2,am3,err
;
;  Assume that t is already normalized in correlation times, i.e. by urms*kf.
;  The following method is commented out, because it wasn't too useful.
;
;t_interval=tt(it4)-tt(it1)
;err=sqrt(total(a(it1:it4)^2)/((it4-it1+1)*t_interval))
;
;  possibility of overplotting error range
;
if keyword_set(oplot) then begin
  ;print,'DEVICE=',!d.name
  if !d.name eq 'X' then begin
    oplot,tt,(tt-tt+1.)*am,col=col,li=li
    if li eq -1 then begin
      if keyword_set(accum) then oplot,tt,pc_accumulate(a),col=55,li=2
      oplot,tt,(tt-tt+1.)*(am+err),col=188,li=2
      oplot,tt,(tt-tt+1.)*(am-err),col=188,li=2
    endif
  endif else begin
    oplot,tt,(tt-tt+1.)*am,li=li
    if li eq -1 then begin
      if keyword_set(accum) then oplot,tt,pc_accumulate(a),li=3
      oplot,tt,(tt-tt+1.)*(am+err),li=1
      oplot,tt,(tt-tt+1.)*(am-err),li=1
    endif
  endelse
endif
;
;  Take square root with /sqr keyword (used for rms value).
;  In that case, the squared values are to be supplied.
;  Prevent taking square roots of negative values and return 0 instead.
;
if keyword_set(sqr) then begin
  am=sqrt(am >  0.)
  err=.5*err/am
endif
;
END
