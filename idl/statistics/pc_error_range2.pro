;$Id: pc_error_range.pro 9840 2008-09-05 07:29:37Z ajohan $
pro pc_error_range2,tt,a,mean=am,error=err,oplot=oplot,accum=accum,col=col
;
;  calculate averages for each stretch of time series that has the
;  length of one turnover time and use the standard deviation of
;  the resulting measurements as error. At the moment this requires
;  that time is already in units of turnover times.
;
;  Location:
;    pencil-code/idl/statistics/pc_error_range.pro
;
;  determine 3 ranges for which separate errors are calculated
;
nt=n_elements(tt)
it1=0
it4=nt-1
;
;  calculate how many points does one correlationt time correspond to
;
di=min(where(tt-tt(0) gt 1.))+1
nt2=fix((it4-it1+1)/di)*di
print,di,(it4-it1+1),nt2
;
;  reformat
;
ndi=nt2/di
a2=reform(a(it4-nt2+1:it4),di,ndi)
am_di=total(a2,1)/di
;
am=mean(am_di)
err=sqrt(mean((am_di-am)^2))
;
;  possibility of overplotting error range
;
if keyword_set(oplot) then begin
  default,col,122
  print,'DEVICE=',!d.name
  if !d.name eq 'X' then begin
    oplot,tt,(tt-tt+1.)*am,col=col
    if keyword_set(accum) then oplot,tt,pc_accumulate(a),col=55,li=2
    oplot,tt,(tt-tt+1.)*(am+err),col=188,li=2
    oplot,tt,(tt-tt+1.)*(am-err),col=188,li=2
  endif else begin
    ;oplot,tt,(tt-tt+1.)*am,li=3
    if keyword_set(accum) then oplot,tt,pc_accumulate(a),li=3
    oplot,tt,(tt-tt+1.)*(am+err),li=1
    oplot,tt,(tt-tt+1.)*(am-err),li=1
  endelse
endif
;
END
