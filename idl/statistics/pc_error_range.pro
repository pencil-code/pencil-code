pro pc_error_range,tt,a,mean=am,error=err,oplot=oplot
;
;  calculate averages for each third of time series and use
;  maximum departure from full average as error estimate
;
nt=n_elements(tt)
it1=0
it2=it1+(nt-it1)/3
it3=it2+(nt-it1)/3
it4=nt-1
;
am=total(a(it1:it4))/(it4-it1+1)
am1=total(a(it1:it2))/(it2-it1+1)
am2=total(a(it2:it3))/(it3-it2+1)
am3=total(a(it3:it4))/(it4-it3+1)
aml=am<am1<am2<am3
amu=am>am1>am2>am3
err=(amu-am) > (am-aml)
;
;  possibility of overplotting error range
;
if keyword_set(oplot) then begin
  oplot,tt,(tt-tt+1.)*am,col=122
  oplot,tt,(tt-tt+1.)*(am+err),col=188
  oplot,tt,(tt-tt+1.)*(am-err),col=188
endif
;
END
