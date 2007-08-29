pro pc_errorrange,tt,a,mean=am,error=err,oplot=oplot
;
;  calculate error of array
;
default,col1,222
nt=n_elements(tt)
;
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
print,'xxx',am,aml,amu,err
;
if keyword_set(oplot) then begin
  oplot,tt,(tt-tt+1.)*am,col=122
  oplot,tt,(tt-tt+1.)*(am+err),col=188
  oplot,tt,(tt-tt+1.)*(am-err),col=188
;  plot,x,am,yr=[min(aml),max(amu)],_extra=_extra
;  polyfill,[x,reverse(x)],[aml,reverse(amu)],col=col1
;  oplot,x,am
endif
END
