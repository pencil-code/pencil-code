;$Id: pc_error_range.pro 11565 2009-08-29 06:03:58Z AxelBrandenburg $
pro pc_error_range_array,tt,a,mean=am,error=err,oplot=oplot,accum=accum,col=col,li=li, ia=ia
;
;  Calculate averages for each third of time series and use
;  maximum departure from full average as error estimate.
;
;  Location:
;    pencil-code/idl/statistics/pc_error_range.pro
;
default,li,-1
default,col,122
;
;  determine 3 ranges for which separate errors are calculated
;
if not keyword_set(ia) then ia=0
nt=n_elements(tt[ia:*])
it1=ia
it2=it1+(nt-it1)/3
it3=it2+(nt-it1)/3
it4=nt-1+ia
;
;  calculate the 3 averages
;
am=total(a(it1:it4,*),1)/(it4-it1+1)
am1=total(a(it1:it2,*),1)/(it2-it1+1)
am2=total(a(it2:it3,*),1)/(it3-it2+1)
am3=total(a(it3:it4,*),1)/(it4-it3+1)
aml=am<am1<am2<am3
amu=am>am1>am2>am3
err=(amu-am) > (am-aml)
;print,'am,am1,am2,am3,err=',am,am1,am2,am3,err
;
END
