power,'_np','_np',k=k,spec1=spec1,spec2=spec2,i=n,tt=t,/noplot
;
default,w,.1
default,yr,[1,1e3]
yr=[1,1e4]
;
nt=n_elements(t)
for it=0,nt-1 do begin
  plot_oo,k,spec1[*,it],xr=[1,max(k)],yr=yr
  wait,w
endfor
;
END
