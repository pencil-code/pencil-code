pc_read_pvar,obj=fp,/proc_id

tmp=fp.iproc

i=where(tmp ne -1,n) & if (n ne 0) then tmp=tmp[i]

plot,histogram(tmp,nbin=max(tmp)),ps=10,xs=3,ys=3,charsize=2.,$
     title='Load balance',xtitle='proc #',ytitle='Npar'
   
end
