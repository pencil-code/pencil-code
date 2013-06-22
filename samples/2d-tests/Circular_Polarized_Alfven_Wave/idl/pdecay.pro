pc_read_ts,o=ts
plot,ts.t,ts.uzpt,/nodata,yr=[-1,1]*.1
oplot,ts.t,ts.bzpt,col=122
oplot,ts.t,ts.uzpt
oplot,ts.t,ts.bzpt,col=122,li=2
;oplot,ts.t,ts.uzpt,col=55
END
