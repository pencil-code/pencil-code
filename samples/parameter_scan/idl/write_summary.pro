;$Id$
;
;  read
;
pc_read_ts,o=ts
pc_read_param,o=param,/param2
@parameters
;
;  determine good range
;
tend=max(ts.t)
good=where(ts.t gt .9*tend)
uxbmm=mean(ts.uxbm(good))
plot,ts.t,ts.uxbm,yst=3
oplot,ts.t(good),ts.t(good)*0.+uxbmm,col=122,thick=6,li=2
;
fo='(e8.1,e11.3,2x,a)'
B0=param.B_ext(0)
;
;  file should be defined via parameters
;
openw,1,file
print,B0,uxbmm,run,fo=fo
printf,1,B0,uxbmm,run,fo=fo
close,1
;
;  append file to the same file on the parent directory
;
spawn,'cat '+file+'>>../'+file
END
