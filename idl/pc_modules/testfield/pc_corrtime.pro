;$Id$
;
;  plot correlation function <u(t).u(t')>/<u^2>
;
;  In the run directory, need to set
;  in print.in
;    u2tm(e11.3)
;  in run.in
;    ltime_integrals=T  !(in run_pars)
;    luut_as_aux=T      !(in hydro_run_pars)
;  in start.in   
;      lwrite_aux=T     !(in init_pars)
;
;  see, e.g., pencil-code/axel/forced/alphaeta/OxJ32j_tau
;
pc_read_ts,o=ts
ratio=ts.u2tm/ts.urms^2
tt=ts.t
plot,tt,ratio
;
;  calculate correlation time
;
default,t1,50.
default,kf,5.
;spawn,'touch parameters.pro'
@parameters.pro
good=where(tt gt t1)
ratiom=mean(ratio(good))
ttgood=tt(good)
oplot,ttgood,ttgood*0.+ratiom,col=122
;
urmsm=mean(ts.urms(good))
St=urmsm*kf/ratiom
print,'ratiom,St=',ratiom,St
;
END
