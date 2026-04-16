;
pc_read_ts,obj=ts

!p.multi=[0,1,4]
window,xsize=600,ysize=700
!p.charsize=2.7

ymin=0.
ymax=11.
plot,ts.t,ts.apmax,ytit="ap",yr=[ymin,ymax]
oplot,ts.t,ts.apmin,col=122

ymax=0.
ymin=min(ts.vpzmin)
plot,ts.t,ts.vpzmin,ytit="vp",yr=[ymin,ymax]
oplot,ts.t,ts.vpzmax,col=122

ymin=0.
ymax=190.
plot,ts.t,ts.imskhm,ytit="mskhm",yr=[ymin,ymax]

np_sworm = ts.npar_found
np_phys  = ts.npar_found*ts.npswarmm
plot,ts.t,np_phys,ytit="np" ;,yr=[ymin,ymax]
oplot,ts.t,np_sworm,col=55


END
