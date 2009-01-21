; $Id$
;
;  horizontal plane; not as interesting as xz cross-section;
;  for which one uses pu_xz.pro
;
n=nz/2+3
;
theta=grange(0.,2*!pi,51)
costheta=cos(theta)
sintheta=sin(theta)
contour,uu(*,*,n,2),x,y,/fil,nlev=30
vel_a,uu(*,*,n,0),uu(*,*,n,1),x,y,nvec=1000,/over,len=.03
oplot,par2.r_int*costheta,par2.r_int*sintheta,thick=5
oplot,par2.r_ext*costheta,par2.r_ext*sintheta,thick=5
oplot,par2.r_int*costheta,par2.r_int*sintheta,thick=3,col=0,li=2
oplot,par2.r_ext*costheta,par2.r_ext*sintheta,thick=3,col=0,li=2
;
END
