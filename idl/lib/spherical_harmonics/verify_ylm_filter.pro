nphi=16*4
ntheta=16*4*4
phi=2*!pi*findgen(nphi)/nphi
theta0=0*!dtor
;
theta=grange(theta0,!pi-theta0,ntheta)
;
dtheta=(!pi-theta0*2.)/ntheta
theta=dtheta*(findgen(ntheta)+.5)
;
;  swap the two lines
;
bb_slice=sin(theta)#sin(phi)
bb_slice=cos(theta)#replicate(1.,n_elements(phi))
;
nl=11 & mmax=10
nl=2 & mmax=2
ylm_filter,bb_slice,theta,phi,nl,mmax,Etot,Em,Pm,modes=m,/quiet,brmm_all=brmm_all
;
iphi=5
plot,theta,bb_slice(*,iphi),yr=[-1,1]*1.14
oplot,theta,brmm_all(*,iphi),col=122,li=2
END
