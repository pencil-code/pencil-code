;$Id: calc_bl.pro,v 1.4 2017/08/09 23:25:57 brandenb Exp $
;
pro ylm_filter,bb_slice,theta,phi,nl,mmax,Etot,Em,Pm, $
  debug=debug,modes=modes,quiet=quiet
;
;  compute the spherical harmonics decomposition
;  input: bb_slice,theta,phi,nl,mmax
;    where bb_slice(itheta,iphi) is a surface slice at one time
;  output: Etot,Em,Pm,/debug
;
ntheta=n_elements(theta)
nphi=n_elements(phi)
;
Brmm_all=fltarr(ntheta,nphi)
;
dtheta=theta[1]-theta[0]
sintheta=sin(theta)
;
dphi=phi[1]-phi[0]
phi=phi-dphi/2.
;
nl=60
nl=120
mmax=28
Em=fltarr(mmax)
Emodd=fltarr(mmax)
Emeven=fltarr(mmax)
modes=findgen(mmax)
for m=0,mmax-1 do begin
;
;  Fourier transform
;
br1c=total(bb_slice*spread(cos(m*phi),[0],[ntheta]),2)/nphi/2.
;br1s=total(bb_slice*spread(sin(m*phi),[0],[ntheta]),2)/nphi/2.
;
Etot=0.
brm=complexarr(ntheta,nphi)
for itheta=0,ntheta-1 do begin
  brm[itheta,*]=fft(bb_slice[itheta,*])
  Etot=Etot+.5*dphi*dtheta*total(sintheta[itheta]*bb_slice[itheta,*]^2)
endfor
;
brl=complexarr(nl,2*nl+1)
;
for l=m,nl-1 do begin
  ;fact=.25*(2*l+1)*factorial(double(l-m))/factorial(double(l+m))
  fact=(2*l+1)*factorial(double(l-m))/factorial(double(l+m))
  brl(l,m)=dtheta*total(sintheta*brm[*,m]*plm(double(theta),double(l),abs(m)))*fact
endfor
if keyword_set(debug) then print,'m=',m,fact,min(plm(double(theta),l-1,abs(m))),max(plm(double(theta),double(l-1),abs(m)))
;
Brmm=fltarr(ntheta,nphi)
Beven=fltarr(ntheta,nphi)
Bodd=fltarr(ntheta,nphi)
for l=m,nl-1 do begin
  brlphi      =cos(m*phi)*float(brl(l,m ))-sin(m*phi)*imaginary(brl(l,m))
  plmfactor=plm(double(theta),double(l),abs(m))
  dBrmm=plmfactor#brlphi
  Brmm=Brmm+dBrmm
  if l mod 2 then Beven=Beven+dBrmm else Bodd=Bodd+dBrmm
endfor
;
for itheta=0,ntheta-1 do begin
  Em[m]=Em[m]+.5*dphi*dtheta*total(sintheta[itheta]*Brmm[itheta,*]^2)
  Emodd[m]=Emodd[m]+.5*dphi*dtheta*total(sintheta[itheta]*Bodd[itheta,*]^2)
  Emeven[m]=Emeven[m]+.5*dphi*dtheta*total(sintheta[itheta]*Beven[itheta,*]^2)
endfor
;
;  sum of all contributions
;
Brmm_all=Brmm_all+Brmm
;
if keyword_set(debug) then begin
  !p.multi=[0,1,3]
  lev=grange(-1,1,30)*.025
  contour,brmm_all,/fill,lev=lev
  contour,bb_slice,/fill,lev=lev
  iphi=0
  plot,theta,Brmm_all[*,iphi],yr=[-1,1]*.02
  oplot,theta,bb_slice[*,iphi],li=1
  oplot,theta,brmm[*,iphi],col=122,li=1,thick=4
  oplot,theta,br1c,col=55,li=2,thick=2
  wait,.04
  !p.multi=0
endif
;
if not keyword_set(quiet) then print,'mode m out of mmax: ',m,mmax
endfor
;
Pm=(Emeven-Emodd)/(Emeven+Emodd)
END
