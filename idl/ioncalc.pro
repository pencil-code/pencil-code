pro saha,yH,lnrho,ss,f,df
;
;  We want to find the root of f
;
  common constants,xHe,TT_ion,TT_ion_,lnrho_e,lnrho_H,lnrho_p, $
                   lnrho_He,lnrho_e_,ss_ion,kappa0,sigmaSB
    xHe=0.1
  ;
  ;  Only valid for one particular unit system:
  ;    unit_length=1e8
  ;    unit_velocity=1e6
  ;    unit_density=1e-7
  ;    unit_temperature=1
  ;
    TT_ion=157821.280907990     
    TT_ion_=8703.37946183766     
    lnrho_e=15.0796006637743     
    lnrho_H=26.3535589434417     
    lnrho_p=26.3527422402800     
    lnrho_He=28.4222860597670     
    lnrho_e_=10.7329728659836     
    ss_ion=5.90480511459068D-005
    kappa0=171073156.311650
    sigmaSB=5.6704d-16
  ;
  ;end common
  ;
  lnTT_=(2./3.)*((ss/ss_ion+(1-yH)*(alog(1-yH)-lnrho_H) $
                  +yH*(2*alog(yH)-lnrho_e-lnrho_p) $
                  +xHe*(alog(xHe)-lnrho_He))/(1+yH+xHe) $
                 +lnrho-2.5)
  f=lnrho_e-lnrho+1.5*lnTT_-exp(-lnTT_)+alog(1.-yH)-2.*alog(yH)
  dlnTT_=((2./3.)*(lnrho_H-lnrho_p-f-exp(-lnTT_))-1)/(1.+yH+xHe)
  df=dlnTT_*(1.5+exp(-lnTT_))-1./(1.-yH)-2./yH
end

function rtsafe,lnrho,ss
;
;  Safe Newton-Raphson root-finding algorithm
;
  common constants,xHe,TT_ion,TT_ion_,lnrho_e,lnrho_H,lnrho_p, $
                   lnrho_He,lnrho_e_,ss_ion,kappa0,sigmaSB
  maxit=1000
  yHacc=1e-5
  yHmax=1
  yHmin=0
  dyH=1
  dyHold=dyH
  yH=0.5
  saha,yH,lnrho,ss,f,df
  for i=1,maxit do begin
     if ((((yH-yHmin)*df-f)*((yH-yHmax)*df-f) gt 0) $
         or (abs(2*f) gt abs(dyHold*df))) then begin
        dyHold=dyH
        dyH=0.5*(yHmin-yHmax)
        yH=yHmax+dyH
        if (yHmax eq yH) then return,yH
     endif else begin
        dyHold=dyH
        dyH=f/df
        temp=yH
        yH=yH-dyH
        if (temp eq yH) then return,yH
     endelse
     if (abs(dyH) lt yHacc*yH) then return,yH
     saha,yH,lnrho,ss,f,df
     if (f<0) then yHmax=yH else yHmin=yH
  endfor
  print,'rtsafe: exceeded maximum iterations. maxit,f,yH=',maxit,f,yH
end

pro ioncalc,lnrho,ss,yH,TT
;
;  Calculate ionization degree and temperature
;
  common constants,xHe,TT_ion,TT_ion_,lnrho_e,lnrho_H,lnrho_p, $
                   lnrho_He,lnrho_e_,ss_ion,kappa0,sigmaSB

  yH=0*lnrho

  for n=0,n_elements(lnrho)-1 do yH[n]=rtsafe(lnrho[n],ss[n])

  lnTT_=(2./3.)*((ss/ss_ion+(1-yH)*(alog(1-yH)-lnrho_H) $
                  +yH*(2*alog(yH)-lnrho_e-lnrho_p) $
                  +xHe*(alog(xHe)-lnrho_He))/(1+yH+xHe) $
                 +lnrho-2.5)
  TT=exp(lnTT_)*TT_ion
end
