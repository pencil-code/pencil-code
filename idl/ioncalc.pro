;$Id$
pro saha,yH,lnrho,ss,f,df,par
;
;  We want to find the root of f
;
  @data/pc_constants.pro
  xHe=par.xHe

  lnTT_=(2./3.)*((ss/ss_ion+(1-yH)*(alog(1-yH)-lnrho_H) $
                  +yH*(2*alog(yH)-lnrho_e-lnrho_p) $
                  +xHe*(alog(xHe)-lnrho_He))/(1+yH+xHe) $
                 +lnrho-2.5)
  f=lnrho_e-lnrho+1.5*lnTT_-exp(-lnTT_)+alog(1.-yH)-2.*alog(yH)
  dlnTT_=((2./3.)*(lnrho_H-lnrho_p-f-exp(-lnTT_))-1)/(1.+yH+xHe)
  df=dlnTT_*(1.5+exp(-lnTT_))-1./(1.-yH)-2./yH

end
;*******************************************************************
function rtsafe,lnrho,ss,par
;
;  Safe Newton-Raphson root-finding algorithm
;
  @data/pc_constants.pro
  yHacc=par.yHacc

  maxit=1000
  yHmax=1
  yHmin=0
  dyH=1
  dyHold=dyH
  yH=0.5

  saha,yH,lnrho,ss,f,df,par

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
     saha,yH,lnrho,ss,f,df,par
     if (f<0) then yHmax=yH else yHmin=yH
  endfor

  print,'rtsafe: exceeded maximum iterations. maxit,f,yH=',maxit,f,yH

end
;*******************************************************************
pro ioncalc,lnrho,ss,par,yH,TT
;
;  Calculate ionization degree and temperature
;
  @data/pc_constants.pro
  xHe=par.xHe

  size=size(lnrho)
  dim=size[0]

  if dim eq 0 then yH=rtsafe(lnrho,ss,par)

  if dim eq 1 then begin
    yH=fltarr(size[1])
    for l=0,size[1]-1 do begin
      yH[l]=rtsafe(lnrho[l],ss[l],par)
    endfor
  endif
  if dim eq 2 then begin
    yH=fltarr(size[1],size[2])
    for l=0,size[1]-1 do begin
    for m=0,size[2]-1 do begin
      yH[l,m]=rtsafe(lnrho[l,m],ss[l,m],par)
    endfor
    endfor
  endif
  if dim eq 3 then begin
    yH=fltarr(size[1],size[2],size[3])
    for l=0,size[1]-1 do begin
    for m=0,size[2]-1 do begin
    for n=0,size[3]-1 do begin
      yH[l,m,n]=rtsafe(lnrho[l,m,n],ss[l,m,n],par)
    endfor
    endfor
    endfor
  endif

  lnTT_=(2./3.)*((ss/ss_ion+(1-yH)*(alog(1-yH)-lnrho_H) $
                  +yH*(2*alog(yH)-lnrho_e-lnrho_p) $
                  +xHe*(alog(xHe)-lnrho_He))/(1+yH+xHe) $
                 +lnrho-2.5)
  TT=exp(lnTT_)*TT_ion

end
