;  $Id: thermodynamics.pro,v 1.5 2003-08-09 19:47:27 mee Exp $

if (cs20 ne impossible) then begin
  print,'Using simple equation of state...'
  cs2=cs20*exp(gamma1*llnrho+gamma*sss)
  ppp=rho*cs2/gamma
  cp1tilde=1.
  eee=cs2/gamma1
end else begin
  if (iyH ne 0) then begin 
    print,'Using full ionisation equation of state...'
    yyH=yH(l1:l2,m1:m2,n1:n2)
  end else begin
    print,'Using fixed ionisation equation of state...'
    seedarray=fltarr(nx)
    seedarray[*]=yH0
    yyH=spread(spread(seedarray,1,ny),2,nz)
  endelse  
  if (iTT ne 0) then begin
    TTT=TT(l1:l2,m1:m2,n1:n2)
  end else begin
    yyH_term=yyH*(2.*log(yyH)-lnrho_e-lnrho_p)
    yH_term[where (yyH eq 0.)]=0.
    one_yyH_term=(1.-yyH)*(log(1.-yyH)-lnrho_H)
    one_yyH_term[where (yyH eq 1.)]=0.
    xHe_term=xHe*(log(xHe)-lnrho_He)
    xHe_term[where (xHe eq 0.)]=0.

    lnTTT_=(2./3.)*((sss/ss_ion $
                   + one_yyH_term $
                   + yyH_term $
                   + xHe_term)/(1.+yyH+xHe) $
                   + lnrho-2.5)
    TTT=exp(lnTTT_)*TT_ion
  endelse
  ;
  ;  calculate cs2, TT1, and cp1tilde
  ;
  fff=lnrho_e-llnrho+1.5*alog(TTT/TT_ion)-TT_ion/TTT+alog(1.-yyH)-2.*alog(yyH)
  dlnTT_dy=(lnmHmp-gamma1*(fff+TT_ion/TTT)-1.)/(1.+yyH+xHe)
  dfdy=dlnTT_dy*(1.5+TT_ion/TTT)-1./(1.-yyH)-2./yyH
  dlnTT_dlnrho=gamma1
  dfdlnrho=gamma1*TT_ion/TTT
  dydlnrho=-dfdlnrho/dfdy
  dlnPdlnrho=1.+dydlnrho/(1.+yyH+xHe)+dlnTT_dy*dydlnrho+dlnTT_dlnrho
  dlnTT_dss=gamma1/((1.+yyH+xHe)*ss_ion)
  dfdss=(1.+dfdlnrho)/((1.+yyH+xHe)*ss_ion)
  dydss=-dfdss/dfdy
  dlnPdss=dydss/(1.+yyH+xHe)+dlnTT_dy*dydss+dlnTT_dss
  ;
  ;  calculate sound speed, coefficient cp1tilde in
  ;  the expression (1/rho)*gradp = cs2*(gradlnrho + cp1tilde*gradss)
  ;  and internal energy for calculating thermal energy
  ;
  cs2=(1.+yyH+xHe)*ss_ion*TTT*dlnPdlnrho
  cp1tilde=dlnPdss/dlnPdlnrho
  eee=1.5*(1.+yyH+xHe)*ss_ion*TTT+yyH*ss_ion*TT_ion
  ppp=(1.+yyH+xHe)*exp(llnrho)*ss_ion*TTT
endelse
