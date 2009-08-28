;  $Id$

if (not lionization and not lionization_fixed) then begin
  print,'Using simple equation of state...'

  if (par.lcalc_cp) then cp=k_B/(mu*m_H) else cp=1.  
  TT0=cs20/(cp * gamma_m1)
                                                    
  cs2=cs20*exp(gamma_m1*(llnrho-lnrho0)+gamma*sss)
  ppp=rho*cs2/gamma
  cp1tilde=1.
  eee=cs2/(gamma*gamma_m1)
  TTT=TT0*exp(gamma*sss+gamma_m1*(llnrho-lnrho0))
endif else begin
  xHe=par.xHe
  if (lionization_fixed) then begin 
    print,'Using fixed ionisation equation of state...'
    yH0=par.yH0
    xH2=par.xH2
    nabla_ad=(gamma-1.)/gamma
    yyH=reform(spread(spread(spread(yH0,0,nx),1,ny),2,nz))
    lnTT=lnTTss*ss+lnTTlnrho*lnrho+lnTT0
    TT=exp(lnTT)
    TTT=reform(TT(l1:l2,m1:m2,n1:n2))
    cs2=gamma*(1.+yH0+xHe-xH2)*ss_ion*TTT
    cp1tilde=nabla_ad/(1.+yH0+xHe-xH2)/ss_ion
    ee=1.5*(1.+yH0+xHe-xH2)*ss_ion*TTT+yH0*ss_ion*TT_ion
    pp=(1.+yH0+xHe-xH2)*exp(llnrho)*TTT*ss_ion
  end else begin
    print,'Using full ionisation equation of state...'
    if (iyH ne 0 and ilnTT ne 0) then begin
      yyH=reform(yH(l1:l2,m1:m2,n1:n2))
      TTT=reform(exp(lnTT(l1:l2,m1:m2,n1:n2)))
      ;
      ;  calculate cs2, TT1, and cp1tilde
      ;
      fff=lnrho_e-llnrho+1.5*alog(TTT/TT_ion)-TT_ion/TTT+alog(1.-yyH)-2.*alog(yyH)
      dlnTT_dy=((2./3.)*(lnrho_H-lnrho_p-fff-TT_ion/TTT)-1)/(1.+yyH+xHe)
      dfdy=dlnTT_dy*(1.5+TT_ion/TTT)-1./(1.-yyH)-2./yyH
      dlnTT_dlnrho=(2./3.)
      dfdlnrho=(2./3.)*TT_ion/TTT
      dydlnrho=-dfdlnrho/dfdy
      dlnPdlnrho=1.+dydlnrho/(1.+yyH+xHe)+dlnTT_dy*dydlnrho+dlnTT_dlnrho
      dlnTT_dss=(2./3.)/((1.+yyH+xHe)*ss_ion)
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
    end else begin
      print,"Errr... not implemented calculation of ionization fraction in IDL (yet)"
      stop
    endelse
  endelse  
endelse
