;
;  calculate cs2, TT1, and cp1tilde
;
  fff=lnrho_ion-llnrho+1.5*alog(TTT/TT_ion)-TT_ion/TTT+alog(1.-yyH)-2.*alog(yyH)
  dlnTT_dy=(lnmHmp-gamma1*(fff+TT_ion/TTT)-1.)/(1.+yyH+fHe)
  dfdy=dlnTT_dy*(1.5+TT_ion/TTT)-1./(1.-yyH)-2./yyH
  dlnTT_dlnrho=gamma1
  dfdlnrho=gamma1*TT_ion/TTT
  dydlnrho=-dfdlnrho/dfdy
  dlnPdlnrho=1.+dydlnrho/(1.+yyH+fHe)+dlnTT_dy*dydlnrho+dlnTT_dlnrho
  dlnTT_dss=gamma1/((1.+yyH+fHe)*ss_ion)
  dfdss=(1.+dfdlnrho)/((1.+yyH+fHe)*ss_ion)
  dydss=-dfdss/dfdy
  dlnPdss=dydss/(1.+yyH+fHe)+dlnTT_dy*dydss+dlnTT_dss
;
;  calculate sound speed, coefficient cp1tilde in
;  the expression (1/rho)*gradp = cs2*(gradlnrho + cp1tilde*gradss)
;  and internal energy for calculating thermal energy
;
  cs2=(1.+yyH+fHe)*ss_ion*TTT*dlnPdlnrho
  cp1tilde=dlnPdss/dlnPdlnrho
  eee=1.5*(1.+yyH+fHe)*ss_ion*TTT+yyH*ss_ion*TT_ion
  ppp=1.5*llnrho*ss_ion*TTT
