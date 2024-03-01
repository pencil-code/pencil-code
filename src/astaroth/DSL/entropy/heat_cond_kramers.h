heat_conduction_kramers() {

      cv1 = 1.0/AC_cv    
      rho1 = exp(-value(VTXBUF_LNRHO))    // v
      lnTT = AC_lnTT0+cv1*value(VTXBUF_ENTROPY)+(AC_gamma-1.0)*(value(VTXBUF_LNRHO)-AC_lnrho0)  // v
      glnrho = gradient(VTXBUF_LNRHO)     // v
    
      glnTT  = cv1*gradient(VTXBUF_ENTROPY) + (AC_gamma-1.0)*glnrho   //  v
      del2lnTT = cv1*laplace(VTXBUF_ENTROPY) + (AC_gamma-1.0)*laplace(VTXBUF_LNRHO)  // v
    
      Krho1 = AC_hcond0_kramers * pow(rho1,(2.0*AC_nkramers+1.0)) * pow(exp(lnTT),(6.5*AC_nkramers))   // = K/rho   v
    
      g2=dot(-2.0*AC_nkramers*glnrho+(6.5*AC_nkramers+1.)*glnTT,glnTT)   // v
      return Krho1*(del2lnTT+g2)    // v
}
