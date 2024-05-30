heat_conduction(int step_num) {

// Kramers opacity-based heat conduction.

      cv1 = 1./cv    
      rho1 = exp(-value(LNRHO))    // v
      lnTT = lnTT0+cv1*value(SS)+(gamma-1.)*(value(LNRHO)-lnrho0)  // v
      glnrho = gradient(LNRHO)     // v
    
      glnTT  = cv1*gradient(SS) + (gamma-1.)*glnrho   //  v
      del2lnTT = cv1*laplace(SS) + (gamma-1.)*laplace(LNRHO)  // v
    
      Krho1 = hcond0_kramers * pow(rho1,(2.*nkramers+1.)) * pow(exp(lnTT),(6.5*nkramers))   // = K/rho   v

      reduce_max(step_num==0,Krho1/cv1,AC_maxchi)

      g2 = dot(-2.*nkramers*glnrho+(6.5*nkramers+1.)*glnTT,glnTT)   // v
      return Krho1*(del2lnTT+g2)    // v
}
