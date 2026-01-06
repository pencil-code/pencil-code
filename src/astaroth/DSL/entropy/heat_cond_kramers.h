// Kramers opacity-based heat conduction.
    if (lheatc_kramers){
      //cv1 = 1./cv
      rho1 = exp(-value(LNRHO))    // v
      lnTT = lnTT0+cv1*value(SS)+(gamma-1.)*(value(LNRHO)-lnrho0)  // v
      glnrho = gradient(LNRHO)     // v
    
      glnTT  = cv1*gradient(SS) + (gamma-1.)*glnrho   //  v
      del2lnTT = cv1*laplace(SS) + (gamma-1.)*laplace(LNRHO)  // v
    
      Krho1 = hcond0_kramers * pow(rho1,(2.*nkramers+1.)) * pow(exp(lnTT),(6.5*nkramers))   // = K/rho   v

      g2 = dot(-2.*nkramers*glnrho+(6.5*nkramers+1.)*glnTT,glnTT)   // v
      rhs += Krho1*(del2lnTT+g2)    // v

      if (lchit_total && lgravz && chi_t != 0.) {
        g2 = dot(glnrho+glnTT,gss)
        rhs += chi_t*(chit_prof_stored[vertexIdx.z-NGHOST]*(del2ss+g2) + gss.z*dchit_prof_stored[vertexIdx.z-NGHOST])

        if (step_num == 0 && ldt && lcourant_dt)
        {
       	   //reduce_max(Krho1/cv1+chi_t*chit_prof_stored[vertexIdx.z-NGHOST],maxchi)
       	   chitot += Krho1/cv1+chi_t*chit_prof_stored[vertexIdx.z-NGHOST]
        }
      } else {
        if (step_num == 0 && ldt && lcourant_dt)
        {
           chitot += Krho1/cv1
        }
      }
    }
