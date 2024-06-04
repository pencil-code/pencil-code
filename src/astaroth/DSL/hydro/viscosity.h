if (lvisc_nu_const) {  
   rhs += nu * (veclaplace(UU) + (1./3.) * gradient_of_divergence(UU) + 2. * stress_tensor(UU) * glnrho)
   + zeta * gradient_of_divergence(UU)
   }
if (lvisc_hyper3_nu_const) { 
   rhs += nu_hyper3 * (del6v(UU) + u_dot_grad_vec(gij5(UU),glnrho))
		   }
#if LSHOCK
if (lvisc_nu_shock) {
   rhs += nu_shock *  value(SHOCK) * ( divergence(UU) * glnrho + gradient_of_divergence(UU))
	}
#endif
