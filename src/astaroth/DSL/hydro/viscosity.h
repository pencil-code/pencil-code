#if LVISC_NU_CONST {
rhs += nu * (veclaplace(UU) + (1./3.) * gradient_of_divergence(UU) + 2. * stress_tensor(UU) * gradient(LNRHO))
   + zeta * gradient_of_divergence(UU)}
#if LVISC_HYPER3_NU_CONST {
   rhs += nu_hyper3 * (del6s(UU) + u_dot_grad_vec(gij5(UU),gradient(LNRHO))
   }

