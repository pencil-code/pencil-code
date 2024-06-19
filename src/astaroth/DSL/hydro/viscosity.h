//checked 17.6.

if (lvisc_nu_const){
   rhs +=  nu * ( veclaplace(UU) + (1./3.) * gradient_of_divergence(UU) + 2. * stress_tensor(UU) * glnrho )  // stress_tensor -> traceless strain tensor!
         + zeta * gradient_of_divergence(UU) * rho1
} // v otherwise
if (lvisc_hyper3_nu_const){
   rhs += nu_hyper3 * ( del6v(UU) + u_dot_grad_vec(gij5(UU), glnrho) )   //p%uij5glnrho   v
}
#if LSHOCK
if (lvisc_nu_shock){
   divu = divergence(UU)
   rhs += nu_shock*(value(SHOCK)*( divu * glnrho + gradient_of_divergence(UU) ) + divu*gradient(SHOCK))   // v
}
#endif
