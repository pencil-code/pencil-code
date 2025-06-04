//checked 17.6.

if (lvisc_nu_const){
   rhs +=  nu * ( laplace(UU) + (1./3.) * gradient_of_divergence(UU) + 2. * traceless_rateof_strain(UU) * glnrho )
         + zeta * gradient_of_divergence(UU) * rho1
}
if (lvisc_hyper3_nu_const){
   rhs += nu_hyper3 * ( del6(UU) + u_dot_grad_vec(gij5(UU), glnrho) )   //p%uij5glnrho
}
#if LSHOCK
if (lvisc_nu_shock){
   divu = divergence(UU)
   rhs += nu_shock*(SHOCK*( divu * glnrho + gradient_of_divergence(UU) ) + divu*gradient(SHOCK))
   reduce_max(nu_shock*SHOCK,AC_maxnu)
}
#endif
