// checked 17.6.
if (ldiff_hyper3lnrho){
      rhs += diffrho_hyper3*del6(LNRHO)
}

#if LSHOCK
if (ldiff_shock){

 del2rho = laplace(RHO)         // \nabla^2 RHO  or  \nabla^2 LNRHO !
 gshock = gradient(SHOCK)

 if (ldensity_nolog){
      rhs += diffrho_shock * (value(SHOCK) * del2rho + dot(gshock,glnrho))
 }
 else
 {
//     if (ldiffusion_nolog){ 
//              call dot_mn(p%gshock,p%grho,tmp)
//              rhs =  exp(-value(LNRHO)) * diffrho_shock * (value(SHOCK) * laplace(exp(LNRHO)) + dot(gradient(SHOCK),glnrho))
//  	}
      rhs += diffrho_shock * (value(SHOCK) * (del2rho + norm2(glnrho)) + dot(gshock,glnrho))
 }
}
#endif
