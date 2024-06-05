#if LSHOCK
if (ldiff_shock) {
/*        if (ldensity_nolog) then
          call dot_mn(p%gshock,p%grho,tmp)
          fdiff = fdiff + diffrho_shock * (p%shock * p%del2rho + tmp)
        else
          if (ldiffusion_nolog) then
            call dot_mn(p%gshock,p%grho,tmp)
            fdiff = fdiff + p%rho1 * diffrho_shock * (p%shock * p%del2rho + tmp)
	  else
*/
            
            rhs = diffrho_shock * value(SHOCK) * (laplace(LNRHO) + norm2(glnrho) + dot(gradient(SHOCK),glnrho))
}
#endif
