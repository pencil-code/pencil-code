#if LSHOCK
if (ldiff_shock) {
if (ldensity_nolog){
          rhs = diffrho_shock * (value(SHOCK) * laplace(RHO) + dot(gradient(SHOCK),glnrho))
}
else
{
          /*if (ldiffusion_nolog){ 
            call dot_mn(p%gshock,p%grho,tmp)
            rhs =  exp(-value(LNRHO)) * diffrho_shock * (value(SHOCK) * laplace(exp(LNRHO)) + dot(gradient(SHOCK),glnrho))
	}else*/
            rhs = diffrho_shock * (value(SHOCK) * (laplace(LNRHO) + norm2(glnrho)) + dot(gradient(SHOCK),glnrho))

}
}
#endif
