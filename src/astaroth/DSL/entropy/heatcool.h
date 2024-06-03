cool=0.0

#define tanh 


for i in 0:AC_ncool-1{
    if (AC_lncoolT[i] <= lnTT && lnTT < AC_lncoolT[i+1]){
	    cool=cool+exp(AC_lncoolH[i]+lnrho+lnTT*AC_coolB[i])
    }
}
heat = AC_GammaUV*0.5*(1.0+!tanh(AC_cUV*(AC_T0UV-exp(lnTT))))
heatcool=heat-cool
