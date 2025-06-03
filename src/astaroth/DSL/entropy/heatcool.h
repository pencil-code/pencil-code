//checked 18.6.
interstellar_cool=0.0
for i in 0:ncool-1
{
    if (lncoolT[i] <= lnTT && lnTT < lncoolT[i+1]) {
      interstellar_cool = exp(lncoolH[i]+lnrho+lnTT*coolB[i])     // +=?
    }
}
heat = GammaUV*0.5*(1.0+tanh(cUV*(T0UV-TT)))
	
rhs += heat - interstellar_cool   // division by TT deferred to heat_ss.h
