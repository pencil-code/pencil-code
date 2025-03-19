//checked 18.6.
cool=0.0
for i in 0:ncool-1
{
    if (lncoolT[i] <= lnTT && lnTT < lncoolT[i+1]) {
      cool = exp(lncoolH[i]+lnrho+lnTT*coolB[i])     // +=?
    }
}
heat = GammaUV*0.5*(1.0+tanh(cUV*(T0UV-TT)))
	
rhs += heat - cool   // division by TT deferred to heat_ss.h
