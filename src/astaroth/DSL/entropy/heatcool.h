//checked 18.6.
interstellar_cool=0.0
for i in 0:ncool-1
{
    if (lncoolt[i] <= lnTT && lnTT < lncoolt[i+1]) {
      interstellar_cool = exp(lncoolh[i]+lnrho+lnTT*coolb[i])     // +=?
    }
}
heat = gammauv*0.5*(1.0+tanh(cuv*(t0uv-TT)))
	
rhs += heat - interstellar_cool   // division by TT deferred to heat_ss.h
