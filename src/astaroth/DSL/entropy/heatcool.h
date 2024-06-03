cool=0.0
for i in 0:ncool-1
{
    if (lncoolT[i] <= lnTT && lnTT < lncoolT[i+1]) {
      cool += exp(lncoolH[i]+lnrho+lnTT*coolB[i])
    }
}
heat=0.
//heat = GammaUV*0.5*(1.0+tanh(cUV*(T0UV-exp(lnTT))))
heatcool=heat-cool
