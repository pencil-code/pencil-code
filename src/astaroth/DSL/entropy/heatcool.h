cool=0.0
for i in range(ncool){
    if (lncoolT[i] <= lnTT && lnTT < lncoolT[i+1])
    {cool=cool+exp(lncoolH[i]+lnrho+lnTT*coolB[i])}
}
heat = GammaUV*0.5*(1.0+tanh(cUV*(T0UV-exp(lnTT))))
heatcool=heat-cool
rhs += heatcool
